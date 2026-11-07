"""Build a new, traceable ToxCast CSV snapshot from the pinned EPA release."""
import argparse
import csv
import hashlib
import io
import json
import re
import shutil
import sqlite3
import zipfile
from collections import Counter
from collections import defaultdict
from contextlib import ExitStack, contextmanager
from decimal import Decimal, InvalidOperation
from itertools import zip_longest
from pathlib import Path
from tempfile import TemporaryDirectory
from uuid import NAMESPACE_URL, uuid5

from openpyxl import load_workbook
from download_toxcast import DEFAULT_LOCK, RELEASE, fingerprint, load_lock, verify

POLICY_VERSION = 1
MISSING = {"", "NA", "N/A", "#N/A", "NAN", "NULL", "NONE"}
SUMMARY_FIELDS = ("toxcast_mc_result_count", "toxcast_sc_result_count", "toxcast_endpoint_count")
LINEAGE = ("chemical_id", "source_dataset", "source_release", "source_file", "source_sheet",
           "source_row_number", "source_sha256", "evidence_kind", "result_id", "activity_status", "potency_unit")
SID_PATTERN = re.compile(r"DTXSID[0-9]{7,}")


def clean(value):
    text = "" if value is None else str(value).strip()
    return "" if text.upper() in MISSING else text


def identifier(value):
    text = clean(value)
    if not text:
        return ""
    try:
        number = Decimal(text)
        if number.is_finite() and number == number.to_integral() and number >= 0:
            return str(int(number))
    except InvalidOperation:
        pass
    raise ValueError(f"Invalid integer identifier: {text!r}")


def activity_status(row, screen):
    raw = clean(row.get("hitc", row.get("hitcall")))
    if not raw:
        return "missing"
    try:
        value = Decimal(raw)
    except InvalidOperation:
        return "invalid"
    if not value.is_finite():
        return "invalid"
    if screen == "mc" and str(row.get("modl", "")).strip().lower() in {"none", "nofit"}:
        return "unfit"
    if value < 0:
        return "negative_direction" if screen == "mc" else "unfit"
    if value > 1 or (screen == "sc" and value not in (0, 1)):
        return "invalid"
    return "active" if value >= (Decimal("0.90") if screen == "mc" else 1) else "inactive"


def potency_unit(row, screen):
    # Keep the actual concentration unit. The release also contains CF (concentration factor).
    # Missing units stay unknown; no values are converted or silently assumed to be uM.
    if screen != "mc":
        return ""
    units = clean(row.get("conc_unit"))
    return units


def cytotox_value_status(row):
    try:
        value = Decimal(clean(row.get("cytotox_median_um")))
        tested = Decimal(clean(row.get("ntested")))
        hits = Decimal(clean(row.get("nhit")))
    except InvalidOperation:
        return "unknown"
    weak_support = tested <= 0 or hits < 2 or hits / tested < Decimal("0.05")
    return "default_1000_um" if value == 1000 and weak_support else "estimated_from_assays"


@contextmanager
def csv_rows(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream, strict=True)
        fields = tuple(reader.fieldnames or ())
        validate_header(fields)
        def rows():
            for number, row in enumerate(reader, 2):
                if None in row or any(v is None for v in row.values()):
                    raise ValueError(f"Malformed CSV {path}, logical row {number}")
                yield number, row
        yield fields, rows()


def validate_header(fields):
    if not fields or any(not x for x in fields) or len(set(fields)) != len(fields):
        raise ValueError("Missing or duplicate table headers")


def init_db(path):
    db = sqlite3.connect(path)
    db.executescript("""
        PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF; PRAGMA temp_store=FILE;
        CREATE TABLE chemicals (id TEXT PRIMARY KEY, baseline INTEGER, tox21 INTEGER);
        CREATE TABLE sids (sid TEXT, id TEXT, PRIMARY KEY(sid,id));
        CREATE TABLE new_chemicals (id TEXT PRIMARY KEY, sid TEXT UNIQUE);
        CREATE TABLE results (screen TEXT, rid TEXT, chemical_id TEXT, aeid TEXT, spid TEXT,
          status TEXT, signature TEXT, ac50 TEXT, flagged INTEGER, unit TEXT, PRIMARY KEY(screen,rid));
        CREATE TABLE assays (aeid TEXT PRIMARY KEY, name TEXT, assay_source TEXT, usability TEXT, documented INTEGER);
        CREATE TABLE targets (aeid TEXT, kind TEXT, target_id TEXT, name TEXT,
          PRIMARY KEY(aeid,kind,target_id));
    """)
    return db


def load_baseline(db, base, review):
    manifest = json.loads((base / "dsstox_import.json").read_text())
    audit = json.loads((review / "key_review_manifest.json").read_text())
    for name, digest in audit["input_sha256"].items():
        if fingerprint(base / name) != digest:
            raise ValueError(f"Baseline fingerprint mismatch: {name}")
    for name, digest in audit["output_sha256"].items():
        if fingerprint(review / name) != digest:
            raise ValueError(f"Key-review fingerprint mismatch: {name}")
    if audit["dsstox_release_sha256"] != manifest["source_sha256"]:
        raise ValueError("Key review belongs to a different DSSTox release")
    with csv_rows(review / "key_updates.csv") as (_, rows):
        updates = {row["chemical_id"]: row for _, row in rows}
    changed = set()
    with csv_rows(base / "master_chemicals.csv") as (original_fields, original):
        with csv_rows(review / "master_chemicals_reviewed.csv") as (fields, reviewed):
            if fields != original_fields or set(SUMMARY_FIELDS) & set(fields):
                raise ValueError("Reviewed master has incompatible or already-imported columns")
            batch = []
            for left, right in zip_longest(original, reviewed):
                if left is None or right is None:
                    raise ValueError("Review changed the master row count")
                _, old = left
                _, row = right
                differences = {k for k in fields if old[k] != row[k]}
                if differences:
                    if differences - {"standard_inchikey", "identity_status"} or row["chemical_id"] not in updates:
                        raise ValueError("Unapproved modification in reviewed master")
                    changed.add(row["chemical_id"])
                cid = row["chemical_id"]
                if not cid:
                    raise ValueError("Empty master chemical ID")
                has_tox21 = any(row[k] != "" for k in fields if k.startswith(("NR-", "SR-")))
                batch.append((cid, 1, int(has_tox21)))
                if len(batch) == 10000:
                    db.executemany("INSERT INTO chemicals VALUES (?,?,?)", batch)
                    batch.clear()
            db.executemany("INSERT INTO chemicals VALUES (?,?,?)", batch)
    if len(changed) != audit["derived_keys_applied"] or changed != set(updates):
        raise ValueError("Reviewed-key changes do not reconcile")
    with csv_rows(base / "dsstox_records.csv") as (_, rows):
        batch = []
        source_count = 0
        for _, row in rows:
            if row["source_sha256"] != manifest["source_sha256"]:
                raise ValueError("DSSTox source lineage mismatch")
            batch.append((clean(row["DTXSID"]), row["chemical_id"]))
            source_count += 1
            if len(batch) == 10000:
                db.executemany("INSERT OR IGNORE INTO sids VALUES (?,?)", batch)
                batch.clear()
        db.executemany("INSERT OR IGNORE INTO sids VALUES (?,?)", batch)
    if db.execute("SELECT COUNT(*) FROM sids LEFT JOIN chemicals ON sids.id=chemicals.id WHERE chemicals.id IS NULL").fetchone()[0]:
        raise ValueError("DSSTox crosswalk contains dangling chemical IDs")
    if source_count != manifest["linked_rows"] or db.execute("SELECT COUNT(*) FROM chemicals").fetchone()[0] != audit["master_rows"]:
        raise ValueError("Baseline counts do not reconcile")
    db.commit()
    print(f"Verified baseline: {audit['master_rows']:,} chemicals, {len(changed)} reviewed keys", flush=True)
    return fields, audit


class Importer:
    def __init__(self, db, staging, source_hashes):
        self.db, self.staging, self.source_hashes = db, staging, source_hashes
        self.stats = Counter()
        self.tables = []
        self.identity_cache = {}
        self.chemical_records = set()
        self.stack = ExitStack()
        self.chemical_writer = self.writer("toxcast_chemicals.csv", (
            "chemical_id", "match_status", "dsstox_substance_id", "chid", "chnm", "casn",
            "source_file", "source_sheet", "source_row_number", "source_sha256"))
        self.review_writer = self.writer("toxcast_review.csv", (
            "reason", "chemical_id", "dsstox_substance_id", "chid", "aeid", "spid",
            "source_file", "source_sheet", "source_row_number", "source_sha256"))

    def writer(self, name, fields):
        path = self.staging / name
        path.parent.mkdir(parents=True, exist_ok=True)
        stream = self.stack.enter_context(path.open("w", newline="", encoding="utf-8"))
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        return writer

    def resolve(self, row, provenance):
        sid = clean(row.get("dsstox_substance_id"))
        if sid not in self.identity_cache:
            if not SID_PATTERN.fullmatch(sid):
                result = ("", "missing_or_invalid_dtxsid")
            else:
                ids = [r[0] for r in self.db.execute("SELECT id FROM sids WHERE sid=?", (sid,))]
                if len(ids) == 1:
                    result = (ids[0], "matched_dtxsid")
                elif len(ids) > 1:
                    result = ("", "ambiguous_dtxsid")
                else:
                    cid = str(uuid5(NAMESPACE_URL, f"researchchem:dsstox:substance:{sid}"))
                    self.db.execute("INSERT INTO chemicals VALUES (?,0,0)", (cid,))
                    self.db.execute("INSERT INTO new_chemicals VALUES (?,?)", (cid, sid))
                    result = (cid, "new_substance")
            self.identity_cache[sid] = result
        cid, status = self.identity_cache[sid]
        record = (cid, status, sid, clean(row.get("chid")), clean(row.get("chnm")), clean(row.get("casn")))
        if record not in self.chemical_records:
            self.chemical_records.add(record)
            self.chemical_writer.writerow(dict(zip(("chemical_id", "match_status", "dsstox_substance_id", "chid", "chnm", "casn"), record)) | provenance)
            self.stats["chemical_crosswalk_rows"] += 1
        return cid, status

    def review(self, reason, row, cid, provenance):
        self.review_writer.writerow({"reason": reason, "chemical_id": cid,
            **{k: clean(row.get(k)) for k in ("dsstox_substance_id", "chid", "aeid", "spid")}, **provenance})
        self.stats["review_rows"] += 1

    def process(self, name, sheet, fields, rows, output_name, kind, archive_hash):
        validate_header(fields)
        if set(fields) & set(LINEAGE):
            raise ValueError(f"Reserved lineage field in source {name}/{sheet}")
        writer = self.writer(output_name, (*fields, *LINEAGE, *(("cytotox_value_status",) if kind == "cytotox" else ())))
        screen = kind if kind in ("mc", "sc") else ""
        if screen and not {"aeid", "spid", "dsstox_substance_id"}.issubset(fields):
            raise ValueError(f"Missing observation identifiers: {name}/{sheet}")
        is_result = screen and (screen == "mc" or "s2id" in fields)
        if is_result and not ({"hitc", "hitcall"} & set(fields)):
            raise ValueError("Missing activity field")
        source_hash = self.source_hashes.get(name, archive_hash)
        count, duplicate_count = 0, 0
        for number, row in rows:
            if set(row) != set(fields):
                raise ValueError(f"Malformed source table: {name}/{sheet}:{number}")
            provenance = dict(source_file=name, source_sheet=sheet, source_row_number=number, source_sha256=source_hash)
            cid, match = ("", "")
            chemical_row = row
            if "dtxsid" in fields:
                chemical_row = {**row, "dsstox_substance_id": row["dtxsid"],
                                "chnm": row.get("name", ""), "casn": row.get("casrn", "")}
            if "dsstox_substance_id" in chemical_row:
                cid, match = self.resolve(chemical_row, provenance)
                if not cid:
                    self.review(match, chemical_row, cid, provenance)
            status, rid, unit = "", "", ""
            if is_result:
                native_id = "m5id" if screen == "mc" else "s2id"
                rid = identifier(row.get(native_id))
                if not rid:
                    raise ValueError(f"Missing {native_id}: {name}/{sheet}:{number}")
                aeid = identifier(row["aeid"])
                if not aeid:
                    raise ValueError("Result has no assay identifier")
                if not clean(row["spid"]):
                    self.review("missing_sample_id", row, cid, provenance)
                status = activity_status(row, screen)
                unit = potency_unit(row, screen)
                if status == "invalid":
                    self.review("invalid_hitcall", row, cid, provenance)
                signature_values = [cid, aeid, clean(row["spid"]), status,
                    clean(row.get("hitc", row.get("hitcall")))]
                if screen == "mc":
                    signature_values.extend(clean(row.get(k)) for k in ("ac50", "acc", "bmd", "fitc", "modl", "mc6_flags"))
                    signature_values.append(unit)
                signature = json.dumps(signature_values, separators=(",", ":"))
                flagged = bool(clean(row.get("mc6_flags")))
                inserted = self.db.execute("INSERT OR IGNORE INTO results VALUES (?,?,?,?,?,?,?,?,?,?)", (
                    screen, rid, cid or None, aeid, clean(row["spid"]), status, signature, clean(row.get("ac50")), int(flagged), unit)).rowcount
                if not inserted:
                    old = self.db.execute("SELECT signature FROM results WHERE screen=? AND rid=?", (screen, rid)).fetchone()[0]
                    if old != signature:
                        raise ValueError(f"Conflicting values for {screen} processing ID {rid}")
                    duplicate_count += 1
            if kind == "assays" and "aeid" in fields:
                self.db.execute("INSERT INTO assays VALUES (?,?,?,?,?)", (
                    identifier(row["aeid"]), clean(row.get("assay_component_endpoint_name")),
                    clean(row.get("assay_source_name")), clean(row.get("data_usability")), 1))
            if kind == "targets" and clean(row.get("target_id")) and clean(row.get("target_type")):
                self.db.execute("INSERT OR IGNORE INTO targets VALUES (?,?,?,?)", (
                    identifier(row["aeid"]), clean(row.get("target_type")), identifier(row.get("target_id")),
                    clean(row.get("official_symbol")) or clean(row.get("official_full_name"))))
            evidence = "processed_experiment" if screen else "derived_model" if kind == "derived" else "cytotoxicity_summary_with_defaults" if kind == "cytotox" else "annotation_or_quality"
            extra = {"cytotox_value_status": cytotox_value_status(row)} if kind == "cytotox" else {}
            writer.writerow({**row, **provenance, "chemical_id": cid, "source_dataset": "ToxCast",
                "source_release": RELEASE, "evidence_kind": evidence, "result_id": rid,
                "activity_status": status, "potency_unit": unit, **extra})
            count += 1
            if count % 250000 == 0:
                self.db.commit()
                print(f"  {output_name}: {count:,} source rows", flush=True)
        self.db.commit()
        self.tables.append(dict(output=output_name, source_file=name, source_sheet=sheet, kind=kind,
                                rows=count, duplicate_result_rows=duplicate_count, source_sha256=source_hash))
        print(f"Exported {output_name}: {count:,} rows", flush=True)

    def close(self):
        self.stack.close()


def workbook_tables(stream):
    book = load_workbook(stream, read_only=True, data_only=True)
    try:
        for sheet in book:
            sheet.reset_dimensions()
            iterator = sheet.iter_rows(values_only=True)
            first = next(iterator, None)
            if first is None:
                continue
            fields = tuple(str(v) if v is not None else "" for v in first)
            validate_header(fields)
            def rows(iterator=iterator, fields=fields):
                for number, values in enumerate(iterator, 2):
                    if all(v is None for v in values):
                        continue
                    if len(values) > len(fields) and any(v is not None for v in values[len(fields):]):
                        raise ValueError("Workbook row extends beyond its header")
                    values = tuple(values[:len(fields)]) + (None,) * max(0, len(fields) - len(values))
                    yield number, {k: "" if v is None else str(v) for k, v in zip(fields, values)}
            yield sheet.title, fields, rows()
    finally:
        book.close()


def slug(text):
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def import_sources(importer, source_dir, lock):
    hashes = importer.source_hashes
    archive_hash = hashes["INVITRODB_SUMMARY.zip"]
    with zipfile.ZipFile(source_dir / "INVITRODB_SUMMARY.zip") as archive:
        inventory = [{"member": i.filename, "bytes": i.file_size, "crc32": f"{i.CRC:08x}"} for i in archive.infolist()]
        # The independently provided AUG2025 files are byte-identical to AUG2024-named members.
        for name in hashes:
            if name.endswith(".xlsx"):
                member = name.replace("AUG2025", "AUG2024")
                with archive.open(member) as stream:
                    if hashlib.file_digest(stream, "sha256").hexdigest() != hashes[name]:
                        raise ValueError(f"Standalone metadata differs from archive member: {name}")
        for prefix, kind, output in (
            ("assay_annotations", "assays", "toxcast_assays.csv"),
            ("assay_target_mappings", "targets", "toxcast_targets.csv"),
            ("cytotox", "cytotox", "toxcast_cytotoxicity.csv"),
            ("analytical_qc", "quality", "toxcast_analytical_qc.csv"),
        ):
            names = [n for n in hashes if n.startswith(prefix) and n.endswith(".xlsx")]
            if len(names) != 1:
                raise ValueError(f"Expected one pinned {prefix} workbook")
            name = names[0]
            with (source_dir / name).open("rb") as stream:
                for sheet, fields, rows in workbook_tables(stream):
                    combined = kind == "assays" and sheet == "annotations_combined"
                    actual_kind = kind if combined or kind != "assays" else "assay_detail"
                    filename = output if combined or kind != "assays" else f"assays/{slug(sheet)}.csv"
                    importer.process(name, sheet, fields, rows, filename, actual_kind, archive_hash)
        name = "mc5-6_winning_model_fits-flags_invitrodbv4_3_AUG2024.csv"
        with archive.open(name) as binary, io.TextIOWrapper(binary, encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream, strict=True)
            fields = tuple(reader.fieldnames or ())
            def records():
                for number, row in enumerate(reader, 2):
                    if None in row or any(v is None for v in row.values()):
                        raise ValueError(f"Malformed MC row {number}")
                    yield number, row
            importer.process(name, "", fields, records(), "toxcast_mc_observations.csv", "mc", archive_hash)
        selected = {
            "sc1_sc2_invitrodb_v4_3_AUG2024.xlsx": "sc",
            "methods_invitrodb_v4_3_AUG2024.xlsx": "methods",
            "assay_quality_invitrodb_v4_3_AUG2024.xlsx": "quality",
            "mc7_invitrodb_v4_3_AUG2024.xlsx": "derived",
            **{n: "derived" for n in archive.namelist() if n.startswith("endocrine_models/") and n.endswith(".xlsx")},
        }
        for name, kind in selected.items():
            # Copy only this member into temporary staging: avoids repeated seeks through deflated ZIP data.
            with TemporaryDirectory(prefix=".xlsx-", dir=importer.staging.parent) as temp:
                path = Path(temp) / "source.xlsx"
                with archive.open(name) as src, path.open("wb") as dst:
                    shutil.copyfileobj(src, dst, length=1024**2)
                for sheet, fields, rows in workbook_tables(path):
                    if kind == "sc":
                        filename = "toxcast_sc_observations.csv" if sheet == "sc2" else "toxcast_sc_replicates.csv"
                    elif kind == "derived":
                        filename = f"derived/{slug(Path(name).stem)}_{slug(sheet)}.csv"
                    else:
                        filename = f"{kind}/{slug(sheet)}.csv"
                    importer.process(name, sheet, fields, rows, filename, kind, archive_hash)
        return inventory


def count(db, sql, args=()):
    return db.execute(sql, args).fetchone()[0]


def reconcile_endpoint_catalog(db, staging):
    """Retain historical results whose endpoint annotation is absent from this release."""
    missing = {r[0] for r in db.execute("""SELECT aeid FROM results UNION SELECT aeid FROM targets
      EXCEPT SELECT aeid FROM assays""")}
    examples = {}
    names = defaultdict(set)
    if missing:
        for screen in ("mc", "sc"):
            with csv_rows(staging / f"toxcast_{screen}_observations.csv") as (_, rows):
                for _, row in rows:
                    aeid = identifier(row["aeid"])
                    if aeid in missing:
                        names[aeid].add(clean(row.get("aenm")))
                        examples.setdefault(aeid, {k: row[k] for k in ("source_file", "source_sheet", "source_row_number", "source_sha256")})
        db.executemany("INSERT INTO assays VALUES (?,?,?,?,0)",
            [(aeid, next(iter(names[aeid])) if len(names[aeid]) == 1 else "", "", "") for aeid in sorted(missing)])
    with (staging / "toxcast_assay_review.csv").open("w", newline="", encoding="utf-8") as stream:
        fields = ("aeid", "reason", "source_endpoint_names", "processed_results", "source_file", "source_sheet", "source_row_number", "source_sha256")
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for aeid in sorted(missing, key=int):
            writer.writerow({"aeid": aeid, "reason": "missing_released_assay_annotation",
                "source_endpoint_names": json.dumps(sorted(names[aeid])),
                "processed_results": count(db, "SELECT COUNT(*) FROM results WHERE aeid=?", (aeid,)),
                **examples.get(aeid, {})})
    with (staging / "toxcast_endpoints.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("aeid", "assay_name", "annotation_status"))
        writer.writerows((a, n, "documented" if documented else "missing_released_annotation")
            for a,n,documented in db.execute("SELECT aeid,name,documented FROM assays ORDER BY CAST(aeid AS INTEGER)"))
    if missing:
        print(f"Retained {len(missing)} referenced endpoints missing released annotations; see toxcast_assay_review.csv", flush=True)


def build_summaries(db, staging, reviewed_master, fields):
    db.executescript("""
      CREATE INDEX results_chemical ON results(chemical_id);
      CREATE INDEX results_endpoint ON results(aeid);
      CREATE TABLE coverage AS SELECT chemical_id,
        SUM(screen='mc') AS mc, SUM(screen='sc') AS sc, COUNT(DISTINCT aeid) AS endpoints,
        MAX(status IN ('active','inactive')) AS interpretable,
        MAX(screen='mc' AND status IN ('active','inactive')) AS mc_interpretable,
        MAX(screen='sc' AND status IN ('active','inactive')) AS sc_interpretable
        FROM results WHERE chemical_id IS NOT NULL GROUP BY chemical_id;
      CREATE UNIQUE INDEX coverage_id ON coverage(chemical_id);
    """)
    reconcile_endpoint_catalog(db, staging)
    if count(db, "SELECT COUNT(*) FROM results r LEFT JOIN assays a ON r.aeid=a.aeid WHERE a.aeid IS NULL"):
        raise ValueError("Observation refers to an undocumented assay endpoint")
    if count(db, "SELECT COUNT(*) FROM targets t LEFT JOIN assays a ON t.aeid=a.aeid WHERE a.aeid IS NULL"):
        raise ValueError("Target mapping refers to an undocumented assay endpoint")
    if count(db, "SELECT COUNT(*) FROM results r LEFT JOIN chemicals c ON r.chemical_id=c.id WHERE r.chemical_id IS NOT NULL AND c.id IS NULL"):
        raise ValueError("Observation refers to a missing master chemical")
    with (staging / "master_chemicals.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=(*fields, *SUMMARY_FIELDS))
        writer.writeheader()
        summaries = {r[0]: r[1:] for r in db.execute("SELECT chemical_id,mc,sc,endpoints FROM coverage")}
        with csv_rows(reviewed_master) as (_, rows):
            for _, row in rows:
                writer.writerow({**row, **dict(zip(SUMMARY_FIELDS, summaries.get(row["chemical_id"], (0, 0, 0))))})
        for cid, sid in db.execute("SELECT id,sid FROM new_chemicals ORDER BY sid,id"):
            row = dict.fromkeys(fields, "")
            row.update(chemical_id=cid, identity_status="toxcast_substance_only")
            writer.writerow({**row, **dict(zip(SUMMARY_FIELDS, summaries.get(cid, (0, 0, 0))))})
    # Independent grouping by result identities must reconcile all master summary totals.
    if sum(r[0] for r in db.execute("SELECT mc FROM coverage")) != count(db, "SELECT COUNT(*) FROM results WHERE screen='mc' AND chemical_id IS NOT NULL"):
        raise ValueError("MC coverage reconciliation failed")
    if sum(r[0] for r in db.execute("SELECT sc FROM coverage")) != count(db, "SELECT COUNT(*) FROM results WHERE screen='sc' AND chemical_id IS NOT NULL"):
        raise ValueError("SC coverage reconciliation failed")
    metrics = {
        "baseline_chemicals": count(db, "SELECT COUNT(*) FROM chemicals WHERE baseline=1"),
        "baseline_with_tox21": count(db, "SELECT COUNT(*) FROM chemicals WHERE tox21=1"),
        "expanded_chemicals": count(db, "SELECT COUNT(*) FROM chemicals"),
        "new_chemical_entries": count(db, "SELECT COUNT(*) FROM new_chemicals"),
        "toxcast_source_result_chemicals": count(db, "SELECT COUNT(*) FROM coverage"),
        "toxcast_interpretable_chemicals": count(db, "SELECT COUNT(*) FROM coverage WHERE interpretable=1"),
        "baseline_newly_interpretable": count(db, "SELECT COUNT(*) FROM chemicals c JOIN coverage v ON c.id=v.chemical_id WHERE c.baseline=1 AND c.tox21=0 AND v.interpretable=1"),
        "baseline_new_source_presence": count(db, "SELECT COUNT(*) FROM chemicals c JOIN coverage v ON c.id=v.chemical_id WHERE c.baseline=1 AND c.tox21=0"),
        "tox21_chemicals_with_toxcast_results": count(db, "SELECT COUNT(*) FROM chemicals c JOIN coverage v ON c.id=v.chemical_id WHERE c.tox21=1"),
        "tox21_chemicals_with_interpretable_toxcast": count(db, "SELECT COUNT(*) FROM chemicals c JOIN coverage v ON c.id=v.chemical_id WHERE c.tox21=1 AND v.interpretable=1"),
        "new_entries_with_interpretable_results": count(db, "SELECT COUNT(*) FROM chemicals c JOIN coverage v ON c.id=v.chemical_id WHERE c.baseline=0 AND v.interpretable=1"),
        "processed_results": count(db, "SELECT COUNT(*) FROM results"),
        "unresolved_processed_results": count(db, "SELECT COUNT(*) FROM results WHERE chemical_id IS NULL"),
        "mc_results": count(db, "SELECT COUNT(*) FROM results WHERE screen='mc'"),
        "sc_results": count(db, "SELECT COUNT(*) FROM results WHERE screen='sc'"),
        "mc_flagged_results": count(db, "SELECT COUNT(*) FROM results WHERE screen='mc' AND flagged=1"),
        "assay_endpoints_documented": count(db, "SELECT COUNT(*) FROM assays WHERE documented=1"),
        "assay_endpoints_missing_annotations": count(db, "SELECT COUNT(*) FROM assays WHERE documented=0"),
        "processed_results_missing_assay_annotations": count(db, "SELECT COUNT(*) FROM results r JOIN assays a ON r.aeid=a.aeid WHERE a.documented=0"),
        "assay_endpoints_with_results": count(db, "SELECT COUNT(DISTINCT aeid) FROM results"),
        "mapped_targets": count(db, "SELECT COUNT(*) FROM (SELECT DISTINCT kind,target_id FROM targets)"),
        "assay_target_links": count(db, "SELECT COUNT(*) FROM targets"),
        "results_with_target_annotation": count(db, "SELECT COUNT(*) FROM results r WHERE EXISTS(SELECT 1 FROM targets t WHERE r.aeid=t.aeid)"),
        "sample_disagreement_groups": count(db, "SELECT COUNT(*) FROM (SELECT chemical_id,aeid,screen FROM results WHERE chemical_id IS NOT NULL AND status IN ('active','inactive') GROUP BY chemical_id,aeid,screen HAVING COUNT(DISTINCT status)>1)"),
    }
    for screen, status, amount in db.execute("SELECT screen,status,COUNT(*) FROM results GROUP BY screen,status"):
        metrics[f"{screen}_{status}_results"] = amount
    metrics["fixed_catalog_covered_after"] = metrics["baseline_with_tox21"] + metrics["baseline_newly_interpretable"]
    metrics["expanded_catalog_covered_after"] = metrics["fixed_catalog_covered_after"] + metrics["new_entries_with_interpretable_results"]
    with (staging / "coverage_summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("metric", "value"))
        writer.writerows(sorted(metrics.items()))
    with (staging / "endpoint_coverage.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("aeid", "source_results", "resolved_chemicals", "interpretable_chemicals", "active_chemicals", "target_links"))
        writer.writerows(db.execute("""SELECT r.aeid,COUNT(*),COUNT(DISTINCT r.chemical_id),
          COUNT(DISTINCT CASE WHEN r.status IN ('active','inactive') THEN r.chemical_id END),
          COUNT(DISTINCT CASE WHEN r.status='active' THEN r.chemical_id END),
          (SELECT COUNT(*) FROM targets t WHERE t.aeid=r.aeid)
          FROM results r GROUP BY r.aeid ORDER BY CAST(r.aeid AS INTEGER)"""))
    with (staging / "sample_disagreements.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("chemical_id", "aeid", "screening_type", "active_results", "inactive_results", "samples"))
        writer.writerows(db.execute("""SELECT chemical_id,aeid,screen,SUM(status='active'),SUM(status='inactive'),COUNT(DISTINCT spid)
          FROM results WHERE chemical_id IS NOT NULL AND status IN ('active','inactive')
          GROUP BY chemical_id,aeid,screen HAVING COUNT(DISTINCT status)>1 ORDER BY chemical_id,aeid,screen"""))
    example = db.execute("""SELECT r.rid,r.chemical_id,r.aeid,a.name FROM results r JOIN assays a ON r.aeid=a.aeid
      WHERE r.screen='mc' AND r.chemical_id IS NOT NULL AND r.status='active'
      AND r.ac50 NOT IN ('','NA','NaN') AND r.unit='uM' AND EXISTS(SELECT 1 FROM targets t WHERE t.aeid=r.aeid)
      ORDER BY CAST(r.rid AS INTEGER) LIMIT 1""").fetchone()
    if example:
        with csv_rows(staging / "toxcast_mc_observations.csv") as (_, rows):
            evidence = next(row for _, row in rows if row["result_id"] == example[0])
        targets = list(db.execute("SELECT kind,target_id,name FROM targets WHERE aeid=? ORDER BY kind,target_id", (example[2],)))
        (staging / "trace_example.json").write_text(json.dumps({"observation": evidence, "assay": {"aeid": example[2], "name": example[3]},
            "target_annotations": [dict(zip(("target_type", "target_id", "name"), t)) for t in targets]}, indent=2, sort_keys=True) + "\n")
    write_meeting_summary(staging, metrics, example)
    return metrics


def write_meeting_summary(staging, m, example):
    before = 100 * m["baseline_with_tox21"] / m["baseline_chemicals"]
    fixed = 100 * m["fixed_catalog_covered_after"] / m["baseline_chemicals"]
    expanded = 100 * m["expanded_catalog_covered_after"] / m["expanded_chemicals"]
    gap_closed = 100 * m["baseline_newly_interpretable"] / (m["baseline_chemicals"] - m["baseline_with_tox21"])
    text = f"""# ToxCast v4.3 integration: meeting summary

The verified baseline contains {m['baseline_chemicals']:,} chemical entries and {m['baseline_with_tox21']:,} chemicals with Tox21 results.

| Measure | Result |
| --- | ---: |
| Baseline coverage | {before:.4f}% |
| Original catalog covered after import | {m['fixed_catalog_covered_after']:,} ({fixed:.4f}%) |
| Original chemicals newly gaining interpretable results | {m['baseline_newly_interpretable']:,} |
| Original experimental-coverage gap closed | {gap_closed:.4f}% |
| Expanded catalog | {m['expanded_chemicals']:,} |
| Expanded catalog covered | {m['expanded_catalog_covered_after']:,} ({expanded:.4f}%) |
| New substance entries | {m['new_chemical_entries']:,} |
| ToxCast chemicals with source results | {m['toxcast_source_result_chemicals']:,} |
| ToxCast chemicals with interpretable results | {m['toxcast_interpretable_chemicals']:,} |
| Existing Tox21 chemicals gaining ToxCast detail | {m['tox21_chemicals_with_toxcast_results']:,} |
| Distinct MC / SC processed results | {m['mc_results']:,} / {m['sc_results']:,} |
| Unresolved processed results retained | {m['unresolved_processed_results']:,} |
| Assay endpoints with results | {m['assay_endpoints_with_results']:,} |
| Referenced endpoints missing released annotations | {m['assay_endpoints_missing_annotations']:,} |
| Distinct annotated targets (gene/AOP/key event) | {m['mapped_targets']:,} |
| Chemical/endpoint/screen groups with sample disagreement | {m['sample_disagreement_groups']:,} |

Coverage means at least one active or inactive processed assay result linked to a chemical. Missing, unfit, negative-direction, and invalid calls do not count as interpretable coverage. Source presence includes these records. These figures do not establish overall human toxicity or safety.

ToxCast includes Tox21-origin assays. Chemical overlap is reported; overlapping labels are not counted as independent experimental evidence. The original 12 Tox21 columns and individual observations remain intact. Different assay IDs are not automatically equivalent biological endpoints.

All caution flags, replicate rows, and analytical QC are retained. Assay-to-target links are annotations, not proof of a causal path to an organ outcome. Endocrine/QSAR scores, administered-equivalent doses, and cytotoxicity summary defaults are kept separate from experimental coverage. A 1000 uM cytotoxicity default is not a measured safe concentration.

Some result endpoints lack rows in the released assay annotations. These results remain in the observation tables and coverage counts; `toxcast_endpoints.csv` marks the annotation gap and `toxcast_assay_review.csv` gives source examples. No biological description or target is inferred for an undocumented endpoint.

The three master summary columns count distinct MC results, distinct SC results, and distinct assay endpoints, respectively. Replicate rows and one-to-many target annotations do not inflate these counts. New substance entries do not receive inferred structures.

"""
    if example:
        text += f"## Traceable example\n\nChemical `{example[1]}` → MC processing result `{example[0]}` → assay `{example[2]}` ({example[3]}) → the target annotations in `trace_example.json`. That file retains the original potency, units, flags, sample ID, source logical row, and release fingerprint.\n\n"
    text += "## Next research milestone\n\nReview the coverage gain and quality limitations with the professor, then agree on ChEMBL ingestion and the first DAG retrieval benchmark. No ChEMBL data or agentic retrieval engine is included in this snapshot.\n"
    (staging / "MEETING_SUMMARY.md").write_text(text)


def run(base, review, source_dir, output, lock_file=DEFAULT_LOCK):
    csv.field_size_limit(16 * 1024**2)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    lock = load_lock(lock_file)
    print("Verifying pinned ToxCast sources", flush=True)
    for item in lock["files"]:
        verify(source_dir / item["filename"], item)
    hashes = {x["filename"]: x["sha256"] for x in lock["files"]}
    input_paths = [*sorted(p for p in base.iterdir() if p.is_file()), *sorted(p for p in review.iterdir() if p.is_file()), lock_file]
    input_hashes = {str(p.resolve()): fingerprint(p) for p in input_paths}
    with TemporaryDirectory(prefix=".toxcast-", dir=output.parent) as temp:
        work = Path(temp)
        staging = work / "snapshot"
        staging.mkdir()
        db = init_db(work / "index.sqlite")
        importer = Importer(db, staging, hashes)
        try:
            fields, audit = load_baseline(db, base, review)
            inventory = import_sources(importer, source_dir, lock)
            importer.close()
            metrics = build_summaries(db, staging, review / "master_chemicals_reviewed.csv", fields)
            for path in sorted(base.iterdir()):
                if path.is_file() and path.name != "master_chemicals.csv":
                    shutil.copyfile(path, staging / path.name)
                    if fingerprint(staging / path.name) != input_hashes[str(path.resolve())]:
                        raise ValueError("Baseline companion copy mismatch")
            (staging / "key_review").mkdir()
            for path in sorted(review.iterdir()):
                if path.is_file() and path.name != "master_chemicals_reviewed.csv":
                    shutil.copyfile(path, staging / "key_review" / path.name)
            for item in lock["files"]:
                verify(source_dir / item["filename"], item)
            for path in input_paths:
                if fingerprint(path) != input_hashes[str(path.resolve())]:
                    raise ValueError(f"Input changed during import: {path}")
            # Reconcile physical output counts independently of importer counters.
            for table in importer.tables:
                with csv_rows(staging / table["output"]) as (_, rows):
                    actual = sum(1 for _ in rows)
                if actual != table["rows"]:
                    raise ValueError(f"Output row reconciliation failed: {table['output']}")
            expected_results = sum(t["rows"] - t["duplicate_result_rows"] for t in importer.tables
                if t["kind"] == "mc" or t["output"] == "toxcast_sc_observations.csv")
            if expected_results != metrics["processed_results"]:
                raise ValueError("Source/result accounting failed")
            output_hashes = {str(p.relative_to(staging)): fingerprint(p) for p in sorted(staging.rglob("*")) if p.is_file()}
            manifest = {"release": RELEASE, "policy_version": POLICY_VERSION,
                "input_sha256": input_hashes, "source_files": lock["files"], "archive_inventory": inventory,
                "tables": importer.tables, "metrics": metrics, "review_counts": dict(importer.stats),
                "derived_keys_preserved": audit["derived_keys_applied"], "output_sha256": output_hashes,
                "checks": {"source_results_accounted_for": True, "master_ids_preserved": True,
                  "original_master_values_preserved": True, "assay_references_valid": True,
                  "target_references_valid": True, "baseline_companions_unchanged": True,
                  "inputs_unchanged": True, "master_counts_reconciled": True}}
            (staging / "toxcast_import.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
            if output.exists() or output.is_symlink():
                raise FileExistsError(output)
            staging.rename(output)
        finally:
            importer.close()
            db.close()
    print(json.dumps(metrics, indent=2, sort_keys=True), flush=True)
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-snapshot", required=True, type=Path)
    parser.add_argument("--key-review", required=True, type=Path)
    parser.add_argument("--source-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--source-lock", default=DEFAULT_LOCK, type=Path)
    args = parser.parse_args()
    run(args.base_snapshot, args.key_review, args.source_dir, args.output, args.source_lock)
