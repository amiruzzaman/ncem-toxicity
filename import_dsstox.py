"""Build a new chemical master snapshot from a full DSSTox CCD export.

From the repo root, run:
    python3 import_dsstox.py tox21_initial DSSToxCCDdump.csv tox21_dsstox_complete
"""
import argparse
import csv
import hashlib
import json
import re
import shutil
import sqlite3
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter
from uuid import NAMESPACE_URL, uuid5

MASTER = "master_chemicals.csv"
OBSERVATIONS = "tox21_observations.csv"
MANIFEST = "dsstox_import.json"
POLICY_VERSION = 2
SOURCE_FIELDS = (
    "DTXSID", "PREFERRED_NAME", "CASRN", "DTXCID", "INCHIKEY", "IUPAC_NAME",
    "SMILES", "MOLECULAR_FORMULA", "AVERAGE_MASS", "MONOISOTOPIC_MASS",
    "QSAR_READY_SMILES", "MS_READY_SMILES", "IDENTIFIER",
)
METADATA = (
    "PREFERRED_NAME", "CASRN", "DTXSID", "DTXCID", "MOLECULAR_FORMULA",
    "AVERAGE_MASS", "MONOISOTOPIC_MASS",
)
STAGED_FIELDS = (*METADATA, "SMILES")
KEY_FORMAT = re.compile(r"[A-Z]{14}-[A-Z]{8}[SN]A-[A-Z]")
MISSING = {"", "NA", "N/A", "NULL", "NAN", "NONE"}


def fingerprint(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def clean(value):
    value = value.strip()
    return None if value.upper() in MISSING else value


def read_master(path):
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, strict=True)
        fields = tuple(reader.fieldnames or ())
        required = {"chemical_id", "standard_inchikey", "representative_smiles", "identity_status"}
        if not required.issubset(fields) or len(fields) != len(set(fields)):
            raise ValueError("Invalid master header")
        rows = list(reader)
    ids = set()
    by_key = {}
    for row in rows:
        if None in row or None in row.values():
            raise ValueError("Malformed master row")
        key, chemical_id = row["standard_inchikey"], row["chemical_id"]
        if not chemical_id or chemical_id in ids or (key and key in by_key):
            raise ValueError("Missing or duplicate master identity")
        ids.add(chemical_id)
        if key:
            by_key[key] = row
    return fields, rows, by_key, ids


def source_rows(stream):
    reader = csv.DictReader(stream, strict=True)
    if tuple(reader.fieldnames or ()) != SOURCE_FIELDS:
        raise ValueError("Unexpected DSSTox header")
    for number, row in enumerate(reader, start=2):
        if None in row or None in row.values():
            raise ValueError(f"Malformed DSSTox record {number}")
        key = clean(row["INCHIKEY"])
        yield number, row, key if key and KEY_FORMAT.fullmatch(key) else None


def index_substances(db, source):
    # Resolve unkeyed records even when their keyed counterpart occurs later.
    db.execute("CREATE TABLE substances (sid TEXT PRIMARY KEY, low TEXT, high TEXT)")
    with source.open(newline="", encoding="utf-8-sig") as stream:
        pairs = ((sid, key, key) for _, row, key in source_rows(stream)
                 if key and (sid := clean(row["DTXSID"])))
        db.executemany("""
            INSERT INTO substances VALUES (?, ?, ?)
            ON CONFLICT(sid) DO UPDATE SET
                low = MIN(low, excluded.low), high = MAX(high, excluded.high)
        """, pairs)


def stage_source(db, source, digest, directory, by_key, ids, stats):
    index_substances(db, source)
    columns = ", ".join(f"{field} TEXT" for field in STAGED_FIELDS)
    db.execute(f"CREATE TABLE records (chemical_id TEXT, key TEXT, {columns})")
    placeholders = ", ".join("?" for _ in range(2 + len(STAGED_FIELDS)))
    with (
        source.open(newline="", encoding="utf-8-sig") as stream,
        (directory / "dsstox_records.csv").open("w", newline="", encoding="utf-8") as linked,
        (directory / "dsstox_review.csv").open("w", newline="", encoding="utf-8") as review,
    ):
        linked_writer, review_writer = csv.writer(linked), csv.writer(review)
        fields = ("chemical_id", "source_row_number", "source_sha256", *SOURCE_FIELDS)
        linked_writer.writerow(fields)
        review_writer.writerow(("review_reason", *fields))

        def linked_rows():
            # Source row numbers count logical CSV records, with the header at 1.
            for number, row, key in source_rows(stream):
                stats["input_rows"] += 1
                reason = sid = None
                if key is None:
                    reason = "malformed_key" if clean(row["INCHIKEY"]) else "missing_key"
                    stats[reason] += 1
                    sid = clean(row["DTXSID"])
                    match = db.execute(
                        "SELECT low, high FROM substances WHERE sid = ?", (sid,)
                    ).fetchone()
                    if match and match[0] == match[1]:
                        key = match[0]
                        stats["resolved_by_dtxsid_rows"] += 1

                existing = by_key.get(key)
                if existing is not None:
                    chemical_id = existing["chemical_id"]
                else:
                    if key is not None:
                        identity = f"researchchem:chemical:{key}"
                    elif sid:
                        identity = f"researchchem:dsstox:substance:{sid}"
                    else:
                        identity = f"researchchem:dsstox:record:{digest}:{number}"
                    chemical_id = str(uuid5(NAMESPACE_URL, identity))
                    if chemical_id in ids:
                        raise ValueError(f"Generated chemical ID is already in use for record {number}")

                record = (chemical_id, number, digest, *(row[field] for field in SOURCE_FIELDS))
                stats["linked_rows"] += 1
                linked_writer.writerow(record)
                if reason is not None:
                    review_writer.writerow((reason, *record))
                yield (chemical_id, key, *(clean(row[field]) for field in STAGED_FIELDS))

        db.executemany(f"INSERT INTO records VALUES ({placeholders})", linked_rows())
    db.execute("CREATE INDEX records_id ON records(chemical_id)")


def write_master(db, path, fields, rows, by_key, stats):
    output_fields = (*fields, "dsstox_inchikey", *(f"dsstox_{field.lower()}" for field in METADATA))
    if len(set(output_fields)) != len(output_fields):
        raise ValueError("Master already contains DSSTox fields")
    aggregates = ", ".join(f"MIN({field}), MAX({field})" for field in STAGED_FIELDS)
    query = f"SELECT MIN(key), chemical_id, {aggregates} FROM records GROUP BY chemical_id ORDER BY MIN(key), chemical_id"
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=output_fields)
        writer.writeheader()
        for key, chemical_id, *extrema in db.execute(query):
            values = {}
            for field, low, high in zip(STAGED_FIELDS, extrema[::2], extrema[1::2]):
                values[field] = low if low == high else None
                stats["metadata_disagreements"] += low != high
            metadata = {f"dsstox_{field.lower()}": values[field] for field in METADATA}
            metadata["dsstox_inchikey"] = key
            standard = key is not None and key[23] == "S"
            stats["nonstandard_key_chemicals"] += key is not None and not standard
            stats["provisional_rows"] += key is None
            row = by_key.get(key)
            if row is None:
                stats["new_chemicals"] += 1
                status = "matched_by_structure" if standard else "dsstox_nonstandard_key"
                row = {
                    "chemical_id": chemical_id, "standard_inchikey": key if standard else "",
                    "representative_smiles": values["SMILES"],
                    "identity_status": status if key else "dsstox_unresolved",
                }
                writer.writerow({**row, **metadata})
            else:
                stats["matched_chemicals"] += 1
                row.update(metadata)
        writer.writerows(rows)
    stats["master_rows"] = len(rows) + stats["new_chemicals"]
    stats["keyed_master_rows"] = len(by_key) + stats["new_chemicals"] - stats["provisional_rows"]


def run(snapshot, source, output):
    started = perf_counter()
    digest = fingerprint(source)
    manifest = snapshot / MANIFEST
    if manifest.is_file():
        previous = json.loads(manifest.read_text())
        if previous["source_sha256"] == digest and previous.get("policy_version") == POLICY_VERSION:
            print(f"DSSTox release already imported in {snapshot}")
            return
        raise ValueError("Use the pre-DSSTox snapshot to rebuild with this release and import policy")
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    fields, rows, by_key, ids = read_master(snapshot / MASTER)
    stats = Counter(dict.fromkeys((
        "input_rows", "missing_key", "malformed_key", "linked_rows", "new_chemicals",
        "matched_chemicals", "metadata_disagreements", "provisional_rows",
        "resolved_by_dtxsid_rows", "nonstandard_key_chemicals",
    ), 0))

    # Stage next to the destination so publishing is a single directory rename.
    with TemporaryDirectory(prefix=".dsstox-", dir=output.parent) as temporary:
        work = Path(temporary)
        staging = work / "snapshot"
        staging.mkdir()
        db = sqlite3.connect(work / "records.sqlite")
        try:
            stage_source(db, source, digest, staging, by_key, ids, stats)
            write_master(db, staging / MASTER, fields, rows, by_key, stats)
        finally:
            db.close()
        if stats["input_rows"] != stats["linked_rows"]:
            raise ValueError("DSSTox row reconciliation failed")
        if fingerprint(source) != digest:
            raise ValueError("DSSTox input changed during the import")

        observations = snapshot / OBSERVATIONS
        copied = staging / OBSERVATIONS
        shutil.copyfile(observations, copied)
        if fingerprint(observations) != fingerprint(copied):
            raise ValueError("Tox21 observations copy failed verification")
        report = {"source_sha256": digest, "policy_version": POLICY_VERSION, **stats}
        (staging / MANIFEST).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        if output.exists() or output.is_symlink():
            raise FileExistsError(output)
        staging.rename(output)
    print(json.dumps(report, indent=2, sort_keys=True))
    print(f"Done in {perf_counter() - started:.2f} seconds.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot", type=Path, help="Existing Tox21 snapshot directory")
    parser.add_argument("source", type=Path, help="Full DSSTox CCD dump CSV")
    parser.add_argument("output", type=Path, help="New snapshot directory")
    args = parser.parse_args()
    run(args.snapshot, args.source, args.output)