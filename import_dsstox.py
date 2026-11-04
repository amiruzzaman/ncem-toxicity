"""

# From the repo root, run: python3 import_dsstox.py tox21_initial "/path/to/DSSToxCCDdump.csv" tox21_dsstox

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
from uuid import NAMESPACE_URL, uuid5

MASTER = "master_chemicals.csv"
OBSERVATIONS = "tox21_observations.csv"
MANIFEST = "dsstox_import.json"
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


def stage_source(db, source, digest, directory, by_key, ids, stats):
    columns = ", ".join(f"{field} TEXT" for field in STAGED_FIELDS)
    db.execute(f"CREATE TABLE records (key TEXT, chemical_id TEXT, {columns})")
    placeholders = ", ".join("?" for _ in range(2 + len(STAGED_FIELDS)))
    with (
        source.open(newline="", encoding="utf-8-sig") as stream,
        (directory / "dsstox_records.csv").open("w", newline="", encoding="utf-8") as linked,
        (directory / "dsstox_review.csv").open("w", newline="", encoding="utf-8") as review,
    ):
        reader = csv.DictReader(stream, strict=True)
        if tuple(reader.fieldnames or ()) != SOURCE_FIELDS:
            raise ValueError("Unexpected DSSTox header")
        linked_writer, review_writer = csv.writer(linked), csv.writer(review)
        provenance = ("source_row_number", "source_sha256")
        linked_writer.writerow(("chemical_id", *provenance, *SOURCE_FIELDS))
        review_writer.writerow(("review_reason", *provenance, *SOURCE_FIELDS))

        def linked_rows():
            # Source row numbers count logical CSV records, with the header at 1.
            for number, row in enumerate(reader, start=2):
                if None in row or None in row.values():
                    raise ValueError(f"Malformed DSSTox record {number}")
                stats["input_rows"] += 1
                raw = (number, digest, *(row[field] for field in SOURCE_FIELDS))
                key = clean(row["INCHIKEY"])
                if not key or not KEY_FORMAT.fullmatch(key):
                    reason = "malformed_key" if key else "missing_key"
                    stats[reason] += 1
                    review_writer.writerow((reason, *raw))
                    continue

                existing = by_key.get(key)
                if existing is None:
                    chemical_id = str(uuid5(NAMESPACE_URL, f"researchchem:chemical:{key}"))
                    if chemical_id in ids:
                        raise ValueError(f"Generated chemical ID is already in use for {key}")
                else:
                    chemical_id = existing["chemical_id"]
                stats["linked_rows"] += 1
                linked_writer.writerow((chemical_id, *raw))
                yield (key, chemical_id, *(clean(row[field]) for field in STAGED_FIELDS))

        db.executemany(f"INSERT INTO records VALUES ({placeholders})", linked_rows())
    db.execute("CREATE INDEX records_key ON records(key)")


def write_master(db, path, fields, rows, by_key, stats):
    output_fields = (*fields, "dsstox_inchikey", *(f"dsstox_{field.lower()}" for field in METADATA))
    if len(set(output_fields)) != len(output_fields):
        raise ValueError("Master already contains DSSTox fields")
    aggregates = ", ".join(f"MIN({field}), MAX({field})" for field in STAGED_FIELDS)
    query = f"SELECT key, MIN(chemical_id), {aggregates} FROM records GROUP BY key ORDER BY key"
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
            standard = key[23] == "S"
            stats["nonstandard_key_chemicals"] += not standard
            row = by_key.get(key)
            if row is None:
                stats["new_chemicals"] += 1
                row = {
                    "chemical_id": chemical_id, "standard_inchikey": key if standard else "",
                    "representative_smiles": values["SMILES"],
                    "identity_status": "matched_by_structure" if standard else "dsstox_nonstandard_key",
                }
                writer.writerow({**row, **metadata})
            else:
                stats["matched_chemicals"] += 1
                row.update(metadata)
        writer.writerows(rows)
    stats["master_rows"] = len(rows) + stats["new_chemicals"]


def run(snapshot, source, output):
    digest = fingerprint(source)
    manifest = snapshot / MANIFEST
    if manifest.is_file():
        if json.loads(manifest.read_text())["source_sha256"] == digest:
            print(f"DSSTox release already imported in {snapshot}")
            return
        raise ValueError("Use the pre-DSSTox snapshot to build a different DSSTox release")
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    fields, rows, by_key, ids = read_master(snapshot / MASTER)
    stats = Counter(dict.fromkeys((
        "input_rows", "missing_key", "malformed_key", "linked_rows", "new_chemicals",
        "matched_chemicals", "metadata_disagreements",
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
        if stats["input_rows"] != stats["linked_rows"] + stats["missing_key"] + stats["malformed_key"]:
            raise ValueError("DSSTox row reconciliation failed")
        if fingerprint(source) != digest:
            raise ValueError("DSSTox input changed during the import")

        observations = snapshot / OBSERVATIONS
        copied = staging / OBSERVATIONS
        shutil.copyfile(observations, copied)
        if fingerprint(observations) != fingerprint(copied):
            raise ValueError("Tox21 observations copy failed verification")
        report = {"source_sha256": digest, **stats}
        (staging / MANIFEST).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        if output.exists() or output.is_symlink():
            raise FileExistsError(output)
        staging.rename(output)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot", type=Path, help="Existing Tox21 snapshot directory")
    parser.add_argument("source", type=Path, help="Full DSSTox CCD dump CSV")
    parser.add_argument("output", type=Path, help="New snapshot directory")
    args = parser.parse_args()
    run(args.snapshot, args.source, args.output)
