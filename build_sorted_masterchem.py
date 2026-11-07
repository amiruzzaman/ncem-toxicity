"""Export a sorted, assay-filtered CSV from a published three-source snapshot."""

# Dev notes
# From the repository root, use Python 3.11+ with requirements/toxcast.txt installed.
# python3 build_sorted_masterchem.py \
#     --snapshot tox21_dsstox_toxcast_v4_3 \
#     --output "$HOME/Documents/ChatGPT/ResearchChem/SortedMasterChem-rerun.csv"
# Choose an unused output path for each run.

import argparse
import csv
import json
from collections import defaultdict
from functools import cache
from itertools import zip_longest
from pathlib import Path
from tempfile import TemporaryDirectory

from download_toxcast import fingerprint
from import_toxcast import SUMMARY_FIELDS, activity_status, csv_rows, identifier

TOX21_ASSAYS = (
    "NR-AR", "NR-AR-LBD", "NR-AhR", "NR-Aromatase", "NR-ER", "NR-ER-LBD",
    "NR-PPAR-gamma", "SR-ARE", "SR-ATAD5", "SR-HSE", "SR-MMP", "SR-p53",
)
TOX21_COUNT = "tox21_usable_assay_count"
TOXCAST_COUNT = "toxcast_usable_endpoint_count"
MIN_USABLE = 6
MASTER_NAME = "master_chemicals.csv"
OBSERVATIONS = (
    ("mc", "toxcast_mc_observations.csv"),
    ("sc", "toxcast_sc_observations.csv"),
)
CALL_BITS = {"active": 1, "inactive": 2}


def input_hashes(snapshot):
    manifest = json.loads((snapshot / "toxcast_import.json").read_text(encoding="utf-8"))
    output_hashes = manifest.get("output_sha256") if isinstance(manifest, dict) else None
    names = (MASTER_NAME, *(name for _, name in OBSERVATIONS))
    if not isinstance(output_hashes, dict):
        raise ValueError("Missing output_sha256 in toxcast_import.json")
    for name in names:
        if not output_hashes.get(name):
            raise ValueError(f"Missing manifest hash for {name}")
    return {snapshot / name: output_hashes[name] for name in names}


def verify_inputs(expected_hashes, phase):
    for path, expected in expected_hashes.items():
        if fingerprint(path) != expected:
            raise ValueError(f"{path.name}: SHA-256 mismatch {phase}")


def toxcast_counts(snapshot):
    by_chemical = defaultdict(dict)
    endpoint_id = cache(identifier)
    for screen, name in OBSERVATIONS:
        with csv_rows(snapshot / name) as (fields, rows):
            required = {"chemical_id", "aeid", "activity_status", "result_id"}
            if screen == "mc":
                required.add("modl")
            if not required.issubset(fields) or not {"hitc", "hitcall"}.intersection(fields):
                raise ValueError(f"Missing observation columns in {name}")

            for number, row in rows:
                status = activity_status(row, screen)
                if status != row["activity_status"]:
                    raise ValueError(
                        f"{name}, logical row {number}, result {row['result_id']}: "
                        f"activity_status is {row['activity_status']!r}, expected {status!r}"
                    )
                chemical_id = row["chemical_id"]
                if not chemical_id.strip():
                    continue
                endpoint_calls = by_chemical[chemical_id]
                if status in CALL_BITS:
                    aeid = endpoint_id(row["aeid"])
                    if not aeid:
                        raise ValueError(f"Missing aeid in {name}, logical row {number}")
                    # Seeing both bits makes the endpoint a disagreement across all screens.
                    endpoint_calls[aeid] = endpoint_calls.get(aeid, 0) | CALL_BITS[status]

    return {
        chemical_id: sum(bits in (1, 2) for bits in endpoint_calls.values())
        for chemical_id, endpoint_calls in by_chemical.items()
    }


def qualifying_rows(master, cast_counts):
    seen = set()
    qualified = []
    totals = dict.fromkeys(("tox21_only", "toxcast_only", "both"), 0)
    with csv_rows(master) as (fields, rows):
        required = {"chemical_id", *TOX21_ASSAYS, *SUMMARY_FIELDS}
        if (len(fields) != 27 or not required.issubset(fields)
                or {TOX21_COUNT, TOXCAST_COUNT}.intersection(fields)):
            raise ValueError(f"Unexpected 27-column master schema in {master.name}")

        for number, row in rows:
            chemical_id = row["chemical_id"]
            if not chemical_id.strip() or chemical_id in seen:
                raise ValueError(f"Missing or duplicate chemical_id in {master.name}, logical row {number}")
            seen.add(chemical_id)
            tox21 = sum(row[field] in ("0", "1") for field in TOX21_ASSAYS)
            toxcast = cast_counts.pop(chemical_id, 0)
            if tox21 < MIN_USABLE and toxcast < MIN_USABLE:
                continue

            if tox21 >= MIN_USABLE and toxcast >= MIN_USABLE:
                totals["both"] += 1
            elif tox21 >= MIN_USABLE:
                totals["tox21_only"] += 1
            else:
                totals["toxcast_only"] += 1
            row[TOX21_COUNT] = str(tox21)
            row[TOXCAST_COUNT] = str(toxcast)
            qualified.append(row)

    if cast_counts:
        raise ValueError(f"ToxCast observations reference {len(cast_counts)} IDs outside the master")

    qualified.sort(key=lambda row: (
        -max(int(row[TOX21_COUNT]), int(row[TOXCAST_COUNT])),
        row["chemical_id"],
    ))
    return (*fields, TOX21_COUNT, TOXCAST_COUNT), qualified, {
        "master_rows": len(seen),
        "exported_rows": len(qualified),
        **totals,
    }


def validate_staged_csv(staged, fields, expected_rows):
    with csv_rows(staged) as (actual_fields, rows):
        if actual_fields != fields:
            raise ValueError("Staged CSV header differs from the export schema")
        actual_rows = (row for _, row in rows)
        for number, (expected, actual) in enumerate(zip_longest(expected_rows, actual_rows), 2):
            if expected != actual:
                raise ValueError(f"Staged CSV differs from the validated export at logical row {number}")


def build(snapshot, output):
    snapshot = Path(snapshot).resolve(strict=True)
    output = Path(output)
    output = output.parent.resolve(strict=True) / output.name
    if output.is_relative_to(snapshot):
        raise ValueError("The output path must be outside the input snapshot")
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)

    expected_hashes = input_hashes(snapshot)
    verify_inputs(expected_hashes, "before reading")
    fields, rows, summary = qualifying_rows(snapshot / MASTER_NAME, toxcast_counts(snapshot))

    with TemporaryDirectory(prefix=".sorted-masterchem-", dir=output.parent) as directory:
        staged = Path(directory) / output.name
        with staged.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        validate_staged_csv(staged, fields, rows)
        verify_inputs(expected_hashes, "after reading")
        # Linking on the destination filesystem publishes the complete CSV without replacing a file.
        output.hardlink_to(staged)

    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", type=Path, required=True, help="Published three-source snapshot")
    parser.add_argument("--output", type=Path, required=True, help="New CSV path outside the snapshot")
    args = parser.parse_args()
    try:
        summary = build(args.snapshot, args.output)
    except (OSError, ValueError, csv.Error) as error:
        parser.error(str(error))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
