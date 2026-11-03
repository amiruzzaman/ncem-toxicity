"""
Summarize Tox21 labels without changing input.
Prints the report to standard output.
use this in cmdln: 
python3 summary_tox21.py tox21.csv
||
if you desire a report:
python3 summary_tox21.py tox21.csv > tox21_summary.txt
"""

import argparse
import csv
import hashlib
import io
import sys
from decimal import Decimal, InvalidOperation
from pathlib import Path

ENDPOINTS = (
    "NR-AR", "NR-AR-LBD", "NR-AhR", "NR-Aromatase", "NR-ER",
    "NR-ER-LBD", "NR-PPAR-gamma", "SR-ARE", "SR-ATAD5", "SR-HSE",
    "SR-MMP", "SR-p53",
)

MISSING = {"", "na", "n/a", "nan", "null", "none"}

def classify(val):
    ## Classifies a value as "active", "inactive", "missing", or "invalid" based on its content
    val = val.strip().lower()
    if val in MISSING:
        return "missing"
    try:
        val = Decimal(val)
    except InvalidOperation:
        return "invalid"
    if not val.is_finite():
        return "invalid"
    ## inifinity check just incase.

    if val == 1:
        return "active"
    elif val == 0:
        return "inactive"
    else:
        return "invalid"

def summarize(path):

    raw = path.read_bytes()
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig"), newline=""), strict=True)
    headers = reader.fieldnames

    """
    Checks if the CSV file has headers and if there are any duplicate headers. Raises a ValueError if either condition is not met.
    """

    if not headers:
        raise ValueError("No headers found in CSV file")
    if len(headers) != len(set(headers)):
        raise ValueError("Duplicate headers found in CSV file")

    required = (*ENDPOINTS, "mol_id", "smiles")
    ## req columns are the endpoints plus the two identifiers
    absent = [name for name in required if name not in headers]

    if absent:
        raise ValueError(f"Missing required columns: {', '.join(absent)}")

    counts = {name: dict(active=0, inactive=0, invalid=0, missing=0) for name in ENDPOINTS}
    ## counts is a dictionary mapping each endpoint to a dictionary of its label counts
    missing_by_column = dict.fromkeys(headers, 0)

    total = 0
    for row in reader:
        if None in row or any (value is None for value in row.values()):
            raise ValueError(f"Wrong number of fields near CSV line {reader.line_num}.")
        total += 1
        for name in headers:
            if row[name].strip().lower() in MISSING:
                missing_by_column[name] += 1
        for name in ENDPOINTS:
            counts[name][classify(row[name])] += 1
    ## The above loop iterates through each row in the CSV file, counting the occurrences of each label category for each endpoint and tracking missing values by column for reporting
    lines = [
        "Tox21 data summary:",
        f"Input SHA-256: {hashlib.sha256(raw).hexdigest()}",
        f"Rows: {total:,}",
        f"Columns: {len(headers)}",
        f"Endpoints: {len(ENDPOINTS)}",
        "",
        "Labels 1 = active; 0 = inactive in the specified assay.",
        "These labels do not establish overall human toxicity or safety",
        "Missing: blank, NA, N/A, NaN, null, none (case-insensitive).",
        "Invalid: any other label. Numeric forms such as 0.0 or 1.0 are accepted.",
        "",
        f"{'Endpoint':<16} {'Active':>8} {'Inactive':>9} {'Missing':>8} {'Invalid':>8} {'Missing %':>10}",
    ]

    for name, values in counts.items():
        assert sum(values.values()) == total
        percentage = f"{100 * values['missing'] / total:.2f}%" if total else "N/A"
        lines.append(
            f"{name:<16} {values['active']:>8,} {values['inactive']:>9,}  "
            f"{values['missing']:>8,} {values['invalid']:>8,} {percentage:>10}"
        )

    totals = {category: sum(values[category] for values in counts.values())
              for category in ("active", "inactive", "missing", "invalid")}
    ## totals is a dictionary mapping each label category to the total count across all endpoints
    lines.extend([
        "",
        f"Total endpoint cells: {total * len(ENDPOINTS):,}",
        "Endpoint totals: " + "; ".join(f"{k}={v:,}" for k, v in totals.items()),
        "Check: each endpoint's four counts sum to the row count: PASS",
        "",
        "Missing values by column (including identifiers):",
    ])

    lines.extend(f"{name}: {missing_by_column[name]:,}" for name in headers)
    if total == 0:
        lines.append("No rows found in the input file.")
    if totals["invalid"] > 0:
        lines.append(
            f"Warning: {totals['invalid']:,} invalid endpoint values found. "
            "These may be due to typos or unexpected formatting."
        )
    return "\n".join(lines) + "\n", totals["invalid"] > 0

def main():
    # reads csv filename supplied on cmdline. ex: python summary_tox21.py tox21.csv
    
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "csv_path", 
        type=Path, 
        help="Path to Tox21 CSV file.",
    )

    args = parser.parse_args()

    try:
        report, invalid = summarize(args.csv_path)
    except (OSError, UnicodeError, csv.Error, ValueError) as error:
        print(f"Error processing file {args.csv_path}: {error}", file=sys.stderr)
        return 1

    print(report, end="")
    return 2 if invalid else 0

if __name__ == "__main__":
    sys.exit(main())

