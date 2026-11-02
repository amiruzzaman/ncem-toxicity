"""
summmarize tox21 labels without changing input

saves printed report to tox21_summary.txt
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

    if val.lower() in MISSING:
        return "missing"
    try:
        val = Decimal(val)
    except InvalidOperation:
        return "invalid"
    if val == 1:
        return "active"
    elif val == 0:
        return "inactive"
    else:
        return "invalid"

"""def summarize(path):

    raw = path.readbytes()
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8"), newline="", strict=True))
    headers = reader.fieldnames

    if not headers:
        raise ValueError("No headers found in CSV file")
    if len(headers) != len(set(headers)):
        raise ValueError("Duplicate headers found in CSV file")  """  