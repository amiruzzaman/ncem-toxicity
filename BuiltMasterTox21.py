"""builds initial masterchem csv from Tox21."""
import argparse
import csv
import hashlib
import io
from collections import defaultdict
from decimal import Decimal
from pathlib import Path
from uuid import NAMESPACE_URL, uuid5

from rdkit import Chem, RDLogger, rdBase

ENDPOINTS = (
    "NR-AR", "NR-AR-LBD", "NR-AhR", "NR-Aromatase", "NR-ER", "NR-ER-LBD",
    "NR-PPAR-gamma", "SR-ARE", "SR-ATAD5", "SR-HSE", "SR-MMP", "SR-p53",
)
MISSING = {"", "na", "n/a", "nan", "null", "none"}


def write_csv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
## helper func to write to csv, takes path, fields, and rows as input and writes 
# the rows to the csv file at the specified path with the given fields as headers.

def build(source, output):
    raw = source.read_bytes()
    fingerprint = hashlib.sha256(raw).hexdigest()
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig"), newline=""), strict=True)
    required = {*ENDPOINTS, "mol_id", "smiles"}

    if not required.issubset(reader.fieldnames or []):
        raise ValueError("Input is missing required Tox21 columns")

    chemicals = {}
    observations = {}
    values = defaultdict(set)
    RDLogger.DisableLog("rdApp.warning")
    for row in reader:
        if None in row or any(row[column] is None for column in required):
            raise ValueError(f"Malformed CSV row at line {reader.line_num}")
        
        source_id = row["mol_id"].strip()
        smiles = row["smiles"]

        if not source_id:
            raise ValueError(f"Missing mol_id at line {reader.line_num}")
        
        # uses chem identifiers to properly identify and assign chemicals and observations, 
        # ensuring that each chemical is uniquely identified by its InChIKey or a generated UUID if the InChIKey is not available.

        molecule = Chem.MolFromSmiles(smiles) if smiles.strip() else None
        key = Chem.MolToInchiKey(molecule) if molecule is not None else ""
        identity = key or f"unresolved:{source_id}:{smiles}"
        chemical_id = str(uuid5(NAMESPACE_URL, f"researchchem:chemical:{identity}"))

        # use full struc key if available, otherwise use a generated UUID 
        # based on the source_id and smiles to ensure uniqueness.
        chemicals.setdefault(chemical_id, {
            "chemical_id": chemical_id, "standard_inchikey": key,
            "representative_smiles": smiles,
            "identity_status": "matched_by_structure" if key else "needs_review",
        })
        # keep unresolved structs under source id for later review incase.

        for endpoint in ENDPOINTS:
            raw_label = row[endpoint].strip()
            if raw_label.lower() in MISSING:
                continue

            number = Decimal(raw_label)
            if not number.is_finite() or number not in (0, 1):
                raise ValueError(f"Invalid {endpoint} label for {source_id}: {raw_label!r}")
            
            result = int(number)
            observation_key = (chemical_id, source_id, endpoint, result)
            observations.setdefault(observation_key, {
                "chemical_id": chemical_id, "source_dataset": "Tox21",
                "source_record_id": source_id, "source_smiles": smiles,
                "source_sha256": fingerprint, "endpoint": endpoint, "result": result,
            })

            values[chemical_id, endpoint].add(result)

    conflicts = 0
    for chemical_id, chemical in chemicals.items():
        for endpoint in ENDPOINTS:
            results = values.get((chemical_id, endpoint), set())
            conflicts += len(results) > 1
            chemical[endpoint] = "conflict" if len(results) > 1 else next(iter(results), "")

    output.mkdir(parents=True, exist_ok=False)
    write_csv(output / "master_chemicals.csv",
              ("chemical_id", "standard_inchikey", "representative_smiles",
               "identity_status", *ENDPOINTS), chemicals.values())
    
    write_csv(output / "tox21_observations.csv",
              ("chemical_id", "source_dataset", "source_record_id", "source_smiles",
               "source_sha256", "endpoint", "result"), observations.values())
    
    unresolved = sum(row["identity_status"] == "needs_review" for row in chemicals.values())

    print(f"RDKit {rdBase.rdkitVersion} | {len(chemicals):,} chemicals | "
          f"{len(observations):,} observations | {conflicts:,} conflicting cells | "
          f"{unresolved:,} unresolved chemical rows")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path, help="New directory for the two CSVs")
    arguments = parser.parse_args()
    build(arguments.source, arguments.output)
