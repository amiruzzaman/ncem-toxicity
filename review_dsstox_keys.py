"""Audit a DSSTox key-review queue and write a separate reviewed master.

Usage: python3 review_dsstox_keys.py tox21_dsstox_complete tox21_key_review
Source CSVs are read only. Keys are derived only from original, explicit SMILES.

"""
import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from tempfile import TemporaryDirectory

from rdkit import Chem, rdBase
from rdkit.Chem import rdinchi, rdMolDescriptors

from import_dsstox import KEY_FORMAT, SOURCE_FIELDS, clean, fingerprint

POLICY_VERSION = 1
INPUTS = (
    "master_chemicals.csv", "dsstox_records.csv", "dsstox_review.csv",
    "tox21_observations.csv", "dsstox_import.json",
)
RECORD_FIELDS = ("chemical_id", "source_row_number", "source_sha256", *SOURCE_FIELDS)
REVIEW_FIELDS = ("review_reason", *RECORD_FIELDS)
ANNOTATIONS = (
    "review_category", "next_action", "candidate_standard_inchikey",
    "derived_standard_inchi", "inchi_return_code", "inchi_message",
    "computed_formula", "formula_matches", "roundtrip_matches",
    "local_identifier_key_matches", "existing_key_chemical_ids", "key_update_applied",
)
ACTIONS = {
    "no_source_structure": "Retain substance identity; obtain a curated structure if available",
    "missing_structure_with_dtxcid": "Check the EPA compound record for the missing structure",
    "repeat_or_substance_group": "Retain variable/repeating structure; review its substance-group definition",
    "variable_or_generic_structure": "Retain generic structure; review defined constituents or variable groups",
    "extended_smiles": "Review extended SMILES semantics before deriving a key",
    "smiles_parse_failure": "Review source representation; do not repair valence or aromaticity automatically",
    "inchi_failure": "Review source representation and InChI conversion error",
    "inchi_warning": "Review InChI warning and source structure before accepting a key",
    "roundtrip_changed": "Review structural information lost in the InChI round trip",
    "formula_unverified": "Reconcile source formula with the parsed original structure",
    "disconnected_structure": "Review disconnected components and coordination chemistry",
    "existing_key_match": "Review an identity merge; retain original IDs until a crosswalk is approved",
    "duplicate_recovery_key": "Review candidate identities sharing a derived key",
    "key_derived_from_source_smiles": "Derived key recorded in the reviewed master; original source key remains missing",
}
LOW_PRIORITY = {
    "no_source_structure", "repeat_or_substance_group", "variable_or_generic_structure",
    "key_derived_from_source_smiles",
}


def read_csv(path, expected=None):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream, strict=True)
        fields = tuple(reader.fieldnames or ())
        if not fields or len(fields) != len(set(fields)):
            raise ValueError(f"Missing or duplicate header in {path.name}")
        if expected is not None and fields != tuple(expected):
            raise ValueError(f"Unexpected header in {path.name}")
        for row in reader:
            if None in row or None in row.values():
                raise ValueError(f"Malformed {path.name} row near line {reader.line_num}")
            yield row


def write_csv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def valid_key(value):
    return bool(value and KEY_FORMAT.fullmatch(value))


def classify_structure(row):
    """Generate candidates without treating a parseable generic structure as a molecule."""
    result = dict.fromkeys(ANNOTATIONS, "")
    result["key_update_applied"] = "0"

    def finish(category):
        result["review_category"] = category
        result["next_action"] = ACTIONS[category]
        return result

    smiles = clean(row["SMILES"])
    if not smiles:
        category = "missing_structure_with_dtxcid" if clean(row["DTXCID"]) else "no_source_structure"
        return finish(category)
    # CXSMILES can encode a family or polymer even when its base SMILES is finite.
    if re.search(r"(?:Sg(?:H|D)?|LN):", smiles):
        return finish("repeat_or_substance_group")
    if "*" in smiles or re.search(r"RG:|LOG:|(?:[|,])m:", smiles):
        return finish("variable_or_generic_structure")
    if "|" in smiles:
        return finish("extended_smiles")

    parameters = Chem.SmilesParserParams()
    parameters.parseName = False
    with rdBase.BlockLogs():
        try:
            molecule = Chem.MolFromSmiles(smiles, parameters)
        except (ValueError, RuntimeError) as error:
            result["inchi_message"] = str(error)
            return finish("smiles_parse_failure")
        if molecule is None or molecule.GetNumAtoms() == 0:
            return finish("smiles_parse_failure")
        if any(atom.GetAtomicNum() == 0 or atom.HasQuery() for atom in molecule.GetAtoms()):
            return finish("variable_or_generic_structure")
        if molecule.GetStereoGroups() or Chem.GetMolSubstanceGroups(molecule):
            return finish("extended_smiles")

        result["computed_formula"] = rdMolDescriptors.CalcMolFormula(molecule)
        result["formula_matches"] = str(int(result["computed_formula"] == clean(row["MOLECULAR_FORMULA"])))
        try:
            inchi, code, message, _, _ = rdinchi.MolToInchi(molecule)
            result.update(derived_standard_inchi=inchi, inchi_return_code=str(code), inchi_message=message)
            key = rdinchi.InchiToInchiKey(inchi) if inchi.startswith("InChI=1S/") else ""
            if valid_key(key) and key[23] == "S":
                result["candidate_standard_inchikey"] = key
            if code not in (0, 1) or not result["candidate_standard_inchikey"]:
                return finish("inchi_failure")
            if code != 0 or message.strip():
                return finish("inchi_warning")
            restored, back_code, back_message, _ = rdinchi.InchiToMol(inchi)
        except (ValueError, RuntimeError) as error:
            result["inchi_message"] = str(error)
            return finish("inchi_failure")
        same = restored is not None and back_code == 0 and not back_message.strip()
        if same:
            same = Chem.MolToSmiles(molecule, isomericSmiles=True) == Chem.MolToSmiles(restored, isomericSmiles=True)
        result["roundtrip_matches"] = str(int(same))
        if not same:
            return finish("roundtrip_changed")
        if result["formula_matches"] != "1":
            return finish("formula_unverified")
        if len(Chem.GetMolFrags(molecule)) != 1:
            return finish("disconnected_structure")
    return finish("key_derived_from_source_smiles")


def load_review(snapshot, manifest):
    rows = list(read_csv(snapshot / "dsstox_review.csv", REVIEW_FIELDS))
    ids, numbers = set(), set()
    for row in rows:
        identifier, number = row["chemical_id"], row["source_row_number"]
        if not identifier or identifier in ids or number in numbers:
            raise ValueError("Review queue must have unique, nonempty chemical IDs and source rows")
        if not number.isdecimal() or int(number) < 2:
            raise ValueError("Invalid source row number in review queue")
        if row["source_sha256"] != manifest["source_sha256"]:
            raise ValueError("Review source fingerprint disagrees with import manifest")
        key = clean(row["INCHIKEY"])
        reason = "malformed_key" if key else "missing_key"
        if valid_key(key) or row["review_reason"] != reason:
            raise ValueError("Review reason disagrees with original source key")
        ids.add(identifier)
        numbers.add(number)
    expected = manifest["missing_key"] + manifest["malformed_key"]
    if len(rows) != expected:
        raise ValueError("Review row count disagrees with import manifest")
    return rows


def check_source_records(snapshot, rows, manifest):
    by_number = {row["source_row_number"]: row for row in rows}
    identifiers = {
        field: {clean(row[field]) for row in rows} - {None}
        for field in ("DTXSID", "DTXCID")
    }
    matches = defaultdict(set)
    found, seen_numbers = set(), set()
    count = 0
    for record in read_csv(snapshot / "dsstox_records.csv", RECORD_FIELDS):
        count += 1
        number = record["source_row_number"]
        if number in by_number:
            if number in seen_numbers or any(record[field] != by_number[number][field] for field in RECORD_FIELDS):
                raise ValueError(f"Review row {number} is duplicated or differs from linked source record")
            found.add(number)
            seen_numbers.add(number)
        if record["source_sha256"] != manifest["source_sha256"]:
            raise ValueError("Linked record fingerprint disagrees with import manifest")
        key = clean(record["INCHIKEY"])
        if valid_key(key):
            for field in identifiers:
                identifier = clean(record[field])
                if identifier in identifiers[field]:
                    matches[field, identifier].add(key)
    if found != set(by_number) or count != manifest["linked_rows"]:
        raise ValueError("Linked source records do not reconcile with the review queue or manifest")
    return matches


def check_master(snapshot, rows, results, manifest):
    by_id = {row["chemical_id"]: row for row in rows}
    candidate_keys = {r["candidate_standard_inchikey"] for r in results.values()} - {""}
    key_matches = defaultdict(set)
    found = set()
    count = 0
    with (snapshot / "master_chemicals.csv").open(newline="", encoding="utf-8-sig") as stream:
        fields = tuple(next(csv.reader(stream)))
    required = {"chemical_id", "standard_inchikey", "identity_status", "representative_smiles", "dsstox_inchikey", "dsstox_dtxsid"}
    if not required.issubset(fields):
        raise ValueError("Master is missing identity columns")
    for row in read_csv(snapshot / "master_chemicals.csv", fields):
        count += 1
        identifier = row["chemical_id"]
        if identifier in by_id:
            source = by_id[identifier]
            if identifier in found or row["identity_status"] != "dsstox_unresolved":
                raise ValueError("Review master identity is duplicated or no longer unresolved")
            if clean(row["standard_inchikey"]) or clean(row["dsstox_inchikey"]):
                raise ValueError("Review master entry already has a key")
            if (clean(row["representative_smiles"]) != clean(source["SMILES"])
                    or clean(row["dsstox_dtxsid"]) != clean(source["DTXSID"])):
                raise ValueError("Review structure or substance ID disagrees with the master")
            found.add(identifier)
        for field in ("standard_inchikey", "dsstox_inchikey"):
            key = row[field]
            if key in candidate_keys:
                key_matches[key].add(identifier)
    if found != set(by_id) or count != manifest["master_rows"]:
        raise ValueError("Master does not reconcile with review identities or import manifest")
    return fields, key_matches


def choose_updates(rows, results, identifier_matches, key_matches):
    proposed = Counter(
        r["candidate_standard_inchikey"] for r in results.values()
        if r["review_category"] == "key_derived_from_source_smiles"
    )
    updates = {}
    for row in rows:
        identifier = row["chemical_id"]
        result = results[identifier]
        local = {
            field: sorted(identifier_matches.get((field, clean(row[field])), ()))
            for field in ("DTXSID", "DTXCID")
        }
        local = {field: keys for field, keys in local.items() if keys}
        result["local_identifier_key_matches"] = json.dumps(local, sort_keys=True) if local else ""
        key = result["candidate_standard_inchikey"]
        result["existing_key_chemical_ids"] = ";".join(sorted(key_matches.get(key, ())))
        if result["review_category"] == "key_derived_from_source_smiles":
            if key_matches.get(key):
                result["review_category"] = "existing_key_match"
            elif proposed[key] > 1:
                result["review_category"] = "duplicate_recovery_key"
            elif local and any(keys != [key] for keys in local.values()):
                raise ValueError("Generated key conflicts with a local substance or compound identifier")
            else:
                result["key_update_applied"] = "1"
                updates[identifier] = key
            result["next_action"] = ACTIONS[result["review_category"]]
    return updates


def write_reviewed_master(snapshot, path, fields, updates):
    changed = 0

    def reviewed_rows():
        nonlocal changed
        for row in read_csv(snapshot / "master_chemicals.csv", fields):
            if row["chemical_id"] in updates:
                row["standard_inchikey"] = updates[row["chemical_id"]]
                row["identity_status"] = "key_derived_from_source_smiles"
                changed += 1
            yield row

    write_csv(path, fields, reviewed_rows())
    if changed != len(updates):
        raise ValueError("Reviewed master update count did not reconcile")


def run(snapshot, output):
    snapshot, output = Path(snapshot), Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if output.resolve().is_relative_to(snapshot.resolve()):
        raise ValueError("Write the review beside the source snapshot, not inside it")
    input_hashes = {name: fingerprint(snapshot / name) for name in INPUTS}
    manifest = json.loads((snapshot / "dsstox_import.json").read_text())
    rows = load_review(snapshot, manifest)
    results = {row["chemical_id"]: classify_structure(row) for row in rows}
    identifier_matches = check_source_records(snapshot, rows, manifest)
    master_fields, key_matches = check_master(snapshot, rows, results, manifest)
    updates = choose_updates(rows, results, identifier_matches, key_matches)
    counts = Counter(r["review_category"] for r in results.values())
    priority_count = sum(count for category, count in counts.items() if category not in LOW_PRIORITY)
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".key-review-", dir=output.parent) as temporary:
        staging = Path(temporary) / "review"
        staging.mkdir()
        fields = (*REVIEW_FIELDS, *ANNOTATIONS)
        annotated = [{**row, **results[row["chemical_id"]]} for row in rows]
        write_csv(staging / "key_review.csv", fields, annotated)
        write_csv(staging / "priority_review.csv", fields,
                  (r for r in annotated if r["review_category"] not in LOW_PRIORITY))
        update_fields = (
            "chemical_id", "DTXSID", "PREFERRED_NAME", "source_row_number", "source_sha256",
            "SMILES", "candidate_standard_inchikey", "derived_standard_inchi",
            "review_category", "formula_matches", "roundtrip_matches", "key_update_applied",
        )
        write_csv(staging / "key_updates.csv", update_fields,
                  ({field: row[field] for field in update_fields} for row in annotated if row["key_update_applied"] == "1"))
        write_reviewed_master(snapshot, staging / "master_chemicals_reviewed.csv", master_fields, updates)
        if input_hashes != {name: fingerprint(snapshot / name) for name in INPUTS}:
            raise ValueError("Source snapshot changed during review")
        report = {
            "review_policy_version": POLICY_VERSION,
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "rdkit_version": rdBase.rdkitVersion, "inchi_version": rdinchi.GetInchiVersion(),
            "dsstox_import_policy_version": manifest.get("policy_version"),
            "dsstox_release_sha256": manifest["source_sha256"],
            "input_sha256": input_hashes,
            "master_rows": manifest["master_rows"], "review_rows": len(rows),
            "derived_keys_applied": len(updates), "remaining_missing_keys": len(rows) - len(updates),
            "priority_review_rows": priority_count,
            "rows_with_local_identifier_key_matches": sum(bool(r["local_identifier_key_matches"]) for r in results.values()),
            "review_categories": dict(sorted(counts.items())),
            "checks": {
                "source_snapshot_unchanged": True, "review_matches_linked_source": True,
                "review_matches_master": True, "chemical_ids_preserved": True,
                "master_row_count_preserved": True,
            },
            "output_sha256": {p.name: fingerprint(p) for p in sorted(staging.iterdir())},
        }
        (staging / "key_review_manifest.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        if output.exists() or output.is_symlink():
            raise FileExistsError(output)
        staging.rename(output)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot", type=Path, help="Combined DSSTox snapshot containing its review queue")
    parser.add_argument("output", type=Path, help="New directory for reviewed CSVs and manifest")
    args = parser.parse_args()
    report = run(args.snapshot, args.output)
    print(json.dumps({key: report[key] for key in (
        "review_rows", "derived_keys_applied", "remaining_missing_keys", "priority_review_rows", "review_categories",
    )}, indent=2))
