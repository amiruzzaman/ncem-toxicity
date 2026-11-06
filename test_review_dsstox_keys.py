"""Scientific identity and provenance checks for the key-review workflow."""
import csv
import json
import tempfile
import unittest
from pathlib import Path

from review_dsstox_keys import (
    INPUTS, RECORD_FIELDS, REVIEW_FIELDS, classify_structure, fingerprint, run, write_csv,
)

ETHANOL_KEY = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
DIGEST = "a" * 64


def source_row(number=2, smiles="CCO", formula="C2H6O", **changes):
    row = dict.fromkeys(REVIEW_FIELDS, "")
    row.update(
        review_reason="missing_key", chemical_id=f"chemical-{number}",
        source_row_number=str(number), source_sha256=DIGEST,
        DTXSID=f"DTXSID{number}", DTXCID=f"DTXCID{number}",
        PREFERRED_NAME=f"Test substance {number}", CASRN="NOCAS_TEST",
        SMILES=smiles, MOLECULAR_FORMULA=formula,
    )
    row.update(changes)
    return row


def make_snapshot(directory, reviews, existing_ethanol=False):
    snapshot = directory / "source"
    snapshot.mkdir()
    records = [{field: row[field] for field in RECORD_FIELDS} for row in reviews]
    fields = (
        "chemical_id", "standard_inchikey", "representative_smiles", "identity_status",
        "NR-ER", "dsstox_inchikey", "dsstox_dtxsid", "dsstox_preferred_name",
    )
    master = [dict(zip(fields, (
        row["chemical_id"], "", row["SMILES"], "dsstox_unresolved", "",
        "", row["DTXSID"], row["PREFERRED_NAME"],
    ))) for row in reviews]
    master.append(dict(zip(fields, (
        "control", ETHANOL_KEY if existing_ethanol else "OTMSDBZUPAUEDD-UHFFFAOYSA-N",
        "CCO" if existing_ethanol else "CC", "matched_by_structure", "conflict", "", "", "",
    ))))
    write_csv(snapshot / "master_chemicals.csv", fields, master)
    write_csv(snapshot / "dsstox_records.csv", RECORD_FIELDS, records)
    write_csv(snapshot / "dsstox_review.csv", REVIEW_FIELDS, reviews)
    write_csv(snapshot / "tox21_observations.csv", ("chemical_id", "endpoint", "result"), (
        dict(chemical_id="control", endpoint="NR-ER", result=value) for value in (0, 1)
    ))
    (snapshot / "dsstox_import.json").write_text(json.dumps({
        "source_sha256": DIGEST, "policy_version": 2, "master_rows": len(master),
        "linked_rows": len(records), "missing_key": len(reviews), "malformed_key": 0,
    }))
    return snapshot


def csv_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


class StructureTests(unittest.TestCase):
    def test_explicit_structure_recovers_standard_key(self):
        result = classify_structure(source_row())
        self.assertEqual(result["review_category"], "key_derived_from_source_smiles")
        self.assertEqual(result["candidate_standard_inchikey"], ETHANOL_KEY)
        self.assertEqual(result["formula_matches"], "1")
        self.assertEqual(result["roundtrip_matches"], "1")

    def test_generic_and_repeat_structures_never_become_finite_molecules(self):
        examples = (
            ("CCO |Sg:n:0,1,2::ht|", "repeat_or_substance_group"),
            ("CCCC |LN:1:1.7|", "repeat_or_substance_group"),
            ("CC*", "variable_or_generic_structure"),
            ("CCO |lp:2:2|", "extended_smiles"),
        )
        for smiles, category in examples:
            with self.subTest(smiles=smiles):
                result = classify_structure(source_row(smiles=smiles))
                self.assertEqual(result["review_category"], category)
                self.assertEqual(result["candidate_standard_inchikey"], "")

    def test_prepared_smiles_cannot_replace_missing_original_structure(self):
        result = classify_structure(source_row(smiles="", DTXCID="", MS_READY_SMILES="CCO"))
        self.assertEqual(result["review_category"], "no_source_structure")
        self.assertEqual(result["candidate_standard_inchikey"], "")
        self.assertEqual(classify_structure(source_row(smiles=""))["review_category"],
                         "missing_structure_with_dtxcid")

    def test_formula_disagreement_stays_in_review(self):
        self.assertEqual(classify_structure(source_row(formula="CH4O"))["review_category"],
                         "formula_unverified")

    def test_unparseable_source_is_not_repaired(self):
        result = classify_structure(source_row(smiles="c1cccc1", formula="C5H5"))
        self.assertEqual(result["review_category"], "smiles_parse_failure")

    def test_inchi_warning_stays_in_review(self):
        result = classify_structure(source_row(smiles="CC(F)Br", formula="C2H4BrF"))
        self.assertEqual(result["review_category"], "inchi_warning")
        self.assertTrue(result["inchi_message"])

    def test_roundtrip_loss_stays_in_review(self):
        result = classify_structure(source_row(smiles="[Cu].NC1=NN=NN1.NC1=NN=NN1", formula="C2H6CuN10"))
        self.assertEqual(result["review_category"], "roundtrip_changed")

    def test_disconnected_coordination_structure_requires_review(self):
        result = classify_structure(source_row(
            smiles="[Cr].C1=CC=CC=C1.FC1=C(F)C(F)=C(F)C(F)=C1F", formula="C12H6CrF6"))
        self.assertEqual(result["review_category"], "disconnected_structure")


class WorkflowTests(unittest.TestCase):
    def test_only_accepted_fields_change_and_provenance_is_preserved(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reviews = [source_row(), source_row(3, "CC*", "C2"), source_row(4, "", "", DTXCID="")]
            snapshot = make_snapshot(root, reviews)
            before = {name: fingerprint(snapshot / name) for name in INPUTS}
            output = root / "reviewed"
            report = run(snapshot, output)
            self.assertEqual(report["derived_keys_applied"], 1)
            self.assertEqual(report["remaining_missing_keys"], 2)
            self.assertEqual(before, {name: fingerprint(snapshot / name) for name in INPUTS})
            original = csv_rows(snapshot / "master_chemicals.csv")
            reviewed = csv_rows(output / "master_chemicals_reviewed.csv")
            expected = [dict(row) for row in original]
            expected[0].update(standard_inchikey=ETHANOL_KEY, identity_status="key_derived_from_source_smiles")
            self.assertEqual(reviewed, expected)
            self.assertEqual(csv_rows(output / "key_updates.csv")[0]["chemical_id"], "chemical-2")
            for source, annotated in zip(reviews, csv_rows(output / "key_review.csv")):
                self.assertEqual(source, {field: annotated[field] for field in REVIEW_FIELDS})
            for name, digest in report["output_sha256"].items():
                self.assertEqual(fingerprint(output / name), digest)
            with self.assertRaises(FileExistsError):
                run(snapshot, output)

    def test_existing_identity_and_duplicate_candidates_are_not_merged(self):
        for existing in (False, True):
            with self.subTest(existing=existing), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                reviews = [source_row()] if existing else [source_row(), source_row(3)]
                snapshot = make_snapshot(root, reviews, existing_ethanol=existing)
                output = root / "reviewed"
                report = run(snapshot, output)
                self.assertEqual(report["derived_keys_applied"], 0)
                self.assertEqual(csv_rows(output / "master_chemicals_reviewed.csv"),
                                 csv_rows(snapshot / "master_chemicals.csv"))
                category = "existing_key_match" if existing else "duplicate_recovery_key"
                self.assertTrue(all(row["review_category"] == category
                                    for row in csv_rows(output / "key_review.csv")))
                self.assertEqual(len(csv_rows(output / "priority_review.csv")), len(reviews))

    def test_source_mismatch_fails_without_publishing_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            snapshot = make_snapshot(root, [source_row()])
            changed = source_row(smiles="CCN")
            write_csv(snapshot / "dsstox_review.csv", REVIEW_FIELDS, [changed])
            with self.assertRaisesRegex(ValueError, "differs"):
                run(snapshot, root / "reviewed")
            self.assertFalse((root / "reviewed").exists())


if __name__ == "__main__":
    unittest.main()
