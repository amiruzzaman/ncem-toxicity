"""Header validation and label preservation for the Tox21 builder."""

import csv
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from BuiltMasterTox21 import ENDPOINTS, build


class HeaderTests(unittest.TestCase):
    def setUp(self):
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        root = Path(directory.name)
        self.source = root / "tox21.csv"
        self.output = root / "snapshot"
        self.fields = ("mol_id", "smiles", *ENDPOINTS)

    def write_source(self, fields, rows=()):
        with self.source.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(fields)
            writer.writerows(rows)

    def assert_duplicate_error(self, names):
        with self.assertRaises(ValueError) as caught:
            build(self.source, self.output)
        self.assertEqual(str(caught.exception), f"Duplicate Tox21 columns: {names}")
        self.assertFalse(self.output.exists())

    def test_conflicting_duplicate_assay_headers_are_rejected(self):
        self.write_source(
            (*self.fields, "NR-AR"),
            [("sample", "CCO", "1", *([""] * 11), "0")],
        )
        self.assert_duplicate_error("NR-AR")

    def test_duplicate_identifier_and_extra_headers_are_rejected(self):
        for name in ("mol_id", "smiles", "extra"):
            with self.subTest(name=name):
                self.output = self.output.parent / f"snapshot-{name}"
                fields = self.fields if name in self.fields else (*self.fields, name)
                self.write_source((*fields, name))
                self.assert_duplicate_error(name)

    def test_duplicate_names_are_unique_and_sorted(self):
        self.write_source(
            (*self.fields, "z_extra", "a_extra", "z_extra", "mol_id", "a_extra", "mol_id")
        )
        self.assert_duplicate_error("a_extra, mol_id, z_extra")

    def test_duplicate_error_precedes_missing_required_column(self):
        fields = tuple(name for name in self.fields if name != "smiles")
        self.write_source((*fields, "NR-AR"))
        self.assert_duplicate_error("NR-AR")

    def test_missing_required_column_is_still_rejected(self):
        self.write_source(tuple(name for name in self.fields if name != "smiles"))
        with self.assertRaisesRegex(ValueError, "Input is missing required Tox21 columns"):
            build(self.source, self.output)
        self.assertFalse(self.output.exists())

    def test_valid_input_preserves_labels_and_exact_header_names(self):
        labels = {name: "" for name in ENDPOINTS}
        labels.update({"NR-AR": "1", "SR-p53": "0"})
        self.write_source(
            (*self.fields, "nr-ar", " NR-AR "),
            [("sample", "CCO", *(labels[name] for name in ENDPOINTS), "0", "0")],
        )
        build(self.source, self.output)
        with (self.output / "master_chemicals.csv").open(newline="", encoding="utf-8") as stream:
            master = list(csv.DictReader(stream))
        self.assertEqual(len(master), 1)
        self.assertEqual({name: master[0][name] for name in ENDPOINTS}, labels)
        with (self.output / "tox21_observations.csv").open(newline="", encoding="utf-8") as stream:
            observations = list(csv.DictReader(stream))
        self.assertEqual(
            [(row["endpoint"], row["result"]) for row in observations],
            [("NR-AR", "1"), ("SR-p53", "0")],
        )


if __name__ == "__main__":
    unittest.main()
