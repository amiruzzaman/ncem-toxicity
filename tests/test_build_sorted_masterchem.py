import csv
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import build_sorted_masterchem as builder
from download_toxcast import fingerprint
from import_toxcast import SUMMARY_FIELDS


MASTER_FIELDS = (
    "chemical_id", "standard_inchikey", "representative_smiles", "identity_status",
    *builder.TOX21_ASSAYS,
    "dsstox_inchikey", "dsstox_preferred_name", "dsstox_casrn", "dsstox_dtxsid",
    "dsstox_dtxcid", "dsstox_molecular_formula", "dsstox_average_mass",
    "dsstox_monoisotopic_mass", *SUMMARY_FIELDS,
)
OBSERVATION_FIELDS = (
    "chemical_id", "aeid", "result_id", "hitc", "hitcall", "modl", "activity_status",
)


def chemical(chemical_id, usable, **metadata):
    row = dict.fromkeys(MASTER_FIELDS, "")
    for index, field in enumerate(builder.TOX21_ASSAYS[:usable]):
        row[field] = str(index % 2)
    if usable < len(builder.TOX21_ASSAYS):
        row[builder.TOX21_ASSAYS[usable]] = "conflict"
    row.update(chemical_id=chemical_id, **metadata)
    return row


def observation(chemical_id, endpoint, call="1", status="active", model="hill"):
    return {
        "chemical_id": chemical_id,
        "aeid": str(endpoint),
        "result_id": f"{chemical_id or 'unlinked'}-{endpoint}",
        "hitc": call,
        "hitcall": "1",
        "modl": model,
        "activity_status": status,
    }


class SortedMasterChemTests(unittest.TestCase):
    def setUp(self):
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        root = Path(directory.name)
        self.snapshot = root / "snapshot"
        self.snapshot.mkdir()
        self.output = root / "SortedMasterChem.csv"

    def write_fixture(self, master, mc=(), sc=()):
        tables = {
            "master_chemicals.csv": (MASTER_FIELDS, master),
            "toxcast_mc_observations.csv": (OBSERVATION_FIELDS, mc),
            "toxcast_sc_observations.csv": (OBSERVATION_FIELDS, sc),
        }
        for name, (fields, rows) in tables.items():
            with (self.snapshot / name).open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=fields)
                writer.writeheader()
                writer.writerows(rows)
        manifest = {"output_sha256": {
            name: fingerprint(self.snapshot / name) for name in tables
        }}
        (self.snapshot / "toxcast_import.json").write_text(
            json.dumps(manifest), encoding="utf-8",
        )

    def test_either_source_qualifies_and_export_preserves_master_values(self):
        tox21 = chemical("c-tox21", 6)
        toxcast = chemical("b-toxcast", 0, dsstox_preferred_name='Name, with "quotes"\nand a line')
        both = chemical("z-both", 6, dsstox_dtxsid="DTXSID1000001")
        high = chemical("a-high", 7)
        fails = chemical("d-fails", 5)
        observations = [
            observation(row["chemical_id"], endpoint, "0", "inactive")
            for row, count in ((toxcast, 6), (both, 6), (fails, 5))
            for endpoint in range(1, count + 1)
        ]
        self.write_fixture([tox21, toxcast, both, fails, high], sc=observations)

        summary = builder.build(self.snapshot, self.output)

        self.assertEqual(summary, {
            "master_rows": 5, "exported_rows": 4,
            "tox21_only": 2, "toxcast_only": 1, "both": 1,
        })
        with self.output.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            self.assertEqual(tuple(reader.fieldnames), (
                *MASTER_FIELDS, builder.TOX21_COUNT, builder.TOXCAST_COUNT,
            ))
            self.assertEqual(list(reader), [
                {**high, builder.TOX21_COUNT: "7", builder.TOXCAST_COUNT: "0"},
                {**toxcast, builder.TOX21_COUNT: "0", builder.TOXCAST_COUNT: "6"},
                {**tox21, builder.TOX21_COUNT: "6", builder.TOXCAST_COUNT: "0"},
                {**both, builder.TOX21_COUNT: "6", builder.TOXCAST_COUNT: "6"},
            ])

    def test_counts_distinct_agreeing_endpoints_across_screens(self):
        mc = [
            observation("c", 1, "0.90"),
            observation("c", 1, "0.90"),
            observation("c", 1, "-0.8", "negative_direction"),
            observation("c", 1, "0", "unfit", model="none"),
            observation("c", 2, "0.899999", "inactive"),
            observation("c", 3),
            observation("c", 4),
            observation("c", 4, "0", "inactive"),
            observation("c", 5, "NA", "missing"),
            observation("", 6),
        ]
        sc = [
            observation("c", "1.0"),
            observation("c", 2, "0", "inactive"),
            observation("c", 2, "0.9", "invalid"),
            observation("c", 3, "0", "inactive"),
            observation("c", 7, "0.9", "invalid"),
        ]
        self.write_fixture([chemical("c", 0)], mc=mc, sc=sc)

        self.assertEqual(builder.toxcast_counts(self.snapshot), {"c": 2})

    def test_rechecks_native_status_for_an_unlinked_observation(self):
        self.write_fixture([chemical("c", 6)], mc=[observation("", 1, "0.89")])
        with self.assertRaisesRegex(ValueError, r"result unlinked-1:.*expected 'inactive'"):
            builder.build(self.snapshot, self.output)
        self.assertFalse(self.output.exists())

    def test_rejects_a_hash_mismatch_before_parsing(self):
        self.write_fixture([chemical("c", 6)])
        (self.snapshot / "master_chemicals.csv").write_text("changed", encoding="utf-8")

        with self.assertRaisesRegex(ValueError, "master_chemicals.csv: SHA-256 mismatch before reading"):
            builder.build(self.snapshot, self.output)
        self.assertFalse(self.output.exists())

    def test_rejects_an_input_changed_after_reading(self):
        self.write_fixture([chemical("c", 6)])
        validate = builder.validate_staged_csv

        def change_input(*args):
            validate(*args)
            path = self.snapshot / "toxcast_sc_observations.csv"
            with path.open("a", encoding="utf-8") as stream:
                stream.write("new bytes\n")
            manifest_path = self.snapshot / "toxcast_import.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["output_sha256"][path.name] = fingerprint(path)
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

        with patch.object(builder, "validate_staged_csv", side_effect=change_input):
            with self.assertRaisesRegex(ValueError, "toxcast_sc_observations.csv: SHA-256 mismatch after reading"):
                builder.build(self.snapshot, self.output)
        self.assertFalse(self.output.exists())

    def test_rejects_existing_output_and_paths_inside_the_snapshot(self):
        self.write_fixture([chemical("c", 6)])
        self.output.write_text("existing content", encoding="utf-8")
        with self.assertRaises(FileExistsError):
            builder.build(self.snapshot, self.output)
        self.assertEqual(self.output.read_text(encoding="utf-8"), "existing content")
        with self.assertRaisesRegex(ValueError, "outside the input snapshot"):
            builder.build(self.snapshot, self.snapshot / self.output.name)

    def test_rejects_duplicate_master_ids(self):
        self.write_fixture([chemical("c", 6), chemical("c", 0)])
        with self.assertRaisesRegex(ValueError, "duplicate chemical_id"):
            builder.build(self.snapshot, self.output)
        self.assertFalse(self.output.exists())


if __name__ == "__main__":
    unittest.main()
