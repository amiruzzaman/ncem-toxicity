# ncem-toxicity

The current pilot combines Tox21 assay results with EPA DSSTox chemical records in CSV files. Future work will bring in ChEMBL and TOXRIC and build an AI retrieval tool that follows directed acyclic graph (DAG) links between compounds, assays, targets, and toxicity endpoints.

## What is in MasterChemCSV?

The current master is `tox21_dsstox_complete/master_chemicals.csv`. Each row represents a chemical identity or a provisional entry awaiting structure review. An assay is a lab test that measures a particular response.

| Columns | What they tell you |
| --- | --- |
| `chemical_id`, `standard_inchikey`, `representative_smiles`, `identity_status` | The project's internal ID, the standard structure identifier when available, a text representation of the structure, and how the identity was recorded. |
| `dsstox_*` | DSSTox identifiers, preferred name, CAS number, molecular formula, and mass fields where available. A nonstandard DSSTox InChIKey stays in `dsstox_inchikey`. |
| Twelve Tox21 assay columns, such as `NR-ER` | The result for that chemical in each assay. |

A value of `1` means active in the specified assay. A value of `0` means inactive in that assay. A blank means no result is available. A value of `conflict` means the source observations disagree. Activity is specific to the assay and its conditions.

For example, chemical ID `94c966d3-e7b1-56b4-9e8c-757804e199c7` has `conflict` in `NR-ER`. Filter `tox21_observations.csv` by that `chemical_id` and endpoint `NR-ER`. Source record `TOX4994` reports `1`, and `TOX6178` reports `0`. The two records support reporting a disagreement for that assay.

## Data sources and matching keys

| Input | What it contributes | How records are matched |
| --- | --- | --- |
| `tox21.csv`, included in this repository | Chemical structures and results for 12 assays | RDKit generates a full Standard InChIKey from each SMILES string. Rows with the same full key share a `chemical_id`. The original `mol_id` stays in the observations as the source record ID. |
| EPA `DSSToxCCDdump.csv` | Chemical identifiers, names, structures, and properties | An exact match on the full supplied InChIKey reuses the Tox21 `chemical_id`. Records with new keys add master rows. |

An InChIKey is a short identifier derived from a chemical structure. Names and CAS numbers are kept as metadata. For a DSSTox record with a missing or unusable key, the importer checks whether its DSSTox substance ID (`DTXSID`) points to exactly one usable key in the same export. All records with missing or malformed keys go into the review queue. An unresolved record receives a provisional ID. Use `chemical_id` to look up source records and observations from the master.

Use the full `DSSToxCCDdump.csv` inside `DSSTox_CCD_dump_12092025_CSVs.zip`. The numbered CSVs are partitions of that same file. If the historical archive is unavailable through EPA, ask the lab for its archived copy. [Tox21 Data Challenge](https://tripod.nih.gov/tox21/challenge/data.jsp) · [EPA computational toxicology downloads](https://www.epa.gov/comptox-tools/downloadable-computational-toxicology-data)

| Input | SHA-256 fingerprint |
| --- | --- |
| `tox21.csv` | `7d7e7facd853a63e79ddce4e9c3fcb7a0d83a1c300b603031c0f2c64fbe77761` |
| `DSSToxCCDdump.csv` | `e69f56b35ce9d626810c6df2b5c79c48d0c3bb50ccde399da30687c6dd6d6a1b` |

## Set up and build

Use **Python 3.11 or later**. The Tox21 builder needs RDKit in the Python environment used to run it. The verified build used **RDKit 2026.09.1**. Check which version your interpreter loads with the command below. If the import fails, install RDKit in that environment. [RDKit installation guide](https://www.rdkit.org/docs/Install.html)

```bash
python3 -c "from rdkit import rdBase; print(rdBase.rdkitVersion)"
```

Place the extracted `DSSToxCCDdump.csv` alongside the scripts. Allow several GB of free space for the input, generated files, and temporary staging. Run these commands from the repository root. The summary command is optional. Both builders require fresh output directory names, so use unused names if you already have a snapshot.

```bash
python3 summary_tox21.py tox21.csv
python3 BuiltMasterTox21.py tox21.csv tox21_initial
python3 import_dsstox.py tox21_initial DSSToxCCDdump.csv tox21_dsstox_complete
```

The first build creates the Tox21-only master and observation CSVs. The DSSTox import uses those files to create the combined snapshot. SQLite serves as temporary staging during the import. The published data files are CSVs, with a JSON manifest recording the import.

## Combined output files

| File in `tox21_dsstox_complete/` | What it contains |
| --- | --- |
| `master_chemicals.csv` | The main chemical catalog with selected DSSTox properties and Tox21 assay columns. |
| `tox21_observations.csv` | Individual Tox21 results with the source record ID, original SMILES, and source fingerprint. |
| `dsstox_records.csv` | All linked DSSTox records, including all 13 original fields, source row numbers, and the source fingerprint. |
| `dsstox_review.csv` | DSSTox records with missing or malformed keys, with a review reason and a link to the master. |
| `dsstox_import.json` | The DSSTox input fingerprint, import policy version, and counts for checking the build. |

## Expected counts and future imports

These counts apply to the inputs above and DSSTox import policy version **2**.

| Check | Expected count |
| --- | ---: |
| Original Tox21 input records | 7,831 |
| Tox21-only master rows | 7,829 |
| Tox21 observations | 77,946 |
| Linked DSSTox input records | 1,246,399 |
| Existing Tox21 identities matched by exact InChIKey | 6,006 |
| Combined master rows | **1,248,222** |
| Provisional DSSTox entries included in the combined master | 93,332 |

DSSTox expands the chemical catalog. Tox21 supplies all current assay results, and only 7,829 master rows have any Tox21 result. The Tox21 input summary breaks down active, inactive, and missing values by endpoint. [Tox21 summary](tox21_summary.txt)

For a new DSSTox release, start from the Tox21-only snapshot, build into a fresh directory, and check the new manifest against the new input. For other datasets, document the release and matching rule and keep the original observations with their source identifiers. Use a CSV-aware tool, such as Python's built-in `csv` module, to filter the full master into a smaller CSV for browsing.
