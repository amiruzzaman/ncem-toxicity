# ToxCast v4.3 CSV import

This import adds EPA's August 2025 invitrodb v4.3 processed release to the reviewed Tox21 + DSSTox catalog. The original snapshots and README stay unchanged. ChEMBL and an agentic retrieval engine remain future work.

## Rebuild

Use Python 3.11 or newer with `openpyxl==3.1.5` (`python3 -m pip install -r requirements/toxcast.txt`). Allow roughly 25 GB of free space for the archived inputs, a new snapshot, and temporary staging; a simultaneous verification rebuild needs additional space. Run from the repository root:

```bash
python3 download_toxcast.py --source-dir sources/toxcast/invitrodb_v4_3
python3 import_toxcast.py \
  --base-snapshot tox21_dsstox_complete \
  --key-review tox21_key_review \
  --source-dir sources/toxcast/invitrodb_v4_3 \
  --output tox21_dsstox_toxcast_v4_3
python3 verify_toxcast.py \
  --snapshot tox21_dsstox_toxcast_v4_3 \
  --reviewed-master tox21_key_review/master_chemicals_reviewed.csv \
  --report reports/toxcast/validation-recheck.json
python3 -m unittest -v tests.test_import_toxcast
```

The original full-release verification certificate is retained at `reports/toxcast/validation.json`. Use a separate report filename for new verification runs to preserve that historical evidence.

The downloader uses `toxcast_sources.json`, containing EPA's public file URLs, sizes, SHA-256 fingerprints, and release DOI. It resumes `.part` files where HTTP ranges are supported, restarts safely if the server ignores the range, and verifies existing files before reuse. An interrupted download can be resumed by running the same command. A hash mismatch requires investigating the source or damaged local file rather than silently changing the lock.

The importer rejects existing output directories. Use a fresh name for every rebuild. It verifies baseline and key-review manifests, source hashes, schemas, processing IDs, and assay references. Inputs are checked again before the staged directory is atomically published. Temporary SQLite is a build index; the durable data remains CSV plus JSON manifests.

For a full determinism check, rerun the import with `--output tox21_dsstox_toxcast_v4_3_verify`, then run the verifier above with `--compare tox21_dsstox_toxcast_v4_3_verify`. It compares complete manifests and hashes every corresponding output file. The optional validation report records these checks and implementation hashes outside the immutable snapshots. Generated snapshots and downloads are Git-ignored.

## Outputs and joins

| Output | Grain and purpose |
| --- | --- |
| `master_chemicals.csv` | Original reviewed chemical rows followed by new substance entries; original 24 columns plus three ToxCast coverage counts. |
| `toxcast_mc_observations.csv` | Every released winning-model/flag source row; `m5id` identifies a distinct processed MC result. |
| `toxcast_sc_observations.csv` | Every SC summary row; `s2id` identifies a distinct processed SC result. |
| `toxcast_sc_replicates.csv` | Original SC1 replicate/control rows; these do not multiply processed-result coverage. |
| `toxcast_chemicals.csv` | Distinct supplied chemical identifiers/names/CAS combinations with the resolved internal ID or review status. |
| `toxcast_assays.csv`, `assays/` | Complete endpoint annotations and separate assay/component/endpoint tables. Join endpoints through `aeid`. |
| `toxcast_endpoints.csv`, `toxcast_assay_review.csv` | Complete referenced endpoint catalog and an explicit review of endpoints absent from the released annotations. |
| `toxcast_targets.csv` | One-to-many endpoint-to-gene/AOP/key-event annotations. A join can produce several target rows per result. |
| `methods/`, `quality/` | Processing method assignments/definitions and released assay-performance information. |
| `toxcast_analytical_qc.csv` | Original sample/substance QC, including supplied cautions and predicted properties. |
| `toxcast_cytotoxicity.csv` | Cytotoxicity summaries with an added default-versus-estimate designation. |
| `derived/` | Endocrine reference/model/QSAR scores and administered-equivalent dose estimates, excluded from experimental coverage. |
| `toxcast_review.csv` | Unresolved identities, missing sample identifiers, and invalid calls, with source locators. |
| `coverage_summary.csv`, `endpoint_coverage.csv` | Actual import gains and endpoint coverage. |
| `sample_disagreements.csv` | Chemicals with both active and inactive samples within the same endpoint and screening type. |
| `MEETING_SUMMARY.md`, `trace_example.json` | Meeting brief and one compound/result/assay/target trace with source evidence. |
| `toxcast_import.json` | Source locks, baseline fingerprints, archive inventory, source/output row accounting, output hashes, and validation results. |

Original DSSTox/Tox21 companion files are copied byte for byte. Prior review evidence is retained under `key_review/`; its manifests describe the historical baseline, while `toxcast_import.json` describes this new snapshot.

Every exported ToxCast table keeps its native source columns and appends chemical ID, release, source file/sheet/logical row, fingerprint, and evidence classification. For archive members, `source_sha256` fingerprints the entire pinned archive and `source_file` identifies the member listed in the manifest. Standalone metadata uses its own file hash. Original files labeled AUG2024 inside the archive belong to the August 2025 release; the independently provided AUG2025 annotation files are verified byte-identical before use.

## Matching and interpretation

Chemical matching uses DTXSID mappings from **all** linked DSSTox records. Exactly one existing chemical ID is required. New valid DTXSIDs receive stable substance IDs using the existing `researchchem:dsstox:substance:` UUID convention. New master rows are marked `toxcast_substance_only`, with structures and original-source metadata blank; names and DTXSIDs are available through `toxcast_chemicals.csv`. Names/CAS are never merge keys. Ambiguous or missing DTXSIDs remain unresolved, including control wells that have no curated substance assignment. No structures are reconstructed from QSAR-ready representations or model supplements.

`activity_status` is added without replacing native `hitc`/`hitcall`. MC calls from 0.90 through 1 are active; values from 0 through below 0.90 are inactive. Negative-direction calls, missing calls, unfit models (`none`/`nofit`), and invalid calls are distinct. SC calls use their original binary definition. The signed `hitc` takes precedence over the additional model `hitcall` field when both are provided.

`potency_unit` preserves the actual source concentration unit. Most are uM, but the release also contains CF (concentration factor), and some units are missing. Values are not automatically converted, and missing units remain unknown. Native logarithmic concentration columns are retained as logarithmic quantities. A cross-assay potency comparison requires compatible units, fit/model context, and assay conditions.

MC6 flags, assay usability, analytical QC, and disagreements remain visible; no automatic quality filter removes data. A finite active/inactive call is considered *interpretable for coverage*, not a declaration that the study passes every quality criterion. Cytotoxicity defaults (1000 uM with inadequate supporting burst hits) are labeled separately from assay-based estimates. Derived model estimates never count as new experimental coverage.

The processed release contains result endpoints absent from its annotation workbook. Their native assay IDs and result rows remain intact. The complete endpoint catalog marks them `missing_released_annotation`, and the assay review file supplies source examples and affected result counts. Use a left join from observations to annotations so these records stay visible. No target, protocol, or biological description is inferred to fill that gap.

The three master columns count distinct `m5id`, distinct `s2id`, and distinct `aeid` across MC/SC results linked to each chemical. Coverage calculations do not count replicate rows or target-join multiplicity. Missing/unfit/negative/invalid calls count as source presence but not interpretable coverage. ToxCast includes Tox21-origin assays: chemical overlap does not establish independent observations, and the original 12 Tox21 endpoint labels are never overwritten.

Source: [EPA ToxCast downloads](https://www.epa.gov/comptox-tools/exploring-toxcast-data), [release DOI](https://doi.org/10.23645/epacomptox.6062623.v14), and [EPA processing documentation](https://usepa.r-universe.dev/tcpl/doc/Introduction_Appendices.html).


## Repository and Git policy

Keep CLI scripts and `toxcast_sources.json` at the repository root so existing commands and historical input paths remain valid. Tests live in `tests/`; the ToxCast dependency pin is `requirements/toxcast.txt`. The README documents the original pilot, and `CHANGELOG.md` records subsequent work.

Commit code, tests, dependency pins, the source lock, this guide, the changelog, and the curated `reports/toxcast/validation.json` certificate. The small existing `tox21.csv` and `tox21_summary.txt` remain tracked reference files. Do not commit downloaded releases, generated snapshots, temporary indexes, local meeting/coverage/trace report copies, or legacy exports.

The current master is `tox21_dsstox_toxcast_v4_3/master_chemicals.csv`. Full reports are generated alongside it. Optional copies under `reports/toxcast/` are local sharing artifacts; only `validation.json` is allowlisted in Git. Recheck reports use a different filename to preserve the original certificate. Do not edit published snapshot contents or move their inputs: manifests record hashes and historical paths.

Run the DSSTox suite with `python3 -m unittest -v tests.test_review_dsstox_keys` in an RDKit-equipped environment. With both RDKit and openpyxl installed, `python3 -m unittest discover -s tests -v` runs all 24 tests.
