# Changes

## 2026-11-06

- Added the pinned EPA ToxCast v4.3 downloader, streaming CSV/Excel importer, and independent snapshot verifier. Downloads resume safely; verified snapshots publish atomically.
- Preserved native measurements, units, quality flags, samples, source provenance, and the reviewed Tox21 + DSSTox identities. Added three ToxCast coverage columns to the master.
- Kept single-concentration replicates and modeled supplements separate from experimental coverage. Retained unresolved identities and 17 missing endpoint annotations for review.
- Imported **4,060,981 processed assay results**. **4,519 existing chemicals** gained interpretable results, increasing original-catalog coverage from **0.6272% to 0.9892%**. The expanded catalog has **1,252,222 entries**; its 4,000 new substance entries add no interpretable experimental coverage.
- Verified a byte-identical full rebuild, independent coverage reconciliation, and all **24 regression tests**. The [validation certificate](reports/toxcast/validation.json) records the implementation hashes and verification results.
- Grouped tests and dependency configuration into `tests/` and `requirements/`. Consolidated rebuild, interpretation, and repository guidance in the [ToxCast guide](docs/toxcast-import.md).
- Kept downloaded releases, generated snapshots, local report copies, and legacy exports out of Git. Source locks, dependency pins, code, tests, and the curated validation certificate remain versioned.

ChEMBL ingestion and DAG-based agentic retrieval remain future work. The README remains the original Tox21 + DSSTox pilot guide.
