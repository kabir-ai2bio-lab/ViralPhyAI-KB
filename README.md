# ViralPhyAI-KB

Data processing pipeline for *ViralAI-KB: AI-Ready Biophysical Knowledge
Base for Viral Proteins*. It builds a deduplicated, taxonomically
verified dataset of viral protein molecular dynamics (MD) data from four
public repositories (MDRepo, MDDB, mdCATH, DynoDB), then computes
per-residue and per-residue-pair dynamics features from each selected
simulation's trajectory.

This folder has the code for retrieval, verification, harmonization,
deduplication, final representative selection, and feature extraction,
plus the small dataset-level result files produced for the paper. It
does not include plotting code, exploratory scripts, or the full
per-protein feature files (2.5GB, hosted on Zenodo instead). See the
paper's Methods section for the full procedure.

## Pipeline overview

```
analysis/{mdrepo,mddb,mdcath,dynodb}/collect_metadata.py   (Stage 1, per source)
                    │
                    ▼
analysis/cross_database/dedup.py                           (Stage 2)
                    │
                    ▼
analysis/cross_database/final_selection.py                 (Stage 3)
                    │
                    ▼
analysis/cross_database/download_merged_files.py            (Stage 4)
analysis/dynodb/download_trajectories.py
                    │
                    ▼
features_extraction/convert_mdcath.py                       (Stage 5)
features_extraction/convert_dynodb.py
                    │
                    ▼
features_extraction/run_feature_extraction.py                (Stage 6)
  (orchestrates GetContacts + MD_features.py per representative)
```

Each stage is resumable and safe to re-run: later stages read files
written by earlier ones rather than recomputing in memory, and the
feature-extraction stage skips any record whose output already exists.

## Included results

The outputs of Stages 2-3 are included directly, as actually produced
for the paper (verified to match its reported counts: 3,315 kept records
spanning 1,535 distinct sequences):

- `analysis/cross_database/combined_records.json` (3.3MB): every record
  considered during deduplication, each tagged `status: "kept"` or
  `"removed"` (with removed records annotated with what they
  duplicated). This is the Merged Dataset referenced throughout the
  paper's Results section.
- `analysis/cross_database/cross_database_dedup_report.json`: summary
  counts from the deduplication run (per-source kept/removed, excluded
  heteromeric complexes).
- `features_extraction/final_representative_ids.json`: the 1,535 final
  representatives (one per distinct sequence) by source dataset.
- `analysis/cross_database/viral_protein_cross_database.xlsx` (1.2MB):
  a human-readable multi-sheet workbook (per-source detail sheets, an
  overview, and the list of removed duplicates with their cause). Note:
  its generating script (`build_excel.py`) is not included in this
  release, since it is reporting code rather than data-processing code
  and pulls in extra dependencies (`openpyxl`, live UniProt/RCSB/NCBI
  family-resolution calls) not needed anywhere else in this pipeline,
  but the output file itself is included for direct inspection.

If you just want to inspect the dataset or jump straight to feature
extraction without re-running retrieval, start from these files and skip
to Stage 4. If you re-run Stages 1-3 from scratch, expect somewhat
different counts than the files above, since the source databases grow
over time (the paper's numbers reflect data as retrieved on September 1,
2026) and Stage 1 performs fresh live-network verification.

The full per-protein feature-extraction output (`features/{source}/`,
1,516 x 2 files, 2.5GB total) is not included in this code release. It
is hosted on Zenodo instead; see "Feature extraction output" below.

## Requirements

- Python 3.10 (required specifically for Stage 6: `vmd-python`, used by
  GetContacts, only ships Python 3.10 builds)
- `requests`, `numpy`, `pandas`, `mdtraj`, `h5py`
- For Stage 6 only: a working `vmd-python` environment and a local clone
  of [GetContacts](https://getcontacts.github.io/) at
  `features_extraction/getcontacts/` (not included here; it is
  third-party software, cited in the paper as `\cite{getcontacts}`):

  ```bash
  conda create -n viralmd-features python=3.10
  conda activate viralmd-features
  conda install -c conda-forge vmd-python mdtraj h5py
  pip install requests numpy pandas
  git clone https://github.com/getcontacts/getcontacts.git features_extraction/getcontacts
  ```

Stages 1-5 only need `requests`/`numpy`/`pandas` and live network access
(no GPU, no `vmd-python`); Stage 6 is the only stage that needs the
`viralmd-features` environment above.

## Stage 1: Per-source retrieval, verification, harmonization

```bash
cd analysis/mdrepo  && python3 collect_metadata.py
cd analysis/mddb    && python3 collect_metadata.py
cd analysis/mdcath  && python3 collect_metadata.py
cd analysis/dynodb  && python3 collect_metadata.py
```

Each script retrieves candidate records from its source's public API (or,
for mdCATH, CATH's static domain-boundary files), confirms viral origin
via UniProt / RCSB PDB / NCBI Taxonomy lookups, excludes heteromeric
complexes (records whose chains are not all identical; this resource
targets single-chain and homo-oligomeric proteins only), and writes one
JSON file per confirmed-viral, single-chain-or-homo-oligomeric record to
`simulations/` in its own directory. This is a live-network stage and
will take a while; it is safe to re-run (it just re-verifies).

Numbers reported in the paper as of September 1, 2026: 83,465 records
scanned across all four sources, 6,623 confirmed viral, 798 excluded as
heteromeric complexes. Source databases grow over time, so a fresh run
will likely scan more records and should be expected to produce somewhat
different counts.

## Stage 2: Cross-database deduplication

```bash
cd analysis/cross_database
python3 dedup.py
```

Reads every `../{source}/simulations/*.json` file, harmonizes each
source's record into one common schema (sequence, UniProt ID, organism,
temperature, force field, water model, ligand state), and removes
redundancy in two passes: exact sequence + simulation-setup matching
within each source, then source-reported-provenance-based resolution of
mdCATH/MDDB/DynoDB cross-database mirrors. Writes `combined_records.json`
(every record, each tagged `status: "kept"` or `"removed"`, with removed
records annotated with what they duplicated) and
`cross_database_dedup_report.json` (summary counts).

> The versions of both output files used in the paper are already
> included in this directory (see "Included results" above); running
> this script yourself will overwrite them with a fresh run's results.

## Stage 3: Final representative selection

```bash
cd analysis/cross_database
python3 final_selection.py
```

Reads `combined_records.json`'s kept records, groups them by exact
sequence, and picks one representative per distinct sequence via a soft
filter cascade: prefer apo (unbound) → prefer TIP3P water model → prefer
temperature closest to 310 K (human physiological temperature) among
candidates with a known value → first encountered as a final tie-break.
Writes `../../features_extraction/final_representative_ids.json`,
`{"mdrepo": [...], "mddb": [...], "mdcath": [...], "dynodb": [...]}`,
the worklist every later stage reads.

> This file is already included at `features_extraction/final_representative_ids.json`
> (see "Included results" above); running this script yourself will
> overwrite it with a fresh run's results.

## Stage 4: Fetching structure and trajectory files

```bash
cd analysis/cross_database
python3 download_merged_files.py
```

Downloads the real structure + trajectory files for every **kept**
record (not just the final one-per-sequence winners, for flexibility) into
`merged_downloads/{source}/{record_id}_structure.*` +
`{record_id}_trajectory.*`.

DynoDB is a special case: its trajectory ensembles were released by its
authors after initial metadata collection, so they are fetched
separately:

```bash
cd analysis/dynodb
python3 download_trajectories.py
```

This writes `{pdb_id}_frame.pdb` files into `trajectories/`, for DynoDB
records that won a final-selection slot only. Move/rename each one into
the same `merged_downloads/dynodb/` directory Stage 4's first command
already populated, under the filename Stage 5 expects:

```bash
cd analysis/dynodb
for f in trajectories/*_frame.pdb; do
  pdb_id=$(basename "$f" _frame.pdb)
  cp "$f" "../cross_database/merged_downloads/dynodb/${pdb_id}_trajectory.pdb"
done
```

## Stage 5: Normalizing mdCATH and DynoDB

mdrepo/mddb/dynodb files are already structure+trajectory pairs; mdCATH
(`.h5`) and DynoDB (mismatched structure/trajectory atom counts, native
PDB residue numbering) need conversion into the same
`{record_id}_structure.pdb` + `{record_id}_trajectory.xtc` layout first:

```bash
cd features_extraction
python3 convert_mdcath.py --all
python3 convert_dynodb.py --all
```

## Stage 6: Feature extraction

```bash
cd features_extraction
conda activate viralmd-features
python3 run_feature_extraction.py --cores 1
```

For every final representative (`final_representative_ids.json`), this:

1. Renumbers residues to a contiguous `1..N` range (several sources keep
   native/non-contiguous PDB numbering, which the next step requires).
2. Runs GetContacts to detect non-covalent contacts per frame across 9
   interaction types (van der Waals, 3 hydrogen-bond subtypes,
   hydrophobic, salt bridge, π-cation, π-stacking, T-stacking).
3. Runs `MD_features.py` (adapted from Hou et al., SeqDance,
   `\cite{seqdance2026}`), which computes per-residue secondary
   structure (DSSP), solvent accessible surface area, backbone/side-chain
   dihedral angle distributions, and RMSF via MDTraj
   (`\cite{mdtraj2015}`); computes Cα positional correlation from the
   trajectory's coordinate covariance; and assembles the GetContacts
   contact map with the correlation matrix into one combined array.

Writes, per representative, into `features/{source}/`:
- `{record_id}_res_feature.csv`: one row per residue
- `{record_id}_pair_feature.npz`: one `N×N×10` array (9 GetContacts
  channels + 1 correlation channel)

Useful flags: `--sources mdrepo,mddb,...` to restrict to specific
sources, `--limit N` for a pilot run, `--skip-failed` to skip records
already listed in `run_feature_extraction_failures.json` instead of
retrying them, `--keep-contacts` to retain the (large) intermediate
`*_contacts.tsv` files. `--cores 1` is the tested-safe default; GetContacts'
own multiprocessing has been observed to deadlock with `--cores > 1` in
at least one sandboxed environment; raise it if you confirm it's stable
in yours.

The run is resumable (skips any record whose two output files already
exist) and never aborts the whole batch on one record's failure; it
logs the failure and continues. In the run behind the paper, 1,516 of
1,535 (98.8%) final representatives succeeded; the remaining 19 (all
MDDB records) failed for two reasons: 12 from a residue-indexing edge
case triggered by non-canonical residues (most often a bound ligand
assigned a standard amino-acid name in the source structure), and 7 from
memory limits on the largest, multi-thousand-residue structures in the
dataset.

> `--keep-contacts` leaves behind large intermediate `*_contacts.tsv`
> files (and a failed run also leaves behind an intermediate
> `*_renumbered.pdb`, since cleanup only runs after a record fully
> succeeds); neither is part of the documented output format above, and
> neither should be published. Clean them out of `features/{source}/`
> before archiving if you used that flag.

## Feature extraction output

The full feature set described above (2.5GB: `{record_id}_res_feature.csv`
+ `{record_id}_pair_feature.npz` for all 1,516 successfully processed
representatives) is hosted on Zenodo rather than included in this code
repository:

**DOI: [10.5281/zenodo.23590905](https://doi.org/10.5281/zenodo.23590905)**

It is split into eight archives, grouped by source database (`mdrepo` and
`mddb` further split into 3 roughly-equal-byte-size parts each, to keep
individual uploads well under Zenodo's production gateway's apparent
large-request ceiling; see `zenodo_upload/upload_to_zenodo.py`'s
comments). All parts for a source together cover every representative
for that source, with no overlap. Counts match
`features_extraction/final_representative_ids.json`'s per-source totals:

| Archive | Representatives | Size |
|---|---|---|
| `mdrepo_features_part1.tar.gz` | 276 | 0.40GB |
| `mdrepo_features_part2.tar.gz` | 276 | 0.40GB |
| `mdrepo_features_part3.tar.gz` | 277 | 0.40GB |
| `mddb_features_part1.tar.gz` | 117 | 0.48GB |
| `mddb_features_part2.tar.gz` | 118 | 0.47GB |
| `mddb_features_part3.tar.gz` | 117 | 0.48GB |
| `mdcath_features.tar.gz` | 11 | 2.4MB |
| `dynodb_features.tar.gz` | 324 | 14MB |

SHA-256 checksums for each archive are listed in `checksums.sha256`
alongside them on Zenodo. To reproduce this output yourself instead of
downloading it, run Stage 6 above.

### Publishing an update

If you re-run the pipeline and need to publish a new version of this
dataset, `zenodo_upload/upload_to_zenodo.py` (not included in this code
repository, it lives alongside the archives since it is upload tooling
rather than data-processing code) automates the Zenodo side: it creates a
draft deposition, uploads the archives via Zenodo's bucket API with
checksum verification and automatic retry on transient errors, and sets
metadata, but never publishes automatically, since publishing a Zenodo
record is irreversible (the DOI and uploaded files become permanent; a
later fix requires a new version). Requires a Zenodo personal access
token in the `ZENODO_TOKEN`
environment variable (never hardcoded or passed as a CLI argument).
Test against `sandbox.zenodo.org` first with `--sandbox` before touching
production.

## Citation

If you use this pipeline or the dataset it produces, please cite:

```
Dao, T.M., Shehu, A., Kabir, A. ViralAI-KB: AI-Ready Biophysical
Knowledge Base for Viral Proteins. [venue, year, DOI to be added once
published]
```
