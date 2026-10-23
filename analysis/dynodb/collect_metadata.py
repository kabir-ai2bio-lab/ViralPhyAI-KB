#!/usr/bin/env python3
"""
collect_metadata.py

Comprehensive-review workflow for the DynoDB dataset
(Hu, Yang & Ren, "A Dataset of Small Protein Conformational Ensembles
from All-Atom Molecular Dynamics Simulations", Scientific Data 2026;
https://huggingface.co/datasets/renzilin/DynoDB).

NOTE on trajectory data: the `trajectory_frames/` folder IS real and
public (verified directly -- files follow the `{PDB_ID}_frame.pdb`
convention, a 1,000-frame multi-model PDB, ~30MB each, matching the
paper's own stated methodology; an earlier check that concluded
otherwise was testing the wrong filenames, e.g. `_frames.pdb` plural
and `.xtc`, which 404 regardless of whether the real files exist). This
script still only collects metadata + the static starting structure per
confirmed-viral entry -- the trajectory files themselves are pulled
separately, only for the records that win a final-selection slot, by
data/analysis/dynodb/download_trajectories.py.

  STEP 1. Fetch the full metadata table (8,385 rows, 72 columns) via
          HuggingFace's public dataset-server rows API.
  STEP 2. Verify viral organism per row, in priority order per row:
            a. UniProt_ID present -> batch UniProt lookup (same pattern
               as data/analysis/mdrepo and mddb).
            b. else Organism name present -> NCBI taxonomy lookup by
               name (same pattern as mddb/mdrepo's PDB-only fallback).
            c. else (both null) -> PDB_ID -> RCSB polymer_entity
               organism lookup, same as data/analysis/mdcath, but
               simpler here: DynoDB already provides the as-simulated
               sequence directly (Sequence_y(fixed)), so there is no
               need to resolve a specific chain/entity for sequence
               slicing like mdCATH required -- RCSB is consulted here
               ONLY to get an organism when both other fields are null.
  (new)   For every confirmed-viral row, write one file to
          simulations/{PDB_ID}_{uniprot_id_or_pdb_id}.json containing
          the full raw metadata row plus the raw static PDB structure
          text (best-effort; the closest thing to "full detail" this
          dataset can offer until trajectories are released).

No deduplication is performed, within this dataset or against any
other -- consistent with mdrepo/mddb/mdcath: maximize coverage, one
file per confirmed-viral entry.

Usage:
    pip install requests --break-system-packages
    python3 collect_metadata.py
"""

import json
import os
import time
from typing import Optional

import requests

# ---------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------
HF_DATASET = "renzilin/DynoDB"
ROWS_API = "https://datasets-server.huggingface.co/rows"
STATIC_STRUCTURE_URL = ("https://huggingface.co/datasets/renzilin/DynoDB/resolve/main/"
                         "static_structures/{pdb_id}_static.pdb")

METADATA_CACHE_FILE = "dynodb_metadata.json"
SIMULATIONS_DIR = "simulations"
FAILURES_FILE = "dynodb_resolve_failures.json"
ROW_PAGE_SIZE = 100

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
RCSB_ENTITY_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"

REQUEST_HEADERS = {"User-Agent": "viralplm-analysis-dynodb/1.0"}
UNIPROT_BATCH_SIZE = 50
MAX_RETRIES = 5
RETRY_SLEEP = 3
REQUEST_DELAY = 0.3


def _get_with_retry(url: str, params: dict = None, timeout: int = 20) -> Optional[requests.Response]:
    for attempt in range(MAX_RETRIES):
        try:
            resp = requests.get(url, params=params, headers=REQUEST_HEADERS, timeout=timeout)
            if resp.status_code == 200:
                return resp
            if resp.status_code == 404:
                return resp
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_SLEEP * (attempt + 1))
        except Exception:
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_SLEEP)
    return None


# ---------------------------------------------------------------------
# STEP 1: Fetch the full metadata table
# ---------------------------------------------------------------------
def fetch_full_metadata() -> list:
    if os.path.exists(METADATA_CACHE_FILE):
        print(f"  Loading cached metadata from {METADATA_CACHE_FILE}...")
        with open(METADATA_CACHE_FILE) as f:
            return json.load(f)

    rows = []
    offset = 0
    total = None
    while total is None or offset < total:
        resp = _get_with_retry(ROWS_API, params={
            "dataset": HF_DATASET, "config": "default", "split": "train",
            "offset": offset, "length": ROW_PAGE_SIZE,
        }, timeout=30)
        if resp is None or resp.status_code != 200:
            raise RuntimeError(f"failed to fetch rows at offset {offset} after retries")
        data = resp.json()
        total = data.get("num_rows_total", total)
        batch = [r["row"] for r in data.get("rows", [])]
        rows.extend(batch)
        print(f"    {len(rows)}/{total} rows fetched")
        offset += ROW_PAGE_SIZE
        if not batch:
            break

    with open(METADATA_CACHE_FILE, "w") as f:
        json.dump(rows, f)
    return rows


# ---------------------------------------------------------------------
# STEP 2a: UniProt batch verification (same pattern as mddb/mdrepo)
# ---------------------------------------------------------------------
def uniprot_batch_lookup(accessions: list) -> dict:
    if not accessions:
        return {}
    query = " OR ".join(f"accession:{a}" for a in accessions)
    resp = _get_with_retry(f"{UNIPROT_BASE}/search", params={
        "query": query, "fields": "accession,organism_name,lineage",
        "format": "json", "size": len(accessions),
    }, timeout=30)
    if resp is None or resp.status_code != 200:
        print(f"      [warn] batch UniProt lookup failed for a chunk of {len(accessions)}")
        return {}
    data = resp.json()
    out = {}
    for entry in data.get("results", []):
        acc = entry.get("primaryAccession")
        organism = entry.get("organism", {})
        lineage = organism.get("lineage", [])
        if lineage and isinstance(lineage[0], dict):
            lineage = [l.get("scientificName", "") for l in lineage]
        out[acc] = {
            "organism_name": organism.get("scientificName", ""),
            "tax_id": organism.get("taxonId"),
            "lineage": lineage,
        }
    return out


def taxid_is_virus(tax_id) -> bool:
    if not tax_id:
        return False
    resp = _get_with_retry(f"{NCBI_EUTILS}/efetch.fcgi",
                            params={"db": "taxonomy", "id": tax_id, "rettype": "xml"}, timeout=15)
    return resp is not None and resp.status_code == 200 and "Viruses" in resp.text


def verify_uniprot_ids(accessions: list) -> dict:
    """Returns {accession: is_viral bool}."""
    info = {}
    accs = list(accessions)
    for i in range(0, len(accs), UNIPROT_BATCH_SIZE):
        chunk = accs[i:i + UNIPROT_BATCH_SIZE]
        print(f"    UniProt batch {i // UNIPROT_BATCH_SIZE + 1} ({len(chunk)} accessions)...")
        info.update(uniprot_batch_lookup(chunk))
        time.sleep(0.3)
    return {
        acc: (("Viruses" in d["lineage"]) if d["lineage"] else taxid_is_virus(d["tax_id"]))
        for acc, d in info.items()
    }, info


# ---------------------------------------------------------------------
# STEP 2b: Organism-name fallback (UniProt_ID null, Organism present)
# ---------------------------------------------------------------------
def name_is_virus(organism_name: str) -> bool:
    r = _get_with_retry(f"{NCBI_EUTILS}/esearch.fcgi",
                         params={"db": "taxonomy", "term": organism_name, "retmode": "json"}, timeout=15)
    if r is None or r.status_code != 200:
        return False
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        return False
    r2 = _get_with_retry(f"{NCBI_EUTILS}/efetch.fcgi",
                          params={"db": "taxonomy", "id": ids[0], "rettype": "xml"}, timeout=15)
    return r2 is not None and r2.status_code == 200 and "Viruses" in r2.text


# ---------------------------------------------------------------------
# STEP 2c: PDB fallback (both UniProt_ID and Organism null)
# ---------------------------------------------------------------------
def pdb_id_is_virus(pdb_id: str) -> tuple:
    """Returns (is_viral, organism_name_or_None). Checks every polymer
    entity of the entry; any one being viral is enough."""
    entry_resp = _get_with_retry(RCSB_ENTRY_URL.format(pdb_id=pdb_id))
    if entry_resp is None or entry_resp.status_code != 200:
        return False, None
    entity_ids = entry_resp.json().get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids") or []
    for eid in entity_ids:
        time.sleep(REQUEST_DELAY)
        resp = _get_with_retry(RCSB_ENTITY_URL.format(pdb_id=pdb_id, entity_id=eid))
        if resp is None or resp.status_code != 200:
            continue
        for org in (resp.json().get("rcsb_entity_source_organism") or []):
            lineage_names = {seg.get("name") for seg in (org.get("taxonomy_lineage") or [])}
            if "Viruses" in lineage_names:
                return True, org.get("scientific_name") or org.get("ncbi_scientific_name")
    return False, None


# ---------------------------------------------------------------------
# Static structure fetch (best-effort; trajectories aren't out yet)
# ---------------------------------------------------------------------
def fetch_static_structure(pdb_id: str) -> Optional[str]:
    resp = _get_with_retry(STATIC_STRUCTURE_URL.format(pdb_id=pdb_id), timeout=30)
    if resp is None or resp.status_code != 200:
        return None
    return resp.text


def id_slug(pdb_id: str, uniprot_id: Optional[str]) -> str:
    identifier = uniprot_id or pdb_id
    return f"{pdb_id}_{identifier}"


def main():
    print("STEP 1: Fetching full DynoDB metadata table...")
    rows = fetch_full_metadata()
    print(f"  {len(rows)} row(s).\n")

    print("STEP 2a: Verifying viral organism via UniProt_ID...")
    distinct_uniprot = sorted({r["UniProt_ID"] for r in rows if r.get("UniProt_ID")})
    print(f"  {len(distinct_uniprot)} distinct UniProt_ID(s) to verify...")
    uniprot_viral_map, uniprot_info = verify_uniprot_ids(distinct_uniprot)
    print(f"  {sum(uniprot_viral_map.values())} confirmed viral.\n")

    print("STEP 2b: Verifying viral organism via Organism name (rows with no UniProt_ID)...")
    need_name_check = sorted({r["Organism"] for r in rows
                               if not r.get("UniProt_ID") and r.get("Organism")})
    print(f"  {len(need_name_check)} distinct organism name(s) to verify...")
    name_viral_map = {}
    for i, name in enumerate(need_name_check, start=1):
        name_viral_map[name] = name_is_virus(name)
        if i % 50 == 0:
            print(f"    ...{i}/{len(need_name_check)} organism names checked")
        time.sleep(0.2)
    print(f"  {sum(name_viral_map.values())} confirmed viral.\n")

    print("STEP 2c: Verifying viral organism via PDB_ID (rows with neither field)...")
    need_pdb_check = sorted({r["PDB_ID"] for r in rows
                              if not r.get("UniProt_ID") and not r.get("Organism") and r.get("PDB_ID")})
    print(f"  {len(need_pdb_check)} distinct PDB_ID(s) to verify...")
    pdb_viral_map = {}
    for i, pdb_id in enumerate(need_pdb_check, start=1):
        is_viral, organism = pdb_id_is_virus(pdb_id)
        pdb_viral_map[pdb_id] = (is_viral, organism)
        if i % 50 == 0:
            print(f"    ...{i}/{len(need_pdb_check)} PDB IDs checked")
    print(f"  {sum(v[0] for v in pdb_viral_map.values())} confirmed viral.\n")

    print("Exporting metadata + static structure for every confirmed-viral row...")
    os.makedirs(SIMULATIONS_DIR, exist_ok=True)
    n_written = 0
    n_skipped = 0
    n_not_viral = 0
    n_structure_missing = 0
    n_excluded_complex = 0
    failures = {}

    for i, row in enumerate(rows, start=1):
        uniprot_id = row.get("UniProt_ID")
        organism_field = row.get("Organism")
        pdb_id = row.get("PDB_ID")

        if uniprot_id:
            method = "uniprot"
            is_viral = uniprot_viral_map.get(uniprot_id, False)
            organism = uniprot_info.get(uniprot_id, {}).get("organism_name") or organism_field
        elif organism_field:
            method = "organism_name"
            is_viral = name_viral_map.get(organism_field, False)
            organism = organism_field
        elif pdb_id:
            method = "pdb"
            is_viral, organism = pdb_viral_map.get(pdb_id, (False, None))
        else:
            failures[row.get("PDB_ID", f"row_{i}")] = "no PDB_ID, UniProt_ID, or Organism to verify with"
            continue

        if not is_viral:
            n_not_viral += 1
            continue

        # Exclude heteromeric complexes: DynoDB's own metadata already
        # carries a precomputed Unique_Protein_Chains count (distinct
        # chain sequences, not raw chain count), so no parsing needed.
        # Single-chain rows and homo-oligomers (Unique_Protein_Chains <= 1)
        # are kept; in practice this excludes nothing today (DynoDB's small
        # peptide/protein scope means every confirmed-viral row so far is
        # already single-chain), but is checked explicitly in case that
        # changes as the source grows.
        n_unique_chains = row.get("Unique_Protein_Chains")
        if n_unique_chains is not None and n_unique_chains > 1:
            n_excluded_complex += 1
            continue

        slug = id_slug(pdb_id, uniprot_id)
        out_path = os.path.join(SIMULATIONS_DIR, f"{slug}.json")
        if os.path.exists(out_path):
            n_skipped += 1
            continue

        static_pdb_text = fetch_static_structure(pdb_id)
        time.sleep(REQUEST_DELAY)
        if static_pdb_text is None:
            n_structure_missing += 1

        merged = {
            "pdb_id": pdb_id,
            "uniprot_id": uniprot_id,
            "organism": organism,
            "verification_method": method,
            "sequence": row.get("Sequence_y(fixed)"),
            "trajectory_data_status": "available at trajectory_frames/{pdb_id}_frame.pdb -- not "
                                       "downloaded here, see download_trajectories.py",
            "metadata_row": row,
            "static_structure_pdb": static_pdb_text,
        }
        with open(out_path, "w") as f:
            json.dump(merged, f, indent=2)
        n_written += 1

        if i % 500 == 0:
            print(f"    ...{i}/{len(rows)} rows processed "
                  f"({n_written} written, {n_not_viral} not viral, {n_skipped} already existed)")
            with open(FAILURES_FILE, "w") as f:
                json.dump(failures, f, indent=2)

    with open(FAILURES_FILE, "w") as f:
        json.dump(failures, f, indent=2)

    print(f"\n  Done: {n_written} written, {n_not_viral} confirmed non-viral, "
          f"{n_excluded_complex} excluded as heteromeric complexes, "
          f"{n_skipped} already existed, {len(failures)} failed to resolve, "
          f"{n_structure_missing} written without a static structure (fetch failed), "
          f"out of {len(rows)} total rows.")


if __name__ == "__main__":
    main()
