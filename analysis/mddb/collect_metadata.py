#!/usr/bin/env python3
"""
collect_metadata.py

For the comprehensive-review workflow (data/analysis/), NOT the fetch/dedup
pipeline in data/preprocess/. Runs only STEP1-3 of the MDposit pipeline
(fetch every project, keep records with a UniProt reference, verify viral
organism) -- no sequence flattening, no dedup, no CSV. Every confirmed-
viral project is kept, one file per project, with its full live metadata
written to its own file in simulations/ for manual review.

  STEP 1. Fetch every project record from MDposit.
  STEP 2. Keep only ones that have at least one entry in
          metadata.REFERENCES (a UniProt accession -- MDposit never uses
          PDB ids for this field, so there is no PDB-only fallback here,
          unlike the MDRepo pipeline).
  STEP 3. Keep only ones confirmed VIRAL: verify every distinct UniProt
          accession found in step 2 against the live UniProt API,
          falling back to an NCBI taxonomy lookup by tax id when UniProt
          doesn't return a lineage.
  (new)   For every STEP3 survivor, fetch its full project document
          (GET /projects/{accession}) and write it to
          simulations/{accession}_{uniprot_id}.json

All network calls retry on transient failure AND on non-200 responses
(e.g. 429 Too Many Requests) with backoff, and detail requests are
throttled -- this is the exact bug found and fixed in
data/analysis/mdrepo/collect_metadata.py after MDRepo's server started
rate-limiting an unthrottled run.

Usage:
    pip install requests --break-system-packages
    python3 collect_metadata.py
"""

import json
import os
import re
import time
from typing import Optional

import requests

# ---------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------
API_BASE = "https://mdposit.mddbr.eu/api/rest/current"
PAGE_LIMIT = 100  # API-enforced maximum
REQUEST_HEADERS = {"User-Agent": "viralplm-analysis-mddb/1.0"}
REQUEST_TIMEOUT = 30

RAW_CACHE_FILE = "mdposit_raw_cache.json"
VIRAL_LIST_FILE = "mdposit_viral_list.json"  # STEP3 output, before per-project detail fetch
SIMULATIONS_DIR = "simulations"
FAILURES_FILE = "mdposit_detail_failures.json"  # accession -> reason, for diagnosis/retry

UNIPROT_BATCH_SIZE = 50
CHECKPOINT_EVERY_PAGES = 10
DETAIL_CHECKPOINT_EVERY = 50
MAX_RETRIES = 5
RETRY_SLEEP = 3
REQUEST_DELAY = 0.3  # throttle between detail requests -- see module docstring

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# MDposit's metadata.REFERENCES occasionally holds placeholders (e.g.
# "notfound", "noref") instead of a real UniProt accession.
UNIPROT_ACCESSION_RE = re.compile(
    r'^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-\d+)?$'
)

# Listing projection -- only what's needed for STEP2/3; the FULL document
# is fetched separately per confirmed-viral accession.
PROJECT_PROJECTION = {
    "accession": 1,
    "metadata.REFERENCES": 1,
}


# ---------------------------------------------------------------------
# STEP 1: Fetch every project record from MDposit (retry-safe)
# ---------------------------------------------------------------------
def fetch_all_projects() -> list:
    if os.path.exists(RAW_CACHE_FILE):
        print(f"  Loading cached raw records from {RAW_CACHE_FILE}...")
        with open(RAW_CACHE_FILE) as f:
            return json.load(f)

    projections = [json.dumps({k: v}) for k, v in PROJECT_PROJECTION.items()]
    all_records = []
    page = 1
    total = None
    while True:
        params = {"projection": projections, "limit": PAGE_LIMIT, "page": page}
        data = None
        for attempt in range(MAX_RETRIES):
            try:
                resp = requests.get(f"{API_BASE}/projects", params=params,
                                     headers=REQUEST_HEADERS, timeout=REQUEST_TIMEOUT)
                if resp.status_code == 200:
                    data = resp.json()
                    break
                if attempt < MAX_RETRIES - 1:
                    print(f"    [warn] page {page} got http_{resp.status_code}, retrying...")
                    time.sleep(RETRY_SLEEP * (attempt + 1))
            except Exception as e:
                if attempt == MAX_RETRIES - 1:
                    raise
                print(f"    [warn] page {page} failed ({e}), retrying...")
                time.sleep(RETRY_SLEEP)
        if data is None:
            raise RuntimeError(f"page {page} failed after {MAX_RETRIES} attempts")

        batch = data.get("projects", [])
        all_records.extend(batch)
        if total is None:
            total = data.get("filteredCount", 0)
        print(f"    page {page}: {len(all_records)}/{total} records fetched")

        if not batch or len(all_records) >= total:
            break
        page += 1

        if page % CHECKPOINT_EVERY_PAGES == 0:
            with open(RAW_CACHE_FILE, "w") as f:
                json.dump(all_records, f)
            print(f"      [checkpoint saved at page {page}]")

    if total is not None and len(all_records) < total:
        raise RuntimeError(f"Only fetched {len(all_records)}/{total} records -- refusing to "
                            f"cache a truncated listing. Re-run to retry.")

    with open(RAW_CACHE_FILE, "w") as f:
        json.dump(all_records, f)
    return all_records


# ---------------------------------------------------------------------
# STEP 2: Keep only records with at least one UniProt reference
# ---------------------------------------------------------------------
def filter_has_identifier(records: list) -> list:
    kept = []
    for rec in records:
        refs = (rec.get("metadata") or {}).get("REFERENCES") or []
        if refs:
            kept.append(rec)
    return kept


# ---------------------------------------------------------------------
# STEP 3: Keep only records confirmed viral (retry-safe)
# ---------------------------------------------------------------------
def uniprot_batch_lookup(accessions: list) -> dict:
    if not accessions:
        return {}
    query = " OR ".join(f"accession:{a}" for a in accessions)
    data = None
    for attempt in range(MAX_RETRIES):
        try:
            resp = requests.get(
                f"{UNIPROT_BASE}/search",
                params={
                    "query": query,
                    "fields": "accession,organism_name,lineage",
                    "format": "json",
                    "size": len(accessions),
                },
                timeout=30,
            )
            resp.raise_for_status()
            data = resp.json()
            break
        except Exception as e:
            if attempt == MAX_RETRIES - 1:
                print(f"      [warn] batch UniProt lookup failed for a chunk after "
                      f"{MAX_RETRIES} attempts: {e}")
                return {}
            time.sleep(RETRY_SLEEP)
    if data is None:
        return {}

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


def uniprot_lookup_all(accessions: list) -> dict:
    results = {}
    accessions = list(accessions)
    for i in range(0, len(accessions), UNIPROT_BATCH_SIZE):
        chunk = accessions[i:i + UNIPROT_BATCH_SIZE]
        print(f"    looking up UniProt batch {i // UNIPROT_BATCH_SIZE + 1} "
              f"({len(chunk)} accessions)...")
        results.update(uniprot_batch_lookup(chunk))
        time.sleep(0.3)
    return results


def taxid_is_virus(tax_id: Optional[int]) -> bool:
    if not tax_id:
        return False
    for attempt in range(MAX_RETRIES):
        try:
            resp = requests.get(
                f"{NCBI_EUTILS}/efetch.fcgi",
                params={"db": "taxonomy", "id": tax_id, "rettype": "xml"},
                timeout=15,
            )
            return "Viruses" in resp.text
        except Exception:
            if attempt == MAX_RETRIES - 1:
                return False
            time.sleep(RETRY_SLEEP)
    return False


def filter_viral(records: list) -> list:
    uniprot_accessions = set()
    skipped = set()
    for rec in records:
        refs = (rec.get("metadata") or {}).get("REFERENCES") or []
        for ref in refs:
            if UNIPROT_ACCESSION_RE.match(ref):
                uniprot_accessions.add(ref)
            else:
                skipped.add(ref)
    if skipped:
        print(f"    skipping {len(skipped)} malformed reference(s) (not real UniProt "
              f"accessions): {sorted(skipped)[:10]}{'...' if len(skipped) > 10 else ''}")

    print(f"    verifying {len(uniprot_accessions)} distinct UniProt accession(s)...")
    uniprot_info = uniprot_lookup_all(list(uniprot_accessions))
    viral_accessions = {
        acc for acc, info in uniprot_info.items()
        if (("Viruses" in info["lineage"]) if info["lineage"] else taxid_is_virus(info["tax_id"]))
    }

    kept = []
    for rec in records:
        refs = (rec.get("metadata") or {}).get("REFERENCES") or []
        matched = [r for r in refs if r in viral_accessions]
        if not matched:
            continue
        rec = dict(rec)
        rec["_viral_references"] = matched
        rec["_uniprot_id"] = matched[0]
        rec["_organism"] = uniprot_info.get(matched[0], {}).get("organism_name", "")
        kept.append(rec)

    return kept


# ---------------------------------------------------------------------
# Per-project metadata export
# ---------------------------------------------------------------------
def fetch_project_detail(accession: str) -> tuple:
    """Returns (detail_dict_or_None, reason_str_or_None). Same retry/
    backoff-on-non-200 pattern as data/analysis/mdrepo/collect_metadata.py."""
    last_reason = None
    for attempt in range(MAX_RETRIES):
        try:
            resp = requests.get(f"{API_BASE}/projects/{accession}",
                                 headers=REQUEST_HEADERS, timeout=REQUEST_TIMEOUT)
            if resp.status_code == 200:
                return resp.json(), None
            if resp.status_code == 404:
                return None, "404 not found"
            last_reason = f"http_{resp.status_code}"
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_SLEEP * (attempt + 1))
        except Exception as e:
            last_reason = str(e)
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_SLEEP)
    return None, f"failed after {MAX_RETRIES} attempts: {last_reason}"


def id_slug(rec: dict) -> str:
    """accession + uniprot_id (mddb's REFERENCES are always UniProt, never
    PDB, so there is no pdb_id fallback branch here)."""
    return f"{rec['accession']}_{rec['_uniprot_id']}"


# ---------------------------------------------------------------------
# STEP 4: Exclude heteromeric complexes (chains with genuinely different
# sequences). Single-chain records and homo-oligomers (multiple chains,
# all identical sequence) are kept. detail.metadata.PROTSEQ is already a
# list with one entry per chain, so this only requires the already-fetched
# detail record -- no extra network calls.
# ---------------------------------------------------------------------
def is_single_chain_or_homo_oligomer(detail: dict) -> bool:
    md = (detail or {}).get("metadata") or {}
    seqs = md.get("PROTSEQ")
    if not seqs:
        return True  # no sequence to check; harmonization handles this separately
    if isinstance(seqs, str):
        seqs = [seqs]
    return len(set(seqs)) <= 1


def export_metadata(viral_records: list):
    os.makedirs(SIMULATIONS_DIR, exist_ok=True)
    n_written = 0
    n_skipped = 0
    n_excluded_complex = 0
    failures = {}

    for i, rec in enumerate(viral_records, start=1):
        slug = id_slug(rec)
        out_path = os.path.join(SIMULATIONS_DIR, f"{slug}.json")
        if os.path.exists(out_path):
            n_skipped += 1
            continue

        detail, reason = fetch_project_detail(rec["accession"])
        time.sleep(REQUEST_DELAY)
        if detail is None:
            failures[rec["accession"]] = reason
            continue

        if not is_single_chain_or_homo_oligomer(detail):
            n_excluded_complex += 1
            continue

        merged = {
            "accession": rec["accession"],
            "uniprot_id": rec["_uniprot_id"],
            "viral_uniprot_ids": rec["_viral_references"],
            "organism": rec["_organism"],
            "detail": detail,
        }
        with open(out_path, "w") as f:
            json.dump(merged, f, indent=2)
        n_written += 1

        if i % DETAIL_CHECKPOINT_EVERY == 0:
            print(f"    ...{i}/{len(viral_records)} processed "
                  f"({n_written} written, {n_skipped} already existed, "
                  f"{n_excluded_complex} excluded as heteromeric complexes, {len(failures)} failed)")
            with open(FAILURES_FILE, "w") as f:
                json.dump(failures, f, indent=2)

    with open(FAILURES_FILE, "w") as f:
        json.dump(failures, f, indent=2)

    print(f"\n  Done: {n_written} written, {n_skipped} already existed, "
          f"{n_excluded_complex} excluded as heteromeric complexes, "
          f"{len(failures)} failed, out of {len(viral_records)} confirmed-viral projects.")
    if failures:
        print(f"  Failure reasons logged to {FAILURES_FILE} -- rerun the script to retry "
              f"them (already-written files are skipped, so only failures get re-fetched).")


def main():
    print("STEP 1: Fetching every project record from MDposit...")
    all_records = fetch_all_projects()
    print(f"  {len(all_records)} total record(s).\n")

    print("STEP 2: Filtering to records with at least one UniProt reference...")
    has_identifier = filter_has_identifier(all_records)
    print(f"  {len(has_identifier)} of {len(all_records)} have at least one identifier.\n")

    print("STEP 3: Filtering to confirmed-viral records...")
    if os.path.exists(VIRAL_LIST_FILE):
        print(f"  Loading cached STEP3 result from {VIRAL_LIST_FILE}...")
        with open(VIRAL_LIST_FILE) as f:
            viral = json.load(f)
    else:
        viral = filter_viral(has_identifier)
        with open(VIRAL_LIST_FILE, "w") as f:
            json.dump(viral, f)
    print(f"  {len(viral)} of {len(has_identifier)} confirmed viral.\n")

    print(f"Exporting per-project metadata to {SIMULATIONS_DIR}/ ...")
    export_metadata(viral)


if __name__ == "__main__":
    main()
