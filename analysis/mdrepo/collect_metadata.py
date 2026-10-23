#!/usr/bin/env python3
"""
collect_metadata.py

For the comprehensive-review workflow (data/analysis/), NOT the fetch/dedup
pipeline in data/preprocess/. Runs only STEP1-3 of the MDRepo pipeline
(fetch every simulation, keep records with a uniprot_id/pdb_id, verify
viral organism) -- no sequence flattening, no dedup, no CSV. Every
confirmed-viral simulation is kept, one-per-simulation (not one-per-
sequence), and its full live metadata is written to its own file in
simulations/ for manual review.

  STEP 1. Fetch every simulation record from MDRepo.
  STEP 2. Keep only ones that have a uniprot_id OR a pdb_id.
  STEP 3. Keep only ones confirmed VIRAL (verify via uniprot_id; fall
          back to pdb_id if there's no UniProt mapping).
  (new)   For every STEP3 survivor, fetch its full simulation detail
          (GET /simulation_view/{id}) and write it to
          simulations/{sim_id}_{uniprot_id_or_pdb_id}.json

All network calls retry on transient failure (timeouts, connection
resets) instead of silently dropping records, which is the bug found
and fixed in data/preprocess/mdrepo/fetch_v2.py.

Usage:
    pip install requests --break-system-packages
    python3 collect_metadata.py
"""

import csv
import json
import os
import time
from typing import Optional

import requests

# ---------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------
MDREPO_API_URL = "https://mdrepo.org/api/v1/simulation_list"
MDREPO_DETAIL_URL = "https://mdrepo.org/api/v1/simulation_view/{sim_id}"
PAGE_LIMIT = 100
RAW_CACHE_FILE = "mdrepo_raw_cache.json"
VIRAL_LIST_FILE = "mdrepo_viral_list.json"  # STEP3 output, before per-simulation detail fetch
SIMULATIONS_DIR = "simulations"
FAILURES_FILE = "mdrepo_detail_failures.json"  # sim_id -> reason, for diagnosis/retry
REQUEST_HEADERS = {"User-Agent": "viralplm-analysis-mdrepo/1.0"}
ONLY_PUBLIC = True
UNIPROT_BATCH_SIZE = 50
CHECKPOINT_EVERY_PAGES = 20
DETAIL_CHECKPOINT_EVERY = 100
MAX_RETRIES = 5
RETRY_SLEEP = 3
REQUEST_DELAY = 0.3  # throttle between detail requests -- MDRepo's nginx returns 429
                      # after sustained unthrottled load (observed at ~250 requests
                      # within ~75s in testing)

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"
RCSB_ENTITY = "https://data.rcsb.org/rest/v1/core/polymer_entity"
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


# ---------------------------------------------------------------------
# STEP 1: Fetch every simulation record from MDRepo (retry-safe)
# ---------------------------------------------------------------------
def fetch_all_simulations() -> list:
    if os.path.exists(RAW_CACHE_FILE):
        print(f"  Loading cached raw records from {RAW_CACHE_FILE}...")
        with open(RAW_CACHE_FILE) as f:
            return json.load(f)

    all_records = []
    url = MDREPO_API_URL
    params = {"limit": PAGE_LIMIT, "ordering": "-id"}
    page = 0
    total = None

    while url:
        page += 1
        for attempt in range(MAX_RETRIES):
            try:
                resp = requests.get(url, params=params, headers=REQUEST_HEADERS, timeout=30)
                resp.raise_for_status()
                data = resp.json()
                break
            except Exception as e:
                if attempt == MAX_RETRIES - 1:
                    raise
                print(f"    [warn] page {page} failed ({e}), retrying...")
                time.sleep(RETRY_SLEEP)

        all_records.extend(data.get("results", []))
        total = data.get("count")
        print(f"    page {page}: {len(all_records)}/{total} records fetched")

        next_url = data.get("next")
        url = next_url.replace("http://", "https://", 1) if next_url else None
        params = None

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
# STEP 2: Keep only records with a uniprot_id or a pdb_id
# ---------------------------------------------------------------------
def filter_has_identifier(records: list) -> list:
    kept = []
    for rec in records:
        has_uniprot = bool(rec.get("uniprot"))
        has_pdb = bool((rec.get("pdb") or {}).get("pdb_id"))
        if has_uniprot or has_pdb:
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


def pdb_organism_name(pdb_id: str) -> Optional[str]:
    for attempt in range(MAX_RETRIES):
        try:
            resp = requests.get(f"{RCSB_ENTITY}/{pdb_id}/1", timeout=15)
            if resp.status_code != 200:
                return None
            orgs = resp.json().get("rcsb_entity_source_organism", [])
            if orgs:
                return orgs[0].get("scientific_name") or orgs[0].get("ncbi_scientific_name")
            return None
        except Exception:
            if attempt == MAX_RETRIES - 1:
                return None
            time.sleep(RETRY_SLEEP)
    return None


def name_is_virus(organism_name: str) -> bool:
    if not organism_name:
        return False
    for attempt in range(MAX_RETRIES):
        try:
            r = requests.get(
                f"{NCBI_EUTILS}/esearch.fcgi",
                params={"db": "taxonomy", "term": organism_name, "retmode": "json"},
                timeout=15,
            )
            ids = r.json().get("esearchresult", {}).get("idlist", [])
            if not ids:
                return False
            r2 = requests.get(
                f"{NCBI_EUTILS}/efetch.fcgi",
                params={"db": "taxonomy", "id": ids[0], "rettype": "xml"},
                timeout=15,
            )
            return "Viruses" in r2.text
        except Exception:
            if attempt == MAX_RETRIES - 1:
                return False
            time.sleep(RETRY_SLEEP)
    return False


def filter_viral(records: list) -> list:
    uniprot_accessions = set()
    pdb_only_ids = set()
    for rec in records:
        uniprot_entries = rec.get("uniprot") or []
        if uniprot_entries:
            for u in uniprot_entries:
                acc = u.get("uniprot_id")
                if acc:
                    uniprot_accessions.add(acc)
        else:
            pdb_id = (rec.get("pdb") or {}).get("pdb_id")
            if pdb_id:
                pdb_only_ids.add(pdb_id)

    print(f"    verifying {len(uniprot_accessions)} distinct UniProt accession(s)...")
    uniprot_info = uniprot_lookup_all(list(uniprot_accessions))
    viral_accessions = {
        acc for acc, info in uniprot_info.items()
        if (("Viruses" in info["lineage"]) if info["lineage"] else taxid_is_virus(info["tax_id"]))
    }

    print(f"    verifying {len(pdb_only_ids)} distinct PDB-only identifier(s)...")
    pdb_organism = {}
    for pdb_id in pdb_only_ids:
        org_name = pdb_organism_name(pdb_id)
        if org_name and name_is_virus(org_name):
            pdb_organism[pdb_id] = org_name

    kept = []
    for rec in records:
        uniprot_entries = rec.get("uniprot") or []
        if uniprot_entries:
            matched = [u.get("uniprot_id") for u in uniprot_entries
                       if u.get("uniprot_id") in viral_accessions]
            if not matched:
                continue
            rec = dict(rec)
            rec["_uniprot_id"] = matched[0]
            rec["_pdb_id"] = (rec.get("pdb") or {}).get("pdb_id", "")
            rec["_organism"] = uniprot_info[matched[0]]["organism_name"]
            rec["_verification_method"] = "uniprot"
            kept.append(rec)
        else:
            pdb_id = (rec.get("pdb") or {}).get("pdb_id")
            if pdb_id not in pdb_organism:
                continue
            rec = dict(rec)
            rec["_uniprot_id"] = ""
            rec["_pdb_id"] = pdb_id
            rec["_organism"] = pdb_organism[pdb_id]
            rec["_verification_method"] = "pdb"
            kept.append(rec)

    return kept


# ---------------------------------------------------------------------
# Per-simulation metadata export
# ---------------------------------------------------------------------
def fetch_simulation_detail(sim_id) -> tuple:
    """Returns (detail_dict_or_None, reason_str_or_None). A 404 is a
    definitive miss (don't retry). Anything else non-200 -- most notably
    429 Too Many Requests, which MDRepo's nginx starts returning under
    sustained request volume -- is retried with backoff, same as a
    network exception."""
    last_reason = None
    for attempt in range(MAX_RETRIES):
        try:
            resp = requests.get(
                MDREPO_DETAIL_URL.format(sim_id=sim_id), headers=REQUEST_HEADERS, timeout=20
            )
            if resp.status_code == 200:
                return resp.json(), None
            if resp.status_code == 404:
                return None, "404 not found"
            last_reason = f"http_{resp.status_code}"
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_SLEEP * (attempt + 1))  # back off harder each retry
        except Exception as e:
            last_reason = str(e)
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_SLEEP)
    return None, f"failed after {MAX_RETRIES} attempts: {last_reason}"


def id_slug(rec: dict) -> str:
    """simulation_id + uniprot_id (or pdb_id if no uniprot match)."""
    identifier = rec["_uniprot_id"] or rec["_pdb_id"] or "unknown"
    return f"{rec['id']}_{identifier}"


# ---------------------------------------------------------------------
# STEP 4: Exclude heteromeric complexes (chains with genuinely different
# sequences). Single-chain records and homo-oligomers (multiple chains,
# all identical sequence) are kept. detail.fasta_sequence carries one
# '>N' FASTA block per chain, so this only requires the already-fetched
# detail record -- no extra network calls.
# ---------------------------------------------------------------------
def is_single_chain_or_homo_oligomer(detail: dict) -> bool:
    fasta = (detail or {}).get("fasta_sequence") or ""
    if not fasta:
        return True  # no sequence to check; harmonization drops these separately
    chains, current = [], []
    for line in fasta.strip().split("\n"):
        if line.startswith(">"):
            if current:
                chains.append("".join(current))
                current = []
        else:
            current.append(line.strip())
    if current:
        chains.append("".join(current))
    if not chains:
        return True
    return len(set(chains)) <= 1


def export_metadata(viral_records: list):
    os.makedirs(SIMULATIONS_DIR, exist_ok=True)
    n_written = 0
    n_skipped = 0
    n_excluded_complex = 0
    failures = {}  # sim_id -> reason, for diagnosis and easy manual retry

    for i, rec in enumerate(viral_records, start=1):
        slug = id_slug(rec)
        out_path = os.path.join(SIMULATIONS_DIR, f"{slug}.json")
        if os.path.exists(out_path):
            n_skipped += 1
            continue

        detail, reason = fetch_simulation_detail(rec["id"])
        time.sleep(REQUEST_DELAY)  # throttle so we don't trip the rate limit again
        if detail is None:
            failures[str(rec["id"])] = reason
            continue

        if not is_single_chain_or_homo_oligomer(detail):
            n_excluded_complex += 1
            continue

        merged = {
            "simulation_id": rec["id"],
            "slug": rec.get("slug"),
            "uniprot_id": rec["_uniprot_id"],
            "pdb_id": rec["_pdb_id"],
            "organism": rec["_organism"],
            "verification_method": rec["_verification_method"],
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
          f"{len(failures)} failed, out of {len(viral_records)} confirmed-viral simulations.")
    if failures:
        print(f"  Failure reasons logged to {FAILURES_FILE} -- rerun the script to retry "
              f"them (already-written files are skipped, so only failures get re-fetched).")


def main():
    print("STEP 1: Fetching every simulation record from MDRepo...")
    all_records = fetch_all_simulations()
    all_records = [r for r in all_records
                   if not (ONLY_PUBLIC and (not r.get("is_public") or r.get("is_embargoed")))]
    print(f"  {len(all_records)} total record(s) (after excluding embargoed/non-public).\n")

    print("STEP 2: Filtering to records with a uniprot_id or pdb_id...")
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

    print(f"Exporting per-simulation metadata to {SIMULATIONS_DIR}/ ...")
    export_metadata(viral)


if __name__ == "__main__":
    main()
