#!/usr/bin/env python3
"""
collect_metadata.py

Comprehensive-review workflow for the mdCATH dataset
(https://huggingface.co/datasets/compsciencelab/mdCATH), following the
direction agreed on in chat: identify viral domains WITHOUT downloading
any of the 3.61 TB of actual trajectory data, using only:

  1. mdCATH's own master domain list (mdcath_domains.txt, ~5,398 CATH
     domain IDs, e.g. "3oa1B00" = PDB 3oa1, chain B, domain 00).
  2. CATH's official domain-boundary file (cath-domain-boundaries-
     seqreschopping.txt) to get each domain's residue range within its
     PDB chain.
  3. RCSB's polymer_entity API to get, per PDB chain: the full SEQRES
     sequence, the source organism (with taxonomy_lineage -- "Viruses"
     appears directly in there, no separate NCBI call needed, unlike
     the mdrepo/mddb pipelines), and a UniProt cross-reference (SIFTS).

A domain's "sequence" is the SEQRES sequence sliced at its CATH residue
range -- i.e. the fragment actually simulated -- NOT the full UniProt
canonical protein. This matches how mddb's metadata.SEQUENCES and
mdrepo's fasta_sequence both work (verified byte-for-byte against a real
overlapping example, MD-A00007 / MDR00056799 / mdCATH domain 3oa1B00,
before this script was written).

  STEP 1. Fetch the master domain list.
  STEP 2. Fetch and parse the CATH boundary file.
  STEP 3. For each distinct PDB code referenced, fetch its polymer
          entities once (cached) -- sequence, organism, UniProt xref.
  STEP 4. For each domain: resolve its entity by chain, verify viral via
          taxonomy_lineage, slice the sequence at its boundary range.
  (new)   Write one file per confirmed-viral domain to
          simulations/{domain_id}_{uniprot_id_or_pdb_id}.json

Only the domain list and the boundary file are ever downloaded in full;
everything else is one small JSON call per distinct PDB code (there are
far fewer distinct PDB codes than domains, since many domains share a
PDB entry). No H5 files are ever touched by this script.

No heteromeric-complex filter is applied here (unlike mdrepo/mddb/dynodb's
collect_metadata.py STEP4): a CATH domain is defined as a segment of a
SINGLE PDB chain by construction (see domain_id format above -- one chain
letter per domain, e.g. "B" in "3oa1B00"), so every domain simulated by
mdCATH is inherently single-chain. There is no multi-chain case to exclude.

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
DOMAIN_LIST_URL = "https://raw.githubusercontent.com/compsciencelab/mdCATH/main/mdcath_domains.txt"
BOUNDARY_FILE_URL = ("http://download.cathdb.info/cath/releases/latest-release/"
                      "cath-classification-data/cath-domain-boundaries-seqreschopping.txt")
RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
RCSB_ENTITY_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"

DOMAIN_LIST_FILE = "mdcath_domains.txt"
BOUNDARY_FILE = "cath-domain-boundaries-seqreschopping.txt"
ENTITY_CACHE_FILE = "rcsb_entity_cache.json"  # pdb_id -> list of entity dicts (raw + derived)
SIMULATIONS_DIR = "simulations"
FAILURES_FILE = "mdcath_resolve_failures.json"  # domain_id -> reason

REQUEST_HEADERS = {"User-Agent": "viralplm-analysis-mdcath/1.0"}
MAX_RETRIES = 5
RETRY_SLEEP = 3
REQUEST_DELAY = 0.3  # throttle between RCSB calls -- same defensive default as the
                      # mdrepo/mddb analysis scripts, after mdrepo's server was found
                      # to rate-limit unthrottled bursts.


def _get_with_retry(url: str, timeout: int = 20) -> Optional[requests.Response]:
    for attempt in range(MAX_RETRIES):
        try:
            resp = requests.get(url, headers=REQUEST_HEADERS, timeout=timeout)
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
# STEP 1: Master domain list
# ---------------------------------------------------------------------
def fetch_domain_list() -> list:
    if os.path.exists(DOMAIN_LIST_FILE):
        print(f"  Loading cached domain list from {DOMAIN_LIST_FILE}...")
    else:
        print(f"  Downloading domain list from {DOMAIN_LIST_URL} ...")
        resp = _get_with_retry(DOMAIN_LIST_URL, timeout=30)
        if resp is None or resp.status_code != 200:
            raise RuntimeError("failed to fetch mdcath_domains.txt after retries")
        with open(DOMAIN_LIST_FILE, "w") as f:
            f.write(resp.text)

    with open(DOMAIN_LIST_FILE) as f:
        domains = [line.strip() for line in f if line.strip()]
    return domains


# ---------------------------------------------------------------------
# STEP 2: CATH domain boundary file
# ---------------------------------------------------------------------
def fetch_boundary_map() -> dict:
    if os.path.exists(BOUNDARY_FILE):
        print(f"  Loading cached boundary file from {BOUNDARY_FILE}...")
    else:
        print(f"  Downloading boundary file from {BOUNDARY_FILE_URL} ...")
        for attempt in range(MAX_RETRIES):
            try:
                with requests.get(BOUNDARY_FILE_URL, headers=REQUEST_HEADERS,
                                   timeout=120, stream=True) as resp:
                    resp.raise_for_status()
                    with open(BOUNDARY_FILE, "wb") as f:
                        for chunk in resp.iter_content(chunk_size=1024 * 1024):
                            f.write(chunk)
                break
            except Exception as e:
                if attempt == MAX_RETRIES - 1:
                    raise
                print(f"    [warn] boundary file download failed ({e}), retrying...")
                time.sleep(RETRY_SLEEP)

    boundary_map = {}
    with open(BOUNDARY_FILE) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) != 2:
                continue
            domain_id, ranges_str = parts
            segments = []
            for seg in ranges_str.split(","):
                try:
                    start, end = seg.split("-")
                    segments.append((int(start), int(end)))
                except ValueError:
                    segments = None
                    break
            if segments:
                boundary_map[domain_id] = segments
    return boundary_map


def slice_sequence(full_sequence: str, segments: list) -> Optional[str]:
    """Segments are 1-indexed, inclusive, straightforward SEQRES numbering
    (verified against a real domain before this script was written)."""
    pieces = []
    for start, end in segments:
        if start < 1 or end > len(full_sequence) or start > end:
            return None
        pieces.append(full_sequence[start - 1:end])
    return "".join(pieces)


# ---------------------------------------------------------------------
# STEP 3: Per-PDB-code entity resolution (cached; shared across domains)
# ---------------------------------------------------------------------
def load_entity_cache() -> dict:
    if os.path.exists(ENTITY_CACHE_FILE):
        with open(ENTITY_CACHE_FILE) as f:
            return json.load(f)
    return {}


def save_entity_cache(cache: dict):
    with open(ENTITY_CACHE_FILE, "w") as f:
        json.dump(cache, f)


RCSB_REMOVED_URL = "https://data.rcsb.org/rest/v1/holdings/removed/{pdb_id}"


def resolve_current_pdb_id(pdb_id: str, max_hops: int = 3) -> Optional[str]:
    """CATH's boundary file references PDB IDs from whatever CATH release
    built mdCATH; RCSB's live API only serves currently-active entries, so
    an obsoleted/superseded/merged entry 404s. Follow RCSB's own
    'removed -> id_codes_replaced_by' chain to find the current id.
    Returns None if the entry has no successor (a true dead end, e.g. a
    retracted paper) or the redirect chain is too long."""
    current = pdb_id
    for _ in range(max_hops):
        resp = _get_with_retry(RCSB_ENTRY_URL.format(pdb_id=current))
        if resp is not None and resp.status_code == 200:
            return current
        removed_resp = _get_with_retry(RCSB_REMOVED_URL.format(pdb_id=current.upper()))
        if removed_resp is None or removed_resp.status_code != 200:
            return None
        replaced_by = (removed_resp.json().get("rcsb_repository_holdings_removed") or {}).get(
            "id_codes_replaced_by")
        if not replaced_by:
            return None
        current = replaced_by[0]  # follow the first successor; log the hop via the caller
    return None


def fetch_pdb_entities(pdb_id: str) -> list:
    """Returns a list of entity dicts for this PDB code: each with
    auth_asym_ids, sequence, organism info (incl. taxonomy_lineage),
    uniprot_id, the raw RCSB response, and resolved_pdb_id (the id
    entities were actually fetched from -- may differ from pdb_id if the
    original entry was obsolete and a successor was followed). Empty list
    on failure."""
    entry_resp = _get_with_retry(RCSB_ENTRY_URL.format(pdb_id=pdb_id))
    resolved_pdb_id = pdb_id
    if entry_resp is None or entry_resp.status_code != 200:
        resolved_pdb_id = resolve_current_pdb_id(pdb_id)
        if resolved_pdb_id is None:
            return []
        entry_resp = _get_with_retry(RCSB_ENTRY_URL.format(pdb_id=resolved_pdb_id))
        if entry_resp is None or entry_resp.status_code != 200:
            return []
    entry_data = entry_resp.json()
    entity_ids = entry_data.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids") or []

    entities = []
    for eid in entity_ids:
        time.sleep(REQUEST_DELAY)
        resp = _get_with_retry(RCSB_ENTITY_URL.format(pdb_id=resolved_pdb_id, entity_id=eid))
        if resp is None or resp.status_code != 200:
            continue
        raw = resp.json()

        auth_asym_ids = (raw.get("rcsb_polymer_entity_container_identifiers") or {}).get("auth_asym_ids") or []
        sequence = (raw.get("entity_poly") or {}).get("pdbx_seq_one_letter_code_can", "")

        organisms = raw.get("rcsb_entity_source_organism") or []
        is_viral = False
        organism_name = None
        for org in organisms:
            lineage_names = {seg.get("name") for seg in (org.get("taxonomy_lineage") or [])}
            if "Viruses" in lineage_names:
                is_viral = True
                organism_name = org.get("scientific_name") or org.get("ncbi_scientific_name")
                break
        if not organism_name and organisms:
            organism_name = organisms[0].get("scientific_name") or organisms[0].get("ncbi_scientific_name")

        uniprot_id = None
        for align in (raw.get("rcsb_polymer_entity_align") or []):
            if align.get("reference_database_name") == "UniProt":
                uniprot_id = align.get("reference_database_accession")
                break

        entities.append({
            "entity_id": eid,
            "auth_asym_ids": auth_asym_ids,
            "sequence": sequence,
            "organism_name": organism_name,
            "is_viral": is_viral,
            "uniprot_id": uniprot_id,
            "raw": raw,
            "resolved_pdb_id": resolved_pdb_id,
            "redirected": resolved_pdb_id != pdb_id,
        })
    return entities


def get_entities_cached(pdb_id: str, cache: dict) -> list:
    if pdb_id in cache:
        return cache[pdb_id]
    entities = fetch_pdb_entities(pdb_id)
    cache[pdb_id] = entities
    return entities


# ---------------------------------------------------------------------
# STEP 4: Resolve each domain and export confirmed-viral ones
# ---------------------------------------------------------------------
def parse_domain_id(domain_id: str) -> tuple:
    """'3oa1B00' -> ('3oa1', 'B', '00')"""
    return domain_id[:4], domain_id[4], domain_id[5:]


def id_slug(domain_id: str, uniprot_id: Optional[str], pdb_id: str) -> str:
    identifier = uniprot_id or pdb_id
    return f"{domain_id}_{identifier}"


def main():
    print("STEP 1: Fetching mdCATH master domain list...")
    domains = fetch_domain_list()
    print(f"  {len(domains)} domain(s).\n")

    print("STEP 2: Fetching CATH domain boundary file...")
    boundary_map = fetch_boundary_map()
    print(f"  {len(boundary_map)} domain boundary entries loaded.\n")

    print("STEP 3-4: Resolving each domain (sequence + viral organism check)...")
    entity_cache = load_entity_cache()
    os.makedirs(SIMULATIONS_DIR, exist_ok=True)

    n_written = 0
    n_skipped = 0
    n_not_viral = 0
    failures = {}
    pdb_ids_seen = set()

    for i, domain_id in enumerate(domains, start=1):
        pdb_id, chain, _num = parse_domain_id(domain_id)

        segments = boundary_map.get(domain_id)
        if not segments:
            failures[domain_id] = "no boundary entry found"
            continue

        if pdb_id not in pdb_ids_seen:
            pdb_ids_seen.add(pdb_id)
            if pdb_id not in entity_cache:
                time.sleep(REQUEST_DELAY)

        entities = get_entities_cached(pdb_id, entity_cache)
        if not entities:
            failures[domain_id] = "could not resolve PDB entry/entities"
            continue

        entity = next((e for e in entities if chain in e["auth_asym_ids"]), None)
        if entity is None:
            failures[domain_id] = f"no entity contains chain {chain}"
            continue

        if not entity["is_viral"]:
            n_not_viral += 1
            continue

        slug = id_slug(domain_id, entity["uniprot_id"], pdb_id)
        out_path = os.path.join(SIMULATIONS_DIR, f"{slug}.json")
        if os.path.exists(out_path):
            n_skipped += 1
            continue

        sequence = slice_sequence(entity["sequence"], segments)
        if sequence is None:
            failures[domain_id] = "boundary range out of bounds for entity sequence"
            continue

        merged = {
            "domain_id": domain_id,
            "pdb_id": pdb_id,
            "chain": chain,
            "cath_residue_range": segments,
            "sequence": sequence,
            "uniprot_id": entity["uniprot_id"],
            "organism": entity["organism_name"],
            "resolved_pdb_id": entity["resolved_pdb_id"],
            "redirected_from_obsolete_entry": entity["redirected"],
            "detail": entity["raw"],
        }
        with open(out_path, "w") as f:
            json.dump(merged, f, indent=2)
        n_written += 1

        if i % 200 == 0:
            print(f"    ...{i}/{len(domains)} domains processed "
                  f"({n_written} written, {n_not_viral} not viral, "
                  f"{n_skipped} already existed, {len(failures)} failed), "
                  f"{len(pdb_ids_seen)} distinct PDB codes seen")
            save_entity_cache(entity_cache)
            with open(FAILURES_FILE, "w") as f:
                json.dump(failures, f, indent=2)

    save_entity_cache(entity_cache)
    with open(FAILURES_FILE, "w") as f:
        json.dump(failures, f, indent=2)

    print(f"\n  Done: {n_written} written, {n_not_viral} confirmed non-viral, "
          f"{n_skipped} already existed, {len(failures)} failed to resolve, "
          f"out of {len(domains)} total domains ({len(pdb_ids_seen)} distinct PDB codes).")
    if failures:
        print(f"  Failure reasons logged to {FAILURES_FILE}.")


if __name__ == "__main__":
    main()
