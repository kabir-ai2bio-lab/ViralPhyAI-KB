#!/usr/bin/env python3
"""
run_feature_extraction.py

Batch driver: runs the two-step feature extraction pipeline (GetContacts ->
MD_features.py) over the 2,110 Final Unique Sequences representatives
listed in final_representative_ids.json.

By this point every source has been normalized to the same
{record_id}_structure.pdb + _trajectory.xtc layout in
merged_downloads/{source}/ (mdcath via convert_mdcath.py, dynodb via
convert_dynodb.py). But residue numbering is NOT reliably contiguous
from 1 even within mdrepo/mddb themselves -- confirmed by a failure
during testing: MDR00027083 has 92 contiguous residues numbered 5..96
(not 1..92), which overflows MD_features.py's get_contact_maps(),
whose `resid - 1` indexing assumes 1-based numbering. So this driver
renumbers every record's topology to a contiguous 1..N range itself,
writing a renumbered structure PDB to a local staging file before
calling GetContacts/MD_features.py, regardless of source -- the same
fix already applied once via convert_mdcath.py/convert_dynodb.py, just
applied universally here since no source turned out to be reliably
exempt.

Resumable: skips any record whose two output files already exist. Never
aborts the whole batch on one record's failure -- logs it and continues.

Usage:
    python run_feature_extraction.py [--sources mdrepo,mddb,dynodb,mdcath]
                                      [--limit N] [--keep-contacts]
                                      [--cores N]
"""

import argparse
import json
import os
import subprocess
import sys
import time

import mdtraj

HERE = os.path.dirname(os.path.abspath(__file__))
MERGED_DIR = os.path.join(HERE, "..", "analysis", "cross_database", "merged_downloads")
FINAL_IDS_FILE = os.path.join(HERE, "final_representative_ids.json")
FEATURES_DIR = os.path.join(HERE, "features")
GETCONTACTS_SCRIPT = os.path.join(HERE, "getcontacts", "get_dynamic_contacts.py")
MD_FEATURES_SCRIPT = os.path.join(HERE, "MD_features.py")
FAILURES_FILE = os.path.join(HERE, "run_feature_extraction_failures.json")

ITYPES = ["vdw", "hb", "hp", "sb", "pc", "ps", "ts"]


def process_one(source, record_id, cores, keep_contacts):
    src_dir = os.path.join(MERGED_DIR, source)
    struct_path = os.path.join(src_dir, f"{record_id}_structure.pdb")
    traj_path = os.path.join(src_dir, f"{record_id}_trajectory.xtc")

    out_dir = os.path.join(FEATURES_DIR, source)
    os.makedirs(out_dir, exist_ok=True)
    out_prefix = os.path.join(out_dir, record_id)
    res_out = f"{out_prefix}_res_feature.csv"
    pair_out = f"{out_prefix}_pair_feature.npz"

    if os.path.exists(res_out) and os.path.exists(pair_out):
        return "skipped"

    if not (os.path.exists(struct_path) and os.path.exists(traj_path)):
        raise FileNotFoundError(f"missing structure/trajectory: {struct_path} / {traj_path}")

    contacts_path = f"{out_prefix}_contacts.tsv"
    renumbered_struct_path = f"{out_prefix}_structure_renumbered.pdb"

    # Renumber residues to a contiguous 1..N range -- required by
    # MD_features.py's `resid - 1` indexing, and not reliably true of the
    # source structures themselves (see module docstring). Renumbering the
    # topology doesn't touch coordinates/atom order, so the original
    # trajectory.xtc is reused unchanged against this renumbered topology.
    ref = mdtraj.load_frame(traj_path, 0, top=struct_path)
    for i, res in enumerate(ref.topology.residues):
        res.resSeq = i + 1
    ref.save_pdb(renumbered_struct_path)

    # Step 1: GetContacts. --cores 1 avoids a multiprocessing/semaphore
    # crash observed in one sandboxed environment during testing; raise
    # this if you confirm parallel cores work fine in your own shell.
    r1 = subprocess.run(
        [sys.executable, GETCONTACTS_SCRIPT,
         "--topology", renumbered_struct_path, "--trajectory", traj_path,
         "--itypes", *ITYPES, "--cores", str(cores),
         "--output", contacts_path],
        capture_output=True, text=True,
    )
    if r1.returncode != 0:
        raise RuntimeError(f"get_dynamic_contacts.py failed: {r1.stderr[-2000:]}")

    # Step 2: MD_features.py (same renumbered topology as Step 1, so the
    # contacts.tsv resid values line up with res_feature.csv's row order)
    r2 = subprocess.run(
        [sys.executable, MD_FEATURES_SCRIPT,
         "-p", renumbered_struct_path, "-t", traj_path, "-i", contacts_path, "-o", out_prefix],
        capture_output=True, text=True,
    )
    if r2.returncode != 0:
        raise RuntimeError(f"MD_features.py failed: {r2.stderr[-2000:]}")

    if not keep_contacts:
        for p in (contacts_path, renumbered_struct_path):
            try:
                os.remove(p)
            except OSError:
                pass

    return "processed"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sources", default="mdrepo,mddb,dynodb,mdcath")
    p.add_argument("--limit", type=int, default=None,
                    help="stop after this many records (across all sources combined), for pilot runs")
    p.add_argument("--cores", type=int, default=1)
    p.add_argument("--keep-contacts", action="store_true",
                    help="keep the intermediate *_contacts.tsv files (large; deleted by default)")
    p.add_argument("--skip-failed", action="store_true",
                    help="skip records already listed in the failures file instead of retrying "
                         "them; their failure entries are left untouched")
    args = p.parse_args()

    with open(FINAL_IDS_FILE) as f:
        all_ids = json.load(f)

    failures = {}
    if os.path.exists(FAILURES_FILE):
        with open(FAILURES_FILE) as f:
            failures = json.load(f)

    sources = args.sources.split(",")
    work = [(src, rid) for src in sources for rid in all_ids.get(src, [])]
    if args.skip_failed:
        n_before = len(work)
        work = [(src, rid) for src, rid in work if f"{src}/{rid}" not in failures]
        print(f"--skip-failed: leaving out {n_before - len(work)} previously failed record(s)")
    if args.limit:
        work = work[: args.limit]

    print(f"Processing {len(work)} record(s) across sources: {sources}")

    n_processed = n_skipped = n_failed = 0
    t0 = time.time()
    for i, (source, record_id) in enumerate(work, 1):
        try:
            status = process_one(source, record_id, args.cores, args.keep_contacts)
            if status == "skipped":
                n_skipped += 1
            else:
                n_processed += 1
                failures.pop(f"{source}/{record_id}", None)
        except Exception as e:
            n_failed += 1
            failures[f"{source}/{record_id}"] = str(e)
            print(f"  FAILED {source}/{record_id}: {e}")

        if i % 10 == 0 or i == len(work):
            elapsed = time.time() - t0
            print(f"  ...{i}/{len(work)} (processed={n_processed} skipped={n_skipped} "
                  f"failed={n_failed}, {elapsed:.0f}s elapsed)")
            with open(FAILURES_FILE, "w") as f:
                json.dump(failures, f, indent=2)

    with open(FAILURES_FILE, "w") as f:
        json.dump(failures, f, indent=2)

    print(f"\nDone. processed={n_processed} skipped={n_skipped} failed={n_failed}")
    if failures:
        print(f"{len(failures)} failure(s) -- see {FAILURES_FILE}")


if __name__ == "__main__":
    main()
