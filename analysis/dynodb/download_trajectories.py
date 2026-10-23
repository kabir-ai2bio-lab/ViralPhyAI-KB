#!/usr/bin/env python3
"""
download_trajectories.py

Downloads the real trajectory ensemble files for DynoDB records that won
a slot in the cross-database Final_Unique_Sequences selection (the ones
that actually matter for downstream use), now that the trajectory_frames/
folder has been confirmed real and public (verified directly: files
follow the {PDB_ID}_frame.pdb convention, each a 1,000-frame multi-model
PDB, ~30MB, matching the paper's own stated methodology).

Reads the winner list from final_representative_ids.json's "dynodb" list
(see final_selection.py), rather than from any Excel report.

Usage:
    pip install requests --break-system-packages
    python3 download_trajectories.py
"""

import json
import os
import time

import requests

HERE = os.path.dirname(os.path.abspath(__file__))
FINAL_IDS_FILE = os.path.join(HERE, "..", "..", "features_extraction", "final_representative_ids.json")
TRAJ_URL = "https://huggingface.co/datasets/renzilin/DynoDB/resolve/main/trajectory_frames/{pdb_id}_frame.pdb"
OUT_DIR = "trajectories"
FAILURES_FILE = "trajectory_download_failures.json"

REQUEST_HEADERS = {"User-Agent": "viralplm-analysis-dynodb-trajectories/1.0"}
MAX_RETRIES = 5
RETRY_SLEEP = 3
REQUEST_DELAY = 0.3


def get_final_dynodb_winners():
    with open(FINAL_IDS_FILE) as f:
        final_ids = json.load(f)
    return final_ids.get("dynodb", [])


def download_one(pdb_id):
    """Writes to a .partial temp file and only renames it to the real
    output path on full success, so out_path either doesn't exist or is
    a genuinely complete file -- never a truncated one that a future
    run's os.path.exists() skip-check would mistake for done. This was
    added after a disk-full mid-write left several 0-to-18MB truncated
    files under the real filenames on a previous run."""
    out_path = os.path.join(OUT_DIR, f"{pdb_id}_frame.pdb")
    if os.path.exists(out_path):
        return "skipped"
    tmp_path = out_path + ".partial"
    url = TRAJ_URL.format(pdb_id=pdb_id)
    last_error = None
    for attempt in range(MAX_RETRIES):
        try:
            resp = requests.get(url, headers=REQUEST_HEADERS, timeout=60, stream=True)
            if resp.status_code == 200:
                with open(tmp_path, "wb") as f:
                    for chunk in resp.iter_content(chunk_size=1024 * 1024):
                        f.write(chunk)
                os.rename(tmp_path, out_path)
                return "downloaded"
            if resp.status_code == 404:
                return "404 not found"
            last_error = f"http_{resp.status_code}"
        except Exception as e:
            last_error = str(e)
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
        if attempt < MAX_RETRIES - 1:
            time.sleep(RETRY_SLEEP * (attempt + 1))
    return f"failed after retries: {last_error}"


def main():
    print("Reading DynoDB final-selection winners from final_representative_ids.json...")
    pdb_ids = get_final_dynodb_winners()
    print(f"  {len(pdb_ids)} DynoDB records won a final unique-sequence slot.\n")

    os.makedirs(OUT_DIR, exist_ok=True)
    n_downloaded = n_skipped = 0
    failures = {}

    for i, pdb_id in enumerate(pdb_ids, start=1):
        result = download_one(pdb_id)
        if result == "downloaded":
            n_downloaded += 1
            time.sleep(REQUEST_DELAY)
        elif result == "skipped":
            n_skipped += 1
        else:
            failures[pdb_id] = result

        if i % 25 == 0:
            print(f"  ...{i}/{len(pdb_ids)} ({n_downloaded} downloaded, {n_skipped} already existed, "
                  f"{len(failures)} failed)")

    with open(FAILURES_FILE, "w") as f:
        json.dump(failures, f, indent=2)

    total_size = sum(os.path.getsize(os.path.join(OUT_DIR, f)) for f in os.listdir(OUT_DIR))
    print(f"\nDone: {n_downloaded} downloaded, {n_skipped} already existed, {len(failures)} failed.")
    print(f"Total size of {OUT_DIR}/: {total_size / 1e9:.2f} GB")
    if failures:
        print(f"Failures logged to {FAILURES_FILE}.")


if __name__ == "__main__":
    main()
