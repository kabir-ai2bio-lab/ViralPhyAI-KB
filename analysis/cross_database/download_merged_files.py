#!/usr/bin/env python3
"""
download_merged_files.py

Downloads the real structure + trajectory files for every record in the
Merged_Dataset (combined_records.json, status == "kept" -- every
simulation that survived dedup.py, not just the one-per-sequence final
winners). All four sources host real files, confirmed by direct
inspection before writing this:

  - mdrepo:  https://assets.mdrepo.org/{slug}/{filename} -- a fully
             public static host, no login needed (verified with a bare
             curl, no cookies, even though the site's own API requires
             a logged-in session for its create_download_instance/zip
             flow -- that flow turned out not to be necessary). Only
             minimal.pdb (solvent-stripped structure) and sampled.xtc
             (downsampled trajectory) are actually hosted here --
             confirmed by testing 5 real kept records: minimal.xtc and
             full.xtc/full.pdb 404 on every single one, while
             minimal.pdb and sampled.xtc 200 on every single one. This
             matches exactly what mdrepo's own site fetches for its
             structure/trajectory preview viewer, and is presumably the
             full extent of what's pushed to this public CDN -- the
             complete (non-downsampled) trajectory stays behind the
             authenticated zip-download flow. A 404 is still logged as
             a per-file miss rather than fatal, in case some records
             genuinely differ.
  - mddb:    https://mdposit.mddbr.eu/api/rest/current/projects/{accession}/files/{filename}
             redirects to the project's actual storage node; public,
             requests' default redirect-following handles it.
             structure.pdb + trajectory.xtc, confirmed present in the
             `files` list of all 2,353 collected records.
  - mdcath:  one H5 file per domain on HuggingFace
             (data/mdcath_dataset_{domain_id}.h5), containing both the
             reference structure and all 5 replicates x 5 temperatures
             together -- there's no separate structure file to fetch.
             Public, no auth. Only 11 kept records, but files run
             hundreds of MB to ~1.6GB each.
  - dynodb:  structure is already embedded in the metadata collected by
             dynodb/collect_metadata.py (static_structure_pdb field) --
             written directly here, no network call. Trajectory from
             trajectory_frames/{pdb_id}_frame.pdb, same source
             download_trajectories.py already uses -- that script only
             covered the 324 Final_Unique_Sequences winners, this
             covers all 370 kept (Merged_Dataset) dynodb records.

Every file write is atomic (temp file + os.rename on success only), so
a crash or disk-full mid-download never leaves a truncated file under
the real filename for a future run's skip-check to mistake as done.

Output layout:
    merged_downloads/{dataset}/{record_id}_structure.{ext}
    merged_downloads/{dataset}/{record_id}_trajectory.{ext}

Usage:
    pip install requests --break-system-packages
    python3 download_merged_files.py
"""

import json
import os
import time

import requests

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # data/analysis
COMBINED_FILE = "combined_records.json"
OUT_DIR = "merged_downloads"
FAILURES_FILE = "merged_download_failures.json"

REQUEST_HEADERS = {"User-Agent": "viralplm-analysis-cross_database-merged-downloads/1.0"}
MAX_RETRIES = 5
RETRY_SLEEP = 3
REQUEST_DELAY = 0.3

MDREPO_ASSET_URL = "https://assets.mdrepo.org/{slug}/{filename}"
MDDB_FILE_URL = "https://mdposit.mddbr.eu/api/rest/current/projects/{accession}/files/{filename}"
MDCATH_H5_URL = "https://huggingface.co/datasets/compsciencelab/mdCATH/resolve/main/data/mdcath_dataset_{domain_id}.h5"
DYNODB_TRAJ_URL = "https://huggingface.co/datasets/renzilin/DynoDB/resolve/main/trajectory_frames/{pdb_id}_frame.pdb"


def load_kept_records(dataset):
    with open(COMBINED_FILE) as f:
        combined = json.load(f)
    return [r for r in combined if r["dataset"] == dataset and r["status"] == "kept"]


def download_url(url, out_path, treat_404_as_skip=True):
    """Atomic GET -> file. Returns "downloaded" / "skipped" / "404" /
    "failed: <reason>". Retries on anything except a definitive 404."""
    if os.path.exists(out_path):
        return "skipped"
    tmp_path = out_path + ".partial"
    last_error = None
    for attempt in range(MAX_RETRIES):
        try:
            resp = requests.get(url, headers=REQUEST_HEADERS, timeout=120, stream=True)
            if resp.status_code == 200:
                with open(tmp_path, "wb") as f:
                    for chunk in resp.iter_content(chunk_size=1024 * 1024):
                        f.write(chunk)
                os.rename(tmp_path, out_path)
                return "downloaded"
            if resp.status_code == 404 and treat_404_as_skip:
                return "404"
            last_error = f"http_{resp.status_code}"
        except Exception as e:
            last_error = str(e)
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
        if attempt < MAX_RETRIES - 1:
            time.sleep(RETRY_SLEEP * (attempt + 1))
    return f"failed: {last_error}"


def run_mdrepo(out_dir, failures):
    records = load_kept_records("mdrepo")
    print(f"mdrepo: {len(records)} kept records")
    n_struct = n_traj = 0
    for i, r in enumerate(records, start=1):
        slug = r["record_id"]
        struct_result = download_url(
            MDREPO_ASSET_URL.format(slug=slug, filename="minimal.pdb"),
            os.path.join(out_dir, f"{slug}_structure.pdb"))
        traj_result = download_url(
            MDREPO_ASSET_URL.format(slug=slug, filename="sampled.xtc"),
            os.path.join(out_dir, f"{slug}_trajectory.xtc"))
        if struct_result == "downloaded":
            n_struct += 1
        if traj_result == "downloaded":
            n_traj += 1
        if struct_result not in ("downloaded", "skipped") or traj_result not in ("downloaded", "skipped"):
            failures.setdefault("mdrepo", {})[slug] = {"structure": struct_result, "trajectory": traj_result}
        time.sleep(REQUEST_DELAY)
        if i % 100 == 0:
            print(f"  ...{i}/{len(records)} ({n_struct} structures, {n_traj} trajectories downloaded)")
    print(f"  done: {n_struct} structures, {n_traj} trajectories newly downloaded.\n")


def run_mddb(out_dir, failures):
    records = load_kept_records("mddb")
    print(f"mddb: {len(records)} kept records")
    n_struct = n_traj = 0
    for i, r in enumerate(records, start=1):
        accession = r["record_id"]
        struct_result = download_url(
            MDDB_FILE_URL.format(accession=accession, filename="structure.pdb"),
            os.path.join(out_dir, f"{accession}_structure.pdb"))
        traj_result = download_url(
            MDDB_FILE_URL.format(accession=accession, filename="trajectory.xtc"),
            os.path.join(out_dir, f"{accession}_trajectory.xtc"))
        if struct_result == "downloaded":
            n_struct += 1
        if traj_result == "downloaded":
            n_traj += 1
        if struct_result not in ("downloaded", "skipped") or traj_result not in ("downloaded", "skipped"):
            failures.setdefault("mddb", {})[accession] = {"structure": struct_result, "trajectory": traj_result}
        time.sleep(REQUEST_DELAY)
        if i % 100 == 0:
            print(f"  ...{i}/{len(records)} ({n_struct} structures, {n_traj} trajectories downloaded)")
    print(f"  done: {n_struct} structures, {n_traj} trajectories newly downloaded.\n")


def run_mdcath(out_dir, failures):
    records = load_kept_records("mdcath")
    print(f"mdcath: {len(records)} kept records (large files -- this may take a while)")
    n = 0
    for i, r in enumerate(records, start=1):
        domain_id = r["record_id"]
        result = download_url(
            MDCATH_H5_URL.format(domain_id=domain_id),
            os.path.join(out_dir, f"{domain_id}.h5"))
        if result == "downloaded":
            n += 1
            print(f"  [{i}/{len(records)}] downloaded {domain_id}.h5")
        elif result != "skipped":
            failures.setdefault("mdcath", {})[domain_id] = result
        time.sleep(REQUEST_DELAY)
    print(f"  done: {n} domain files newly downloaded.\n")


def run_dynodb(out_dir, failures):
    records = load_kept_records("dynodb")
    print(f"dynodb: {len(records)} kept records")

    # Structure is already embedded in the collected metadata -- index
    # the simulations/ folder once by pdb_id to avoid re-downloading it.
    sim_dir = os.path.join(BASE, "dynodb", "simulations")
    by_pdb_id = {}
    for fn in os.listdir(sim_dir):
        pdb_id = fn.split("_", 1)[0]
        by_pdb_id[pdb_id] = fn

    n_struct = n_traj = 0
    for i, r in enumerate(records, start=1):
        pdb_id = r["record_id"]
        struct_path = os.path.join(out_dir, f"{pdb_id}_structure.pdb")
        if not os.path.exists(struct_path):
            fn = by_pdb_id.get(pdb_id)
            if fn:
                rec = json.load(open(os.path.join(sim_dir, fn)))
                pdb_text = rec.get("static_structure_pdb")
                if pdb_text:
                    tmp_path = struct_path + ".partial"
                    with open(tmp_path, "w") as f:
                        f.write(pdb_text)
                    os.rename(tmp_path, struct_path)
                    n_struct += 1
                else:
                    failures.setdefault("dynodb", {}).setdefault(pdb_id, {})["structure"] = "no static_structure_pdb in collected metadata"
            else:
                failures.setdefault("dynodb", {}).setdefault(pdb_id, {})["structure"] = "record not found in simulations/"

        traj_result = download_url(
            DYNODB_TRAJ_URL.format(pdb_id=pdb_id),
            os.path.join(out_dir, f"{pdb_id}_trajectory.pdb"))
        if traj_result == "downloaded":
            n_traj += 1
        elif traj_result != "skipped":
            failures.setdefault("dynodb", {}).setdefault(pdb_id, {})["trajectory"] = traj_result
        time.sleep(REQUEST_DELAY)
        if i % 50 == 0:
            print(f"  ...{i}/{len(records)} ({n_struct} structures written, {n_traj} trajectories downloaded)")
    print(f"  done: {n_struct} structures written, {n_traj} trajectories newly downloaded.\n")


def main():
    failures = {}
    for dataset, runner in [("mdrepo", run_mdrepo), ("mddb", run_mddb),
                             ("mdcath", run_mdcath), ("dynodb", run_dynodb)]:
        out_dir = os.path.join(OUT_DIR, dataset)
        os.makedirs(out_dir, exist_ok=True)
        runner(out_dir, failures)

    with open(FAILURES_FILE, "w") as f:
        json.dump(failures, f, indent=2)

    total_failures = sum(len(v) for v in failures.values())
    if total_failures:
        print(f"{total_failures} record(s) had at least one missing/failed file -- see {FAILURES_FILE}.")
    else:
        print("No failures.")


if __name__ == "__main__":
    main()
