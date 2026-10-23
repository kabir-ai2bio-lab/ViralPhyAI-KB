#!/usr/bin/env python3
"""
convert_dynodb.py

Normalizes one DynoDB record's downloaded files (merged_downloads/dynodb/
{pdb_id}_structure.pdb + {pdb_id}_trajectory.pdb) into the same
{pdb_id}_structure.pdb + {pdb_id}_trajectory.xtc layout used by the other
three sources, so MD_features.py and GetContacts can be run on DynoDB
identically to mdrepo/mddb/mdcath.

Two problems found by direct inspection (not documented anywhere):

1. structure.pdb and trajectory.pdb do NOT describe the same atom set.
   For 7LV2: structure.pdb has 40 atoms, trajectory.pdb has 70 atoms per
   frame -- loading trajectory.pdb with structure.pdb as topology fails
   outright ("xyz must be shape (Any, 40, 3). You supplied (1000, 70, 3)").
   structure.pdb is presumably a different processing pass (e.g. the
   deposited PDB's own model) than whatever the trajectory was actually
   run on. trajectory.pdb, by contrast, is self-consistent: mdtraj loads
   it standalone fine (1000 frames, one shared topology). So this script
   ignores structure.pdb entirely and derives topology + coordinates both
   from trajectory.pdb.

2. Residue numbering is native to the source PDB entry, not renumbered
   from 1 (e.g. 7LV2 runs 179-184, a 6-residue fragment of a larger
   deposited structure -- consistent with DynoDB targeting small peptide
   fragments, Section on Description of the Data). Same fix as mdCATH:
   MD_features.py's get_contact_maps() assumes contiguous 1-based
   numbering, so residues are renumbered 1..N here, once, up front.

Usage:
    python convert_dynodb.py <pdb_id> [--out-dir DIR]
    python convert_dynodb.py --all
        # converts every record listed under dynodb in final_representative_ids.json
"""

import argparse
import json
import os

import numpy as np
import mdtraj

HERE = os.path.dirname(os.path.abspath(__file__))
DYNODB_DIR = os.path.join(HERE, "..", "analysis", "cross_database", "merged_downloads", "dynodb")
FINAL_IDS_FILE = os.path.join(HERE, "final_representative_ids.json")

NORMALIZED_SUFFIX = "_normalized"  # marker file to avoid reconverting


def convert_one(pdb_id, out_dir):
    struct_path = os.path.join(out_dir, f"{pdb_id}_structure.pdb")
    traj_path = os.path.join(out_dir, f"{pdb_id}_trajectory.pdb")
    marker_path = os.path.join(out_dir, f"{pdb_id}{NORMALIZED_SUFFIX}")
    new_traj_path = os.path.join(out_dir, f"{pdb_id}_trajectory.xtc")

    if os.path.exists(marker_path) and os.path.exists(new_traj_path):
        print(f"  skip {pdb_id} (already normalized)")
        return

    if not os.path.exists(traj_path):
        raise FileNotFoundError(f"{pdb_id}: no trajectory.pdb found at {traj_path}")

    traj = mdtraj.load(traj_path)  # self-describing multi-model PDB
    top = traj.topology

    for i, res in enumerate(top.residues):
        res.resSeq = i + 1

    traj.save_xtc(new_traj_path)
    traj[0].save_pdb(struct_path)  # overwrite the mismatched original structure.pdb
    with open(marker_path, "w") as f:
        f.write("ok")

    print(f"  wrote {pdb_id}: {traj.n_frames} frames, {top.n_atoms} atoms, {top.n_residues} residues")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("pdb_id", nargs="?")
    p.add_argument("--all", action="store_true")
    p.add_argument("--out-dir", default=DYNODB_DIR)
    args = p.parse_args()

    if args.all:
        with open(FINAL_IDS_FILE) as f:
            pdb_ids = json.load(f)["dynodb"]
        print(f"Normalizing {len(pdb_ids)} DynoDB records...")
        failures = {}
        for pdb_id in pdb_ids:
            try:
                convert_one(pdb_id, args.out_dir)
            except Exception as e:
                failures[pdb_id] = str(e)
                print(f"  FAILED {pdb_id}: {e}")
        if failures:
            fail_path = os.path.join(args.out_dir, "convert_dynodb_failures.json")
            with open(fail_path, "w") as f:
                json.dump(failures, f, indent=2)
            print(f"{len(failures)} failure(s) -- see {fail_path}")
        else:
            print("No failures.")
    elif args.pdb_id:
        convert_one(args.pdb_id, args.out_dir)
    else:
        p.error("provide a pdb_id or --all")


if __name__ == "__main__":
    main()
