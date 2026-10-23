#!/usr/bin/env python3
"""
convert_mdcath.py

Converts one mdCATH domain's .h5 file (merged_downloads/mdcath/{domain_id}.h5)
into the same {domain_id}_structure.pdb + {domain_id}_trajectory.xtc layout
already used for mdrepo/mddb/dynodb in merged_downloads/, so MD_features.py
and GetContacts can be run on mdCATH domains identically to the other three
sources.

mdCATH schema (confirmed by direct inspection):
    f[domain_id].attrs            -- numChains, numProteinAtoms, numResidues
    f[domain_id]['pdb']           -- scalar bytes, the full solvated system
                                      (protein + water + ions) -- NOT used
                                      here, atom count doesn't match 'coords'
    f[domain_id]['pdbProteinAtoms'] -- scalar bytes, single-model protein-only
                                      reference PDB (the topology), atom
                                      count and order matching 'coords'
                                      exactly (confirmed: both are
                                      numProteinAtoms long); already sliced
                                      to this domain's own atom set -- no
                                      further slicing needed here
    f[domain_id][T][R]['coords']  -- (n_frames, n_atoms, 3) float32, in
                                      Angstrom (confirmed: matches the
                                      Angstrom-scale coordinates embedded
                                      in the 'pdb' text directly)
    f[domain_id][T][R].attrs['numFrames']

T (temperature, K) in {320, 348, 379, 413, 450}; R (replicate) in {0..4}.
No single (T, R) is privileged in the source data -- mdCATH's condition
fields are null in the harmonized schema specifically because there isn't
one representative value (Section on Schema harmonization). To match the
same "closest to 310 K" preference already used everywhere else in this
project's final representative selection, this script defaults to the
closest available temperature (320 K) and the first replicate (0), both
overridable.

Usage:
    python convert_mdcath.py <domain_id>.h5 [--temperature 320] [--replicate 0] [--out-dir DIR]
    python convert_mdcath.py --all [--temperature 320] [--replicate 0]
        # converts every domain in final_representative_ids.json's mdcath list
"""

import argparse
import json
import os

import h5py
import numpy as np
import mdtraj

HERE = os.path.dirname(os.path.abspath(__file__))
MDCATH_DIR = os.path.join(HERE, "..", "analysis", "cross_database", "merged_downloads", "mdcath")
FINAL_IDS_FILE = os.path.join(HERE, "final_representative_ids.json")

ANGSTROM_TO_NM = 0.1
DEFAULT_TEMPERATURE = 320  # closest of {320,348,379,413,450} to the project's 310 K reference
DEFAULT_REPLICATE = 0


def convert_one(h5_path, out_dir, temperature=DEFAULT_TEMPERATURE, replicate=DEFAULT_REPLICATE):
    domain_id = os.path.splitext(os.path.basename(h5_path))[0]
    struct_path = os.path.join(out_dir, f"{domain_id}_structure.pdb")
    traj_path = os.path.join(out_dir, f"{domain_id}_trajectory.xtc")

    if os.path.exists(struct_path) and os.path.exists(traj_path):
        print(f"  skip {domain_id} (already converted)")
        return

    with h5py.File(h5_path, "r") as f:
        g = f[domain_id]
        n_atoms_expected = int(g.attrs["numProteinAtoms"])

        temp_key = str(temperature)
        rep_key = str(replicate)
        if temp_key not in g:
            raise KeyError(f"{domain_id}: temperature {temperature} not present "
                            f"(available: {sorted(g.keys() - {'chain','element','pdb','pdbProteinAtoms','psf','resid','resname','z'})})")
        if rep_key not in g[temp_key]:
            raise KeyError(f"{domain_id}: replicate {replicate} not present at {temperature} K "
                            f"(available: {sorted(g[temp_key].keys())})")

        pdb_bytes = g["pdbProteinAtoms"][()]
        coords_ang = g[temp_key][rep_key]["coords"][:]  # (n_frames, n_atoms, 3), Angstrom

    if coords_ang.shape[1] != n_atoms_expected:
        raise ValueError(f"{domain_id}: coords atom count {coords_ang.shape[1]} != "
                          f"numProteinAtoms {n_atoms_expected}")

    # Write the raw reference structure to a temp path just to load its topology.
    with open(struct_path, "wb") as f_out:
        f_out.write(pdb_bytes)
    top = mdtraj.load(struct_path).topology
    if top.n_atoms != coords_ang.shape[1]:
        raise ValueError(f"{domain_id}: topology has {top.n_atoms} atoms, "
                          f"coords has {coords_ang.shape[1]}")

    # mdCATH domains keep their ORIGINAL residue numbering from the parent
    # PDB chain (e.g. starting at 397, not 1), since a domain is a sliced
    # subregion (see Methods, Sequence reconstruction for mdCATH). But
    # MD_features.py's get_contact_maps() hardcodes
    # `resid = int(atomX.split(':')[2]) - 1` as the 0-based feature-matrix
    # index, which silently assumes contiguous 1-based numbering -- true
    # for mdrepo/mddb/dynodb structures, but not for mdCATH's native
    # numbering (confirmed: residue 397 - 1 = 396 overflows a 64-residue
    # array). Renumbering to a contiguous 1..N range here, once, up front,
    # is far simpler than changing that shared indexing logic, and keeps
    # mdCATH structures usable by the exact same downstream pipeline as
    # the other three sources with no special-casing beyond this script.
    for i, res in enumerate(top.residues):
        res.resSeq = i + 1

    # Build a trajectory from the (now renumbered) topology + this (T, R)'s
    # coordinates, converting Angstrom -> nm for mdtraj's convention.
    traj = mdtraj.Trajectory(xyz=coords_ang.astype(np.float32) * ANGSTROM_TO_NM, topology=top)
    traj[0].save_pdb(struct_path)  # overwrite the raw structure.pdb with the renumbered version
    traj.save_xtc(traj_path)

    print(f"  wrote {domain_id}: {coords_ang.shape[0]} frames, {top.n_atoms} atoms "
          f"({temperature} K, replicate {replicate})")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("h5_file", nargs="?", help="path to one {domain_id}.h5 file")
    p.add_argument("--all", action="store_true",
                    help="convert every domain listed under mdcath in final_representative_ids.json")
    p.add_argument("--temperature", type=int, default=DEFAULT_TEMPERATURE)
    p.add_argument("--replicate", type=int, default=DEFAULT_REPLICATE)
    p.add_argument("--out-dir", default=MDCATH_DIR,
                    help="defaults to merged_downloads/mdcath itself, alongside the source .h5 files")
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    if args.all:
        with open(FINAL_IDS_FILE) as f:
            domain_ids = json.load(f)["mdcath"]
        print(f"Converting {len(domain_ids)} mdCATH domains at {args.temperature} K, replicate {args.replicate}...")
        for domain_id in domain_ids:
            convert_one(os.path.join(MDCATH_DIR, f"{domain_id}.h5"), args.out_dir,
                        args.temperature, args.replicate)
    elif args.h5_file:
        convert_one(args.h5_file, args.out_dir, args.temperature, args.replicate)
    else:
        p.error("provide an h5_file path or --all")


if __name__ == "__main__":
    main()
