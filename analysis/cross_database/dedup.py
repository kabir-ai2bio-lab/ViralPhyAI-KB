#!/usr/bin/env python3
"""
dedup.py

Cross-database redundancy resolution across all four collected sources
(mdrepo, mddb, mdcath, dynodb). Selects unique SIMULATIONS, not unique
sequences -- a sequence with several genuinely different simulations
(different temperature, force field, water model, or ligand state)
keeps all of them. Reducing further to one representative per unique
sequence is a separate, later step (build_excel.py's
select_final_representative / compute_final_selection), out of scope
here.

Rule, applied in a single pass per exact-sequence group:

  - Grouping key is EXACT sequence string equality -- NOT CD-HIT's
    "100%" clustering, which can still merge a full sequence with a
    shorter fragment fully contained in it (that's why e.g. mdrepo went
    from 1,346 exact-unique strings down to 635 CD-HIT clusters
    earlier). Exact string equality is the only grouping that reliably
    keeps "different part of the sequence" cases separate, as required.

  - Within an exact-sequence group, two records are considered "the
    same simulation" (one gets removed) if their core-setup fields
    (temperature, force field, water model, ligand state) are
    compatible. For mdrepo/mddb/dynodb, every field is a single known
    value, so "compatible" means exact equality (see core_setup_key).
    mdCATH is the one dataset whose per-record temperature can't be
    pinned to a single value -- it genuinely runs 5 replicas per domain
    (320/348/379/413/450 K), always with CHARMM22* in TIP3P, always
    apo. Rather than leaving it permanently unmatched or falling back
    to reading free-text descriptions, that known protocol is encoded
    once as a fixed set (FLEXIBLE_CONDITION_SETS) and an mdCATH record
    is matched if some other record in the same sequence group has a
    single value falling inside that set. This is the same underlying
    inference as ordinary exact matching (same sequence + a simulation
    protocol consistent with the other party's), just generalized from
    "equal to one value" to "a member of a known set of values" for the
    one dataset that structurally needs it -- no separate provenance
    pass, no reading COLLECTIONS tags or description text.

    DynoDB does NOT need this treatment: its protocol is a single fixed
    value per field (confirmed directly from the paper's Methods
    section: 300 K, CHARMM36m, TIP3P, always apo -- see
    normalize_dynodb()), so it is matched by plain equality like
    mdrepo/mddb. Its fixed 300 K also never falls inside mdCATH's
    {320,348,379,413,450} set, so a DynoDB record is never matched
    against mdCATH here -- which is the correct outcome: the paper's
    methods state every one of DynoDB's entries gets its own fresh
    100 ns GROMACS run, it never redistributes another dataset's
    trajectory, so sharing a starting structure with mdCATH is not
    evidence of a duplicated simulation. (An earlier version of this
    script removed 20 DynoDB records on exactly that mistaken
    assumption, using DynoDB's self-reported inmdCATH flag as if it
    meant "copied trajectory" -- verified empirically that all 20 have
    DynoDB's fixed 300 K, which never matches mdCATH's set, so under
    this rule they are correctly kept.)

    If a record's own fields are missing/unknown AND it isn't a
    registered "flexible" dataset, it is never auto-removed -- only
    flagged as unverifiable -- since sameness can't be confirmed from
    incomplete data.

This does NOT modify or delete anything in mdrepo/mddb/mdcath/dynodb's
own simulations/ folders. It reads them read-only and writes a new,
separate combined view here.

Usage:
    python3 dedup.py
"""

import json
import os
from collections import defaultdict

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # data/analysis
OUTPUT_FILE = "cross_database_dedup_report.json"
COMBINED_TABLE_FILE = "combined_records.json"


def load_json(path):
    with open(path) as f:
        return json.load(f)


def parse_fasta_chains(fasta_text):
    """mdrepo's fasta_sequence is itself a mini multi-record FASTA string
    (e.g. ">1\\nSEQ1\\n>2\\nSEQ2\\n" for multi-chain systems -- 1,124 of
    3,485 records, 32%, are multi-chain). Same parser as
    data/preprocess/mdrepo/fetch.py's parse_fasta_chains(), so multi-chain
    sequences flatten via "X" the same way the rest of this project
    already does, instead of naively joining lines (which would splice a
    literal '>2' header into the sequence text for every multi-chain
    record)."""
    if not fasta_text:
        return []
    chains, current = [], []
    for line in fasta_text.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current:
                chains.append("".join(current))
                current = []
        else:
            current.append(line)
    if current:
        chains.append("".join(current))
    return [c.upper() for c in chains]


def normalize_mdrepo():
    """Returns (records, n_excluded_as_heteromeric_complex)."""
    d = os.path.join(BASE, "mdrepo", "simulations")
    out = []
    n_excluded = 0
    for fn in os.listdir(d):
        rec = load_json(os.path.join(d, fn))
        detail = rec.get("detail") or {}
        chains = parse_fasta_chains(detail.get("fasta_sequence"))
        # Exclude heteromeric complexes (chains with genuinely different
        # sequences) here, at harmonization -- the earliest point after
        # collection where this can be enforced without re-touching any
        # network calls or downloaded trajectory data. Single-chain
        # records and homo-oligomers (multiple chains, all identical)
        # pass through unchanged.
        if len(set(chains)) > 1:
            n_excluded += 1
            continue
        seq = "X".join(chains) if chains else None
        if not seq:
            continue
        out.append({
            "dataset": "mdrepo",
            "record_id": rec.get("slug") or str(rec.get("simulation_id")),
            "sequence": seq,
            "uniprot_id": rec.get("uniprot_id"),
            "organism": rec.get("organism"),
            "temperature": detail.get("temperature"),
            "force_field": detail.get("forcefield"),
            "water_model": detail.get("water_type"),
            "ligand_state": "apo" if not detail.get("ligands") else "bound",
        })
    return out, n_excluded


def normalize_mddb():
    """Returns (records, n_excluded_as_heteromeric_complex)."""
    d = os.path.join(BASE, "mddb", "simulations")
    out = []
    n_excluded = 0
    for fn in os.listdir(d):
        rec = load_json(os.path.join(d, fn))
        md = (rec.get("detail") or {}).get("metadata") or {}
        seqs = md.get("PROTSEQ")  # PROTSEQ, not SEQUENCES -- see data_quality note
        # Exclude heteromeric complexes -- see normalize_mdrepo() for why
        # this is enforced here, at harmonization.
        if seqs and len(set(seqs)) > 1:
            n_excluded += 1
            continue
        seq = "X".join(seqs) if seqs else None
        if not seq:
            continue
        ff = md.get("FF")
        out.append({
            "dataset": "mddb",
            "record_id": rec.get("accession"),
            "sequence": seq,
            "uniprot_id": rec.get("uniprot_id"),
            "organism": rec.get("organism"),
            "temperature": md.get("TEMP"),
            "force_field": "; ".join(ff) if isinstance(ff, list) else ff,
            "water_model": md.get("WAT"),
            "ligand_state": "apo" if not md.get("INPUT_LIGANDS") else "bound",
        })
    return out, n_excluded


def normalize_mdcath():
    """Returns (records, n_excluded_as_heteromeric_complex). No filtering
    performed: a CATH domain is defined as a segment of a single PDB
    chain by construction, so no mdCATH record can be a heteromeric
    complex -- see collect_metadata.py's module docstring."""
    d = os.path.join(BASE, "mdcath", "simulations")
    out = []
    for fn in os.listdir(d):
        rec = load_json(os.path.join(d, fn))
        seq = rec.get("sequence")
        if not seq:
            continue
        out.append({
            "dataset": "mdcath",
            "record_id": rec.get("domain_id"),
            "sequence": seq,
            "uniprot_id": rec.get("uniprot_id"),
            "organism": rec.get("organism"),
            # Genuinely 5 temperatures/domain (320-450K), not one value --
            # left unset deliberately so this dataset never auto-collapses.
            "temperature": None,
            "force_field": None,
            "water_model": None,
            "ligand_state": None,
        })
    return out, 0


def normalize_dynodb():
    """Returns (records, n_excluded_as_heteromeric_complex)."""
    d = os.path.join(BASE, "dynodb", "simulations")
    out = []
    n_excluded = 0
    for fn in os.listdir(d):
        rec = load_json(os.path.join(d, fn))
        seq = rec.get("sequence")
        if not seq:
            continue
        # Exclude heteromeric complexes -- see normalize_mdrepo(). DynoDB's
        # own metadata already carries a precomputed distinct-chain count,
        # so no sequence parsing is needed (mirrors collect_metadata.py's
        # STEP4, applied here too since simulations/ already has records
        # collected before that check existed).
        n_unique_chains = (rec.get("metadata_row") or {}).get("Unique_Protein_Chains")
        if n_unique_chains is not None and n_unique_chains > 1:
            n_excluded += 1
            continue
        out.append({
            "dataset": "dynodb",
            "record_id": rec.get("pdb_id"),
            "sequence": seq,
            "uniprot_id": rec.get("uniprot_id"),
            "organism": rec.get("organism"),
            # Fixed protocol, confirmed directly from the paper's Methods
            # section (Hu et al., Sci Data 2026, "MD Simulation Pipeline"):
            # "The temperature was maintained at 300 K using the
            # velocity-rescaling (V-rescale) thermostat... CHARMM36m force
            # field (July 2020 release)... solvated with TIP3P water
            # molecules." Ligand state confirmed from the same paper's
            # Technical Validation section: "all ligands, metal ions,
            # cofactors, and chemical modifications were removed" during
            # setup, so every DynoDB system is apo. Unlike mdCATH (which
            # genuinely varies temperature per domain), DynoDB uses one
            # single fixed value for the entire dataset, so these are
            # real constants, not guesses.
            "temperature": 300,
            "force_field": "CHARMM36m",
            "water_model": "TIP3P",
            "ligand_state": "apo",
        })
    return out, n_excluded


def core_setup_key(rec):
    fields = (rec["temperature"], rec["force_field"], rec["water_model"], rec["ligand_state"])
    if any(f is None for f in fields):
        return None  # incomplete -- never comparable
    return fields


# Datasets whose per-record condition fields can't be pinned to a single
# value, but whose real simulation protocol is fully known as a small
# fixed set of possibilities. A record from a "plain" dataset (a real
# core_setup_key) is treated as a match for one of these if its values
# fall inside the corresponding sets. See module docstring for why
# mdCATH is the only entry here and DynoDB is not.
FLEXIBLE_CONDITION_SETS = {
    "mdcath": {
        "temperature": {320, 348, 379, 413, 450},
        "force_field": {"CHARMM22*"},
        "water_model": {"TIP3P"},
        "ligand_state": {"apo"},
    },
}


def matches_flexible_protocol(candidate, flexible_spec):
    """True if `candidate` (a record with a real core_setup_key) has
    values consistent with a flexible-dataset's known protocol."""
    return (candidate["temperature"] in flexible_spec["temperature"]
            and candidate["force_field"] in flexible_spec["force_field"]
            and candidate["water_model"] in flexible_spec["water_model"]
            and candidate["ligand_state"] in flexible_spec["ligand_state"])


def main():
    print("Loading and normalizing all four datasets...")
    mdrepo_records, n_excl_mdrepo = normalize_mdrepo()
    mddb_records, n_excl_mddb = normalize_mddb()
    mdcath_records, n_excl_mdcath = normalize_mdcath()
    dynodb_records, n_excl_dynodb = normalize_dynodb()
    records = mdrepo_records + mddb_records + mdcath_records + dynodb_records
    excluded_complex_by_dataset = {
        "mdrepo": n_excl_mdrepo, "mddb": n_excl_mddb,
        "mdcath": n_excl_mdcath, "dynodb": n_excl_dynodb,
    }
    print(f"  {len(records)} total records with a usable sequence, after excluding "
          f"{sum(excluded_complex_by_dataset.values())} heteromeric complex(es) "
          f"({excluded_complex_by_dataset}).\n")

    print("Grouping by exact sequence string...")
    by_sequence = defaultdict(list)
    for r in records:
        by_sequence[r["sequence"]].append(r)
    multi_groups = {s: recs for s, recs in by_sequence.items() if len(recs) > 1}
    print(f"  {len(by_sequence)} distinct exact sequences, "
          f"{len(multi_groups)} shared by more than one record.\n")

    print("Resolving duplicates within each shared-sequence group...")
    kept, removed = [], []
    cross_dataset_removed_pairs = defaultdict(int)
    same_sequence_kept_different_condition = 0

    for seq, recs in by_sequence.items():
        if len(recs) == 1:
            r = recs[0]
            r["status"] = "kept"
            r["reason"] = "unique sequence"
            kept.append(r)
            continue

        # Step 1: within this exact-sequence group, sub-group records that
        # have a real core_setup_key by exact match. Records without one
        # (unresolved -- currently just mdcath, whose fields are always
        # None) are set aside for step 2.
        by_setup = defaultdict(list)
        unresolved = []
        for r in recs:
            key = core_setup_key(r)
            if key is None:
                unresolved.append(r)
            else:
                by_setup[key].append(r)

        group_representatives = []  # scalar winners, used as step 2 match targets
        for group in by_setup.values():
            # True duplicates: same sequence AND same core setup. Keep the
            # first (deterministic by dict insertion order), remove the rest.
            winner, *dupes = group
            winner["status"] = "kept"
            winner["reason"] = (f"representative of {len(group)} identical (sequence+setup) records"
                                 if len(group) > 1 else
                                 "same sequence as other record(s), but core setup differs -- kept")
            kept.append(winner)
            group_representatives.append(winner)
            if len(group) == 1:
                same_sequence_kept_different_condition += 1
            for d in dupes:
                d["status"] = "removed_duplicate"
                d["reason"] = f"100% identical sequence and core setup to {winner['dataset']}:{winner['record_id']}"
                d["duplicate_of"] = f"{winner['dataset']}:{winner['record_id']}"
                removed.append(d)
                pair = tuple(sorted([winner["dataset"], d["dataset"]]))
                cross_dataset_removed_pairs[f"{pair[0]}<->{pair[1]}"] += 1

        # Step 2: try to resolve each unresolved record against this same
        # group's scalar representatives, using its dataset's known fixed
        # protocol (FLEXIBLE_CONDITION_SETS) if it has one registered.
        for r in unresolved:
            flexible_spec = FLEXIBLE_CONDITION_SETS.get(r["dataset"])
            target = None
            if flexible_spec:
                target = next((o for o in group_representatives
                                if matches_flexible_protocol(o, flexible_spec)), None)
            if target:
                r["status"] = "removed_duplicate"
                r["reason"] = (f"same sequence and a simulation condition consistent with "
                                f"{r['dataset']}'s known fixed protocol, matching "
                                f"{target['dataset']}:{target['record_id']}")
                r["duplicate_of"] = f"{target['dataset']}:{target['record_id']}"
                removed.append(r)
                pair = tuple(sorted([target["dataset"], r["dataset"]]))
                cross_dataset_removed_pairs[f"{pair[0]}<->{pair[1]}"] += 1
            else:
                r["status"] = "kept"
                r["reason"] = "same sequence as other record(s), but core setup " \
                               "incomplete/unverifiable -- kept per safe rule"
                kept.append(r)
                same_sequence_kept_different_condition += 1

    print(f"\nDedup done: {len(kept)} kept, {len(removed)} removed as duplicate simulations.\n")

    combined = kept + removed
    with open(COMBINED_TABLE_FILE, "w") as f:
        json.dump(combined, f, indent=2)

    report = {
        "excluded_heteromeric_complex_by_dataset": excluded_complex_by_dataset,
        "total_input_records": len(records),
        "distinct_exact_sequences": len(by_sequence),
        "sequences_shared_by_multiple_records": len(multi_groups),
        "records_kept": len(kept),
        "records_removed_as_true_duplicates": len(removed),
        "same_sequence_kept_due_to_different_or_unverifiable_condition": same_sequence_kept_different_condition,
        "removed_duplicates_by_dataset_pair": dict(cross_dataset_removed_pairs),
        "kept_by_dataset": {ds: sum(1 for r in kept if r["dataset"] == ds)
                             for ds in ["mdrepo", "mddb", "mdcath", "dynodb"]},
        "removed_by_dataset": {ds: sum(1 for r in removed if r["dataset"] == ds)
                                for ds in ["mdrepo", "mddb", "mdcath", "dynodb"]},
    }
    with open(OUTPUT_FILE, "w") as f:
        json.dump(report, f, indent=2)

    print(json.dumps(report, indent=2))
    print(f"\nWrote {OUTPUT_FILE} and {COMBINED_TABLE_FILE}.")


if __name__ == "__main__":
    main()
