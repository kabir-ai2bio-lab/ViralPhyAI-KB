#!/usr/bin/env python3
"""
final_selection.py

Selects one representative record per distinct sequence from the Merged
Dataset (combined_records.json, status == "kept") produced by dedup.py --
the final representative selection cascade described in the paper's
Methods section, "Final representative selection": prefer apo -> prefer
TIP3P -> closest to 310 K among records that have a temperature -> first
encountered. Every step is a soft filter (falls back to the unfiltered
pool if it would eliminate every candidate), so every distinct sequence
always gets a representative, including ones with no condition data on
any candidate record.

Running this script directly writes final_representative_ids.json into
../../features_extraction/, grouping each winner's record_id by source
dataset. That file is the worklist convert_mdcath.py, convert_dynodb.py,
and run_feature_extraction.py all read to know which records to process.

Usage:
    python3 final_selection.py
"""

import json
import os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
COMBINED_RECORDS = os.path.join(HERE, "combined_records.json")
FINAL_IDS_OUT = os.path.join(HERE, "..", "..", "features_extraction", "final_representative_ids.json")

TARGET_TEMP_K = 310  # human physiological temperature -- convention used throughout this project


def load_kept():
    with open(COMBINED_RECORDS) as f:
        combined = json.load(f)
    return [r for r in combined if r["status"] == "kept"]


def soft_filter(candidates, predicate):
    """Keep only candidates matching predicate, unless that would
    eliminate every candidate, in which case the original list passes
    through unchanged."""
    kept = [c for c in candidates if predicate(c)]
    return kept if kept else candidates


def select_final_representative(records):
    """One representative per unique-sequence group: prefer apo ->
    prefer TIP3P -> closest to 310 K (only among records that have a
    temperature at all) -> first-encountered."""
    pool = records
    pool = soft_filter(pool, lambda r: r.get("ligand_state") == "apo")
    pool = soft_filter(pool, lambda r: r.get("water_model") == "TIP3P")

    with_temp = [r for r in pool if r.get("temperature") is not None]
    if with_temp:
        min_dist = min(abs(r["temperature"] - TARGET_TEMP_K) for r in with_temp)
        pool = [r for r in with_temp if abs(r["temperature"] - TARGET_TEMP_K) == min_dist]

    return pool[0]


def compute_final_selection(kept):
    """One winner per unique-sequence group. Returns a list of
    (winner_record, num_candidates) tuples."""
    by_seq = defaultdict(list)
    for r in kept:
        by_seq[r["sequence"]].append(r)
    return [(select_final_representative(recs), len(recs)) for recs in by_seq.values()]


def main():
    kept = load_kept()
    final_selection = compute_final_selection(kept)

    by_source = defaultdict(list)
    for winner, _n_candidates in final_selection:
        by_source[winner["dataset"]].append(winner["record_id"])

    os.makedirs(os.path.dirname(FINAL_IDS_OUT), exist_ok=True)
    with open(FINAL_IDS_OUT, "w") as f:
        json.dump(dict(by_source), f, indent=2)

    print(f"{len(final_selection)} final representatives selected (from {len(kept)} kept records).")
    for source, ids in sorted(by_source.items()):
        print(f"  {source}: {len(ids)}")
    print(f"Wrote {FINAL_IDS_OUT}")


if __name__ == "__main__":
    main()
