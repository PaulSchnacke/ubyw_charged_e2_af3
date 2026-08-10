#!/usr/bin/env python3
"""Graft precomputed MSAs from an augmented AF3 JSON onto a fresh job JSON.

All three sequence sets we need were already searched in earlier rounds, so no
data-pipeline run is required for the sweep. MSAs depend ONLY on protein
sequence, so an alignment computed for one bond spec / ligand set is valid for
any other job with the same protein chains.

This copies the MSA fields (unpairedMsa, pairedMsa, templates) from the donor
onto the recipient, matched BY SEQUENCE rather than by chain id -- chain letters
differ between the old and new job layouts, and matching by id would silently
attach the wrong alignment. Every recipient protein chain must find an exact
sequence match or the script refuses.

Usage: python graft_msa.py NEW_JOB.json DONOR_data.json OUT.json
"""
import json
import sys


def main(new_job, donor, out):
    d = json.load(open(new_job))
    src = json.load(open(donor))

    # donor alignments keyed by exact sequence
    pool = {}
    for s in src.get("sequences", []):
        if "protein" not in s:
            continue
        p = s["protein"]
        if p.get("unpairedMsa"):
            pool[p["sequence"]] = {k: p.get(k) for k in
                                   ("unpairedMsa", "pairedMsa", "templates")}
    if not pool:
        sys.exit(f"donor {donor} carries no MSAs -- is it an augmented _data.json?")

    grafted, missing = [], []
    for s in d["sequences"]:
        if "protein" not in s:
            continue
        p = s["protein"]
        hit = pool.get(p["sequence"])
        if hit is None:
            missing.append((p["id"], len(p["sequence"])))
            continue
        for k, v in hit.items():
            if v is not None:
                p[k] = v
        grafted.append((p["id"], len(p["sequence"])))

    if missing:
        sys.exit(f"no MSA in donor for chain(s) {missing}; donor has lengths "
                 f"{sorted(len(x) for x in pool)} -- refusing to write a job with "
                 f"partial alignments")

    # the things that must survive the graft
    assert d.get("userCCD") or not any("ligand" in s for s in d["sequences"]), \
        "ligand present but userCCD lost"
    nb = len(d.get("bondedAtomPairs", []))

    json.dump(d, open(out, "w"), indent=1)
    print(f"grafted MSAs onto {len(grafted)} protein chains {grafted}")
    print(f"  bondedAtomPairs preserved: {nb}")
    print(f"  userCCD preserved: {'userCCD' in d}")
    print(f"  seeds: {len(d.get('modelSeeds', []))}  -> {out}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
