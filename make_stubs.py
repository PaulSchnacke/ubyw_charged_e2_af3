#!/usr/bin/env python3
"""Turn a production job into a STUB: single-sequence MSAs, one seed.

WHY. A stub runs inference in ~10 min of GPU instead of 4-6 h of database search,
and it answers the only questions that must be settled before real compute:

  * did AF3 KEEP every declared bond, or silently discard some?
    (grep the log for 'Reducing number of bonds' -- exit 0 either way)
  * what bond lengths and angles came out, i.e. is the chemistry right?
  * how did AF3 renumber the chains? It absorbs a leading single-residue ligand
    into the protein chain, so the covale records, not the input naming, define
    where the atoms ended up.

A stub says NOTHING about placement -- single-sequence AF3 has no coevolution
signal, so poses are unreliable. It is a chemistry and topology check only.

Setting unpairedMsa to just the query and pairedMsa to '' is what makes AF3 skip
the data pipeline; templates must be an explicit empty list.

Usage: python make_stubs.py JOB.json [JOB.json ...] OUTDIR
"""
import json
import os
import sys


def stubify(job):
    d = json.loads(json.dumps(job))          # deep copy
    d["name"] = "stub_" + d["name"]
    d["modelSeeds"] = [1]
    for s in d["sequences"]:
        if "protein" not in s:
            continue
        p = s["protein"]
        p["unpairedMsa"] = f">q_{p['id']}\n{p['sequence']}\n"
        p["pairedMsa"] = ""
        p["templates"] = []
    return d


def main():
    *jobs, outdir = sys.argv[1:]
    os.makedirs(outdir, exist_ok=True)
    for f in jobs:
        d = stubify(json.load(open(f)))
        out = os.path.join(outdir, f"{d['name']}.json")
        with open(out, "w") as fh:
            json.dump(d, fh, indent=2)
        nb = len(d.get("bondedAtomPairs", []))
        print(f"  {d['name']:40s} {nb} bonds -> {out}")


if __name__ == "__main__":
    main()
