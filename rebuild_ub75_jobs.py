#!/usr/bin/env python3
"""Rebuild jobs_ub75/ with valence-correct CCDs. The originals carry BOTH Julian bugs.

WHY. jobs_ub75/ was written AFTER the valence fix but reverted to the pre-fix
components, so the Ub(1-75)-vs-Ub(1-76) UBE2W comparison was run on chemistry known to
be wrong. Both bugs from results/JULIAN_RUN_ANALYSIS.md are present, verified with
ccd_valence.py under the exact external bonds each job declares:

  BUG 1  LIG-1 (the XisoK isopeptide moiety) still carries the spurious C01-H01 bond:
             jobs_ub75 LIG-1 : O01=C01, C01-C02, C01-H01   -> C01 at 5/4 OVER-VALENT
             ccd_v2/LisoK    : O01=C01, C01-C02            -> C01 at 4/4 with the bond
         The acyl carbon reaches five coordination once the isopeptide bond to the
         substrate lysine NZ is added.

  BUG 2  the terminal glycines are bare `GLY`, which retains OXT:
             N  at 4/3 OVER-VALENT   (it also takes the peptide bond from Ub Arg74)
             C  at 5/4 OVER-VALENT   (C, O, OXT, CA, plus the external bond)
         Note BOTH atoms are over-valent, not just the acyl carbon. The consequence is
         the 328-deg sp3 thioester instead of the 357-360 deg planar sp2 measured in
         STUB_RESULTS.md once the no-OXT CCDs are used.

WHAT CHANGES. Same science, correct chemistry:
  * LIG-1              -> ccd_v2/LisoK_userCCD.cif      (no C01-H01)
  * GLY                -> UBG1  for Ub(1-75)            (Gly75 alone, no OXT)
  * GLY + GLY          -> UBGG  for Ub(1-76)            (Gly75-Gly76 dipeptide, no OXT)

Collapsing two GLY ligands into one UBGG dipeptide also removes a bond: the internal
Gly75-Gly76 peptide is now inside the CCD rather than declared, which is one fewer thing
for AF3 to discard.

The `free` jobs (no ubiquitin, no ligand) are unaffected and are copied unchanged.

Usage: python rebuild_ub75_jobs.py OUTDIR
"""
import argparse
import json
import os
import sys

from build_uba1_jobs import check_job, embed_ccd

HERE = os.path.dirname(os.path.abspath(__file__))
CCD_V2 = os.path.join(HERE, "ccd_v2")
UBE2W_CAT_CYS = 91
SUB_LYS = 12                # SUMO2 construct residue 12 = native K11
UB_SPLIT = 74

TAIL = {76: dict(ccd="UBGG", acyl="C2"), 75: dict(ccd="UBG1", acyl="C1")}


def ccd_code(stem):
    """The component id declared INSIDE the CIF, which is what AF3 resolves against.

    Read it rather than assume it matches the filename: ccd_v2/LisoK_userCCD.cif declares
    `_chem_comp.id LIG-1`, so a job asking for ccdCodes ['LisoK'] would reference a
    component that does not exist in its own userCCD.
    """
    path = os.path.join(CCD_V2, f"{stem}_userCCD.cif")
    for line in open(path):
        if line.startswith("_chem_comp.id"):
            return line.split()[1]
    raise ValueError(f"no _chem_comp.id in {path}")


LISOK_CODE = ccd_code("LisoK")


def rebuild(old, tail_len):
    """Rewrite one charged job with valence-correct ligands."""
    d = json.loads(json.dumps(old))
    t = TAIL[tail_len]

    # keep the protein chains exactly as they were
    seqs = [s for s in d["sequences"] if "protein" in s]
    # The corrected LisoK: no spurious C01-H01. NOTE the ccdCode must be the id INSIDE
    # the CIF (`_chem_comp.id LIG-1`), not the filename stem (LisoK_userCCD.cif) -- AF3
    # resolves the ligand by the data_ block name, and a mismatch means it cannot find
    # the component at all.
    seqs.append({"ligand": {"id": "L", "ccdCodes": [LISOK_CODE]}})
    # one tail ligand, no OXT, replacing the one or two bare GLYs
    seqs.append({"ligand": {"id": "G", "ccdCodes": [t["ccd"]]}})
    d["sequences"] = seqs

    d["bondedAtomPairs"] = [
        [["A", SUB_LYS, "NZ"], ["L", 1, "C01"]],          # isopeptide to substrate Lys
        [["U", UB_SPLIT, "C"], ["G", 1, "N1"]],           # extend Ub backbone into tail
        [["G", 1, t["acyl"]], ["B", UBE2W_CAT_CYS, "SG"]],  # the thioester
    ]
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    bad = 0
    for tail in (76, 75):
        src = os.path.join(HERE, "jobs_ub75", f"ube2w_ub{tail}_charged.json")
        old = json.load(open(src))
        d = rebuild(old, tail)
        d["name"] = f"ube2w_ub{tail}_charged_v2"

        problems = check_job(d, CCD_V2, d["name"])
        if problems:
            bad += 1
            print(f"FAIL {d['name']}")
            for p in problems:
                print(f"       {p}")
            continue
        embed_ccd(d, CCD_V2)
        out = os.path.join(a.outdir, f"{d['name']}.json")
        json.dump(d, open(out, "w"), indent=2)
        ligs = [s["ligand"]["ccdCodes"][0] for s in d["sequences"] if "ligand" in s]
        print(f"  {d['name']:28s} ligands={ligs} bonds={len(d['bondedAtomPairs'])}")

    # the free jobs have no ligands and no bonds: nothing to fix
    for tail in (76, 75):
        src = os.path.join(HERE, "jobs_ub75", f"ube2w_ub{tail}_free.json")
        if os.path.exists(src):
            d = json.load(open(src))
            json.dump(d, open(os.path.join(a.outdir, os.path.basename(src)), "w"),
                      indent=2)
            print(f"  {os.path.basename(src):28s} copied unchanged (no ligands, no bonds)")

    if bad:
        sys.exit(f"\n{bad} job(s) failed the gate; nothing written for those")
    print(f"\nwritten to {a.outdir}/")


if __name__ == "__main__":
    main()
