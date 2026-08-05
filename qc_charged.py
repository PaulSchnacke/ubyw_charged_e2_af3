#!/usr/bin/env python3
"""QC for the charged-E2 co-folds: does the LisoK amine reach the THIOESTER?

Usage: python qc_charged.py MODEL_DIR OUT.csv

WHAT CHANGES vs the uncharged rounds. With ubiquitin loaded onto Cys91, the atom
the substrate amine must attack is no longer the cysteine sulfur -- it is
ubiquitin's C-terminal Gly76 carboxyl carbon, the thioester carbonyl:

    Ub Gly76 C(=O)-S-Cys91          nucleophile attacks THIS carbon

So the primary measurement is  L:1:N01 -> U:76:C.

Cys91 SG is still reported, for one reason only: it makes these runs directly
comparable to rounds 2-5, where the bare SG was the target. If the amine sits at
7 A from SG in both the charged and uncharged states, that is informative; if
loading Ub pushes it away, that is the answer.

Reported per model:
    n01_ub_c      amine to the thioester carbonyl carbon   <- THE measurement
    n01_cys_sg    amine to Cys91 SG (comparability with earlier rounds)
    thioester_len Gly76 C to Cys91 SG, i.e. did AF3 honour the bond at all
    isopeptide_ok the XisoK acyl bond to the substrate lysine formed
    ub_iptm/ptm   confidence, including whether Ub is placed confidently
"""
import csv
import glob
import json
import os
import sys

import gemmi
import numpy as np

SUB, ENZ, UB, LIG = "A", "B", "U", "L"
CAT_CYS = 91
UB_CTERM = 76
NUC = "N01"
LINK = "C01"
ATTACK = 4.0            # nucleophilic attack distance on a carbonyl carbon
BOND_MAX = 2.2          # a bondedAtomPairs bond AF3 actually formed


def atoms(cif):
    st = gemmi.read_structure(cif)
    st.setup_entities()
    out = {}
    for ch in st[0]:
        for res in ch:
            for at in res:
                out[(ch.name, res.seqid.num, res.name, at.name)] = \
                    np.array(at.pos.tolist())
    return out


def analyse(cif, conf):
    A = atoms(cif)
    row = {}

    def find(chain, resnum, atom, resname=None):
        for k, v in A.items():
            if k[0] == chain and k[1] == resnum and k[3] == atom:
                if resname is None or k[2] == resname:
                    return v, k[2]
        return None, None

    nuc, _ = find(LIG, 1, NUC)
    if nuc is None:
        # the ligand may be numbered differently; fall back to searching by name
        cand = [(k, v) for k, v in A.items() if k[3] == NUC]
        if not cand:
            return {"error": f"no {NUC} atom found"}
        nuc = cand[0][1]

    ubc, ubres = find(UB, UB_CTERM, "C")
    sg, cysres = find(ENZ, CAT_CYS, "SG")

    if ubc is not None:
        row["n01_ub_c"] = float(np.linalg.norm(nuc - ubc))
        row["ub_cterm_res"] = ubres
    if sg is not None:
        row["n01_cys_sg"] = float(np.linalg.norm(nuc - sg))
        row["cys_res"] = cysres
    if ubc is not None and sg is not None:
        row["thioester_len"] = float(np.linalg.norm(ubc - sg))
        row["thioester_ok"] = row["thioester_len"] <= BOND_MAX

    # the XisoK acyl bond onto the substrate lysine
    link, _ = find(LIG, 1, LINK)
    if link is None:
        cand = [v for k, v in A.items() if k[3] == LINK]
        link = cand[0] if cand else None
    if link is not None:
        nz = [(np.linalg.norm(link - v), k) for k, v in A.items()
              if k[3] == "NZ" and k[0] == SUB and k[2] == "LYS"]
        if nz:
            d, k = min(nz)
            row["isopeptide_len"] = float(d)
            row["isopeptide_res"] = f"{k[0]}/{k[2]}{k[1]}"
            row["isopeptide_ok"] = d <= BOND_MAX

    # SPECIFICITY: is the amine near the thioester carbon rather than merely
    # near something? report the nearest carbonyl carbon of any chain.
    if ubc is not None:
        others = [(np.linalg.norm(nuc - v), f"{k[0]}/{k[2]}{k[1]}/{k[3]}")
                  for k, v in A.items()
                  if k[3] == "C" and not (k[0] == UB and k[1] == UB_CTERM)]
        if others:
            d, w = min(others)
            row["nearest_other_C"] = float(d)
            row["nearest_other_C_id"] = w

    if conf and os.path.exists(conf):
        c = json.load(open(conf))
        for key in ("iptm", "ptm", "ranking_score"):
            if key in c:
                row[key] = c[key]
        # per-chain pair iptm tells you whether Ub itself is confidently placed
        cp = c.get("chain_pair_iptm")
        if cp:
            row["n_chains_scored"] = len(cp)
    return row


def main(model_dir, out_csv):
    cifs = sorted(glob.glob(os.path.join(model_dir, "**", "*model*.cif"),
                            recursive=True))
    if not cifs:
        sys.exit(f"no *model*.cif under {model_dir}")

    rows = []
    for cif in cifs:
        base = os.path.dirname(cif)
        conf = None
        for cand in glob.glob(os.path.join(base, "*summary_confidences*.json")):
            conf = cand
            break
        r = analyse(cif, conf)
        fn = os.path.basename(cif)
        r["model"] = os.path.relpath(cif, model_dir)
        r["job"] = fn.split("__")[0] if "__" in fn else \
            os.path.basename(os.path.dirname(os.path.dirname(cif)))
        r["sample"] = fn.split("__")[1] if "__" in fn else \
            os.path.basename(base)
        rows.append(r)

    keys, seen = [], set()
    for r in rows:
        for k in r:
            if k not in seen:
                seen.add(k)
                keys.append(k)
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)

    print(f"{len(rows)} models -> {out_csv}\n")
    byjob = {}
    for r in rows:
        byjob.setdefault(r["job"], []).append(r)
    for job, rs in sorted(byjob.items()):
        d = [r["n01_ub_c"] for r in rs if r.get("n01_ub_c") is not None]
        s = [r["n01_cys_sg"] for r in rs if r.get("n01_cys_sg") is not None]
        te = [r["thioester_len"] for r in rs if r.get("thioester_len") is not None]
        ok = sum(1 for r in rs if r.get("isopeptide_ok"))
        line = f"{job:38s} n={len(rs):3d}"
        if d:
            line += (f"  N01->UbG76C min {min(d):5.2f} med "
                     f"{float(np.median(d)):5.2f}  <=4A {sum(1 for x in d if x <= ATTACK)}/{len(d)}")
        if s:
            line += f"  [SG med {float(np.median(s)):5.2f}]"
        if te:
            line += f"  thioester {float(np.median(te)):4.2f}"
        line += f"  isopep {ok}/{len(rs)}"
        print(line)

    # Distinguish "no ubiquitin chain present" from "Ub present but the bond did
    # not form". Those mean different things: the first is the wrong input file,
    # the second is AF3 declining to honour bondedAtomPairs.
    no_ub = [r for r in rows if r.get("n01_ub_c") is None]
    bad_te = [r for r in rows if r.get("thioester_len") is not None
              and r["thioester_len"] > BOND_MAX]
    if no_ub:
        print(f"\nNOTE: {len(no_ub)}/{len(rows)} models have no chain "
              f"{UB}/Gly{UB_CTERM} -- these are UNCHARGED co-folds, not charged "
              f"ones. Only the Cys91 SG column is meaningful for them (and it is "
              f"directly comparable to rounds 2-5).")
    if bad_te:
        print(f"\nWARNING: {len(bad_te)}/{len(rows)} models have Gly76 C to "
              f"Cys91 SG > {BOND_MAX} A despite ubiquitin being present -- AF3 "
              f"did not honour the thioester bond there, so those distances are "
              f"not measuring a loaded E2. Check the bondedAtomPairs entry.")
    print(f"\nNote: AF3 has no thioester chemistry. The bond length is enforced "
          f"but the geometry around it is not parameterised as a thioester, so "
          f"this is a geometric proxy for a charged E2, not a faithful one.")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
