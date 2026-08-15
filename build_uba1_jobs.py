#!/usr/bin/env python3
"""Build the UBA1 AF3 jobs: two sites x two tail lengths x UBE2W present/absent.

THE QUESTION. Paul engineered the system to accept Ub(1-75) as well as Ub(1-76);
experimentally it failed. The failure could sit at UBA1 (charging) or at UBE2W
(accepting the shorter donor). The UBE2W half is already modelled (jobs_ub75/).
This is the UBA1 half.

UBA1 has TWO chemically distinct catalytic centres, and they are modelled
differently on purpose:

  ADENYLATION SITE (~K528/R478 region, ATP-binding)
      Ub C-terminal carboxyl + ATP.Mg -> Ub-AMP, a MIXED ANHYDRIDE (acyl
      phosphate). NOT modelled covalently here: no published AF3-compatible CCD
      for a Ub-AMP acyl phosphate exists in this project, and an acyl phosphate is
      exactly the chemistry AF3 handles worst (it has no thioester chemistry
      either -- the acyl carbon pyramidalises to 338 deg). Modelled instead as
      Ub + ATP + Mg with NO bond declared, asking the reach question:
      how close does the Ub C-terminal carbonyl carbon get to the ATP alpha
      phosphate, and does the shorter tail change that?

  CATALYTIC CYSTEINE (Cys632, UniProt: "Glycyl thioester intermediate")
      Ub Gly76 (or Gly75) thioester on Cys632 SG. Modelled BOTH ways:
      non-covalently (reach) and covalently via the established ligand trick,
      because this project's thioester machinery is validated.

DESIGN. 2 sites x {Ub 1-76, Ub 1-75} x {UBE2W absent, UBE2W present} = 8 jobs,
plus 4 covalent-thioester jobs = 12. Every covalent job uses a custom
no-OXT glycine CCD; bare GLY is never used on a reactive C-terminus (see
make_uba1_ccd.py and results/JULIAN_RUN_ANALYSIS.md Bug 2).

THE POLYMER-POLYMER RULE (docs/AF3_POLYMER_BOND_LIMIT.md). AF3 3.0.1 silently
discards covalent bonds between two POLYMER chains -- exit 0, atoms 26 A apart.
So a covalent thioester must express one partner as a LIGAND. Ubiquitin therefore
enters as protein chain U = Ub 1-74, with the terminal glycine(s) as ligand
chains:
    Ub 1-76 thioester:  U(1-74) + UBGG ligand (Gly75-Gly76)  -> C2 bonds SG
    Ub 1-75 thioester:  U(1-74) + UBG1 ligand (Gly75)        -> C1 bonds SG

Usage: python build_uba1_jobs.py OUTDIR [--seeds N] [--ccd DIR]
"""
import argparse
import json
import os
import sys

from ccd_valence import check_ccd

HERE = os.path.dirname(os.path.abspath(__file__))
SEQS = json.load(open(os.path.join(HERE, "seqs.json")))
UBE2W = SEQS["UBE2W"]
UBIQUITIN = SEQS["UBIQUITIN"]          # 1-76, verified == UniProt P0CG47[:76]

# UBA1, human, UniProt P22314 (1058 aa). Catalytic Cys632 is annotated
# "Glycyl thioester intermediate"; ATP-binding residues 478/504/515/528/576-577.
UBA1 = json.load(open(os.path.join(HERE, "seqs_uba1.json")))["UBA1"]
UBA1_CAT_CYS = 632
UBA1_ATP_SITE = [478, 504, 515, 528, 576, 577]

# Ligand atom roles. UBGG is the Gly75-Gly76 dipeptide (thioester carbon C2);
# UBG1 is Gly75 alone (thioester carbon C1). Both have NO OXT.
TAIL = {
    76: dict(ccd="UBGG", n_atom="N1", acyl="C2", nres=2),
    75: dict(ccd="UBG1", n_atom="N1", acyl="C1", nres=1),
}
UB_SPLIT = 74          # Ub 1-74 stays a protein chain, so it keeps its own MSA


def base(name, seeds):
    return {"name": name, "sequences": [], "modelSeeds": list(seeds),
            "dialect": "alphafold3", "version": 1, "bondedAtomPairs": []}


def add_prot(d, cid, seq):
    d["sequences"].append({"protein": {"id": cid, "sequence": seq}})


def add_lig(d, cid, code):
    d["sequences"].append({"ligand": {"id": cid, "ccdCodes": [code]}})


def job_noncovalent(tail_len, with_ube2w, site, seeds, name):
    """Ub + UBA1 with NO bond declared -- the honest reach question.

    site='adenylation' adds ATP and Mg; site='cys' leaves the catalytic cysteine
    free. Ubiquitin is one intact protein chain here (no ligand split needed,
    because there is no covalent bond to express).
    """
    d = base(name, seeds)
    add_prot(d, "A", UBA1)
    add_prot(d, "U", UBIQUITIN[:tail_len])
    if with_ube2w:
        add_prot(d, "B", UBE2W)
    if site == "adenylation":
        add_lig(d, "T", "ATP")
        add_lig(d, "M", "MG")
    d.pop("bondedAtomPairs")            # AF3 rejects an empty bond list
    return d


def job_thioester(tail_len, with_ube2w, seeds, name):
    """Ub thioester on UBA1 Cys632, via the ligand trick.

    Ub 1-74 is a protein chain; the terminal glycine(s) are a ligand, so the
    inter-chain bond has a ligand partner and AF3 will not discard it. Bonds:
        U:74:C      -> tail:N1        extend the backbone
        tail:acyl   -> A:632:SG       the thioester
    """
    t = TAIL[tail_len]
    d = base(name, seeds)
    add_prot(d, "A", UBA1)
    add_prot(d, "U", UBIQUITIN[:UB_SPLIT])
    if with_ube2w:
        add_prot(d, "B", UBE2W)
    add_lig(d, "G", t["ccd"])
    d["bondedAtomPairs"] = [
        [["U", UB_SPLIT, "C"], ["G", 1, t["n_atom"]]],
        [["G", 1, t["acyl"]], ["A", UBA1_CAT_CYS, "SG"]],
    ]
    return d


def external_bonds(d):
    """Count external bonds AF3 will create per ligand atom, for the valence gate."""
    lig_ids = {s["ligand"]["id"] for s in d["sequences"] if "ligand" in s}
    ext = {}
    for (c1, _r1, a1), (c2, _r2, a2) in d.get("bondedAtomPairs", []):
        for cid, atom in ((c1, a1), (c2, a2)):
            if cid in lig_ids:
                ext[atom] = ext.get(atom, 0) + 1
    return ext


def check_job(d, ccddir, label):
    """Gate every job before it is written. Three checks, each for a real bug.

    1. VALENCE of every custom ligand, with the external bonds AF3 will add
       included -- the bug class that produced a pentavalent carbon rendered as a
       cyclopropene (results/JULIAN_RUN_ANALYSIS.md).
    2. Every bonded atom names a chain that EXISTS in the job.
    3. Every bond has at least one LIGAND partner, or AF3 3.0.1 silently discards
       it (docs/AF3_POLYMER_BOND_LIMIT.md).
    """
    problems = []
    prot_ids = {s["protein"]["id"] for s in d["sequences"] if "protein" in s}
    lig_ids = {s["ligand"]["id"] for s in d["sequences"] if "ligand" in s}
    ext = external_bonds(d)

    for s in d["sequences"]:
        if "ligand" not in s:
            continue
        code = s["ligand"]["ccdCodes"][0]
        if code in ("ATP", "MG"):
            continue                    # standard CCDs, no external bonds declared
        path = os.path.join(ccddir, f"{code}_userCCD.cif")
        if not os.path.exists(path):
            problems.append(f"missing CCD {path}")
            continue
        check_ccd(open(path).read(), external=ext, label=f"{label}:{code}")

    for pair in d.get("bondedAtomPairs", []):
        chains = [p[0] for p in pair]
        for cid in chains:
            if cid not in prot_ids | lig_ids:
                problems.append(f"bond names chain {cid}, not in job")
        if not any(c in lig_ids for c in chains):
            problems.append(f"POLYMER-POLYMER bond {pair} -- AF3 will discard it")

    # the catalytic cysteine must really be a cysteine at that position
    for pair in d.get("bondedAtomPairs", []):
        for cid, resi, atom in pair:
            if atom == "SG" and cid in prot_ids:
                seq = next(s["protein"]["sequence"] for s in d["sequences"]
                           if "protein" in s and s["protein"]["id"] == cid)
                if seq[resi - 1] != "C":
                    problems.append(f"residue {cid}:{resi} is {seq[resi-1]}, not CYS")
    return problems


def embed_ccd(d, ccddir):
    """Concatenate the custom CCD blocks this job needs into userCCD."""
    blocks = []
    for s in d["sequences"]:
        if "ligand" not in s:
            continue
        code = s["ligand"]["ccdCodes"][0]
        if code in ("ATP", "MG"):
            continue
        blocks.append(open(os.path.join(ccddir,
                                        f"{code}_userCCD.cif")).read().rstrip() + "\n")
    if blocks:
        d["userCCD"] = "\n".join(blocks)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--seeds", type=int, default=25)
    ap.add_argument("--ccd", default=os.path.join(HERE, "ccd_uba1"))
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    seeds = range(1, a.seeds + 1)

    jobs = {}
    for tail in (76, 75):
        for e2 in (False, True):
            e2tag = "with_ube2w" if e2 else "no_ube2w"
            # non-covalent reach at both sites
            for site, stag in (("adenylation", "aden"), ("cys", "cys")):
                n = f"uba1_{stag}_ub{tail}_{e2tag}"
                jobs[n] = job_noncovalent(tail, e2, site, seeds, n)
            # covalent thioester at Cys632
            n = f"uba1_thio_ub{tail}_{e2tag}"
            jobs[n] = job_thioester(tail, e2, seeds, n)

    bad = 0
    for name, d in sorted(jobs.items()):
        problems = check_job(d, a.ccd, name)
        if problems:
            bad += 1
            print(f"FAIL {name}")
            for p in problems:
                print(f"       {p}")
            continue
        embed_ccd(d, a.ccd)
        with open(os.path.join(a.outdir, f"{name}.json"), "w") as fh:
            json.dump(d, fh, indent=2)
        nprot = sum(1 for s in d["sequences"] if "protein" in s)
        nres = sum(len(s["protein"]["sequence"]) for s in d["sequences"]
                   if "protein" in s)
        nlig = sum(1 for s in d["sequences"] if "ligand" in s)
        nb = len(d.get("bondedAtomPairs", []))
        print(f"  {name:34s} {nprot} chains {nres:>5} res  {nlig} lig  {nb} bonds")

    if bad:
        sys.exit(f"\n{bad} job(s) failed the gate; nothing written for those")
    print(f"\n{len(jobs)} jobs written to {a.outdir}/ ({a.seeds} seeds each)")


if __name__ == "__main__":
    main()
