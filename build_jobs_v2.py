#!/usr/bin/env python3
"""Build the AF3 validation jobs with valence-correct ligands.

Four species at SUMO2 K11 (native numbering; construct residue 12), plus every
SUMO2 lysine for the site sweep:

  1. xisok        XisoK alone on UBE2W        -- reproduces rounds 2-5, tests the
                                                 corrected isopeptide geometry
  2. lyscontrol   plain Lys, no modification  -- the control Paul asked for; the
                                                 earlier lysine-only run showed
                                                 the modification drives docking
  3. charged      Ub thioester on Cys91,      -- the intended pre-transfer state
                  substrate amine free
  4. tetrahedral  amine bonded to the         -- what AF3 built by ACCIDENT last
                  thioester carbon               time; now declared on purpose
  5. product      Ub isopeptide on the        -- Ub G76 joined to the XisoK amine
                  XisoK amine

Every job's bond list is checked against the ligand valences before writing.

Usage: python build_jobs_v2.py OUTDIR --ccd CCD_DIR [--sweep] [--seeds N ...]
"""
import argparse
import json
import os
import sys

from ccd_valence import valences

HERE = os.path.dirname(os.path.abspath(__file__))
SEQS = json.load(open(os.path.join(HERE, "seqs.json")))
UBE2W = SEQS["UBE2W"]
SUMO2_NATIVE = SEQS["SUMO2_NATIVE"]
UBIQUITIN = SEQS["UBIQUITIN"]

# The construct carries an N-terminal Pro that the paper introduced to block
# modification of SUMO2's native N-terminus, so every native position shifts by
# one: native K11 -> construct residue 12.
SUMO2_CONSTRUCT = "P" + SUMO2_NATIVE
OFFSET = 1
UBE2W_CAT_CYS = 91

# ligand atom roles
LIG_LINK = "C01"        # acyl carbon -> substrate Lys NZ
LIG_NUC = "N01"         # free alpha-amine (the nucleophile)
UBGG_N = "N1"           # Gly75 N -> Ub Arg74 C
UBGG_C = "C2"           # Gly76 C -> Cys91 SG (the thioester carbon)


def sumo2_lysines():
    """Native lysine positions in SUMO2, and their construct residue numbers."""
    return [(i + 1, i + 1 + OFFSET) for i, c in enumerate(SUMO2_NATIVE) if c == "K"]


def base(name, seeds):
    return {"name": name, "sequences": [], "modelSeeds": list(seeds),
            "dialect": "alphafold3", "version": 1, "bondedAtomPairs": []}


def add_prot(d, cid, seq):
    d["sequences"].append({"protein": {"id": cid, "sequence": seq}})


def add_lig(d, cid, code):
    d["sequences"].append({"ligand": {"id": cid, "ccdCodes": [code]}})


def job_xisok(pos, seeds, name):
    """XisoK on the substrate lysine, co-folded with UBE2W. No ubiquitin."""
    d = base(name, seeds)
    add_prot(d, "A", SUMO2_CONSTRUCT)
    add_prot(d, "B", UBE2W)
    add_lig(d, "L", "LIG-1")
    d["bondedAtomPairs"] = [[["A", pos, "NZ"], ["L", 1, LIG_LINK]]]
    return d, {"LIG-1": {LIG_LINK: 1}}


def job_lyscontrol(seeds, name):
    """Unmodified substrate + UBE2W. No ligand, no bonds."""
    d = base(name, seeds)
    add_prot(d, "A", SUMO2_CONSTRUCT)
    add_prot(d, "B", UBE2W)
    del d["bondedAtomPairs"]
    return d, {}


def job_charged(pos, seeds, name):
    """Ub thioester on Cys91; substrate XisoK amine left free.

    Ubiquitin is split 1-74 as a protein chain with Gly75-Gly76 supplied as the
    custom UBGG ligand, because AF3 3.0.1 discards polymer-polymer bonds. UBGG
    has no OXT, so its C2 is a genuine sp2 thioester carbon once bonded to SG.
    """
    d = base(name, seeds)
    add_prot(d, "A", SUMO2_CONSTRUCT)
    add_prot(d, "B", UBE2W)
    add_prot(d, "U", UBIQUITIN[:74])
    add_lig(d, "L", "LIG-1")
    add_lig(d, "G", "UBGG")
    d["bondedAtomPairs"] = [
        [["A", pos, "NZ"], ["L", 1, LIG_LINK]],
        [["U", 74, "C"], ["G", 1, UBGG_N]],
        [["G", 1, UBGG_C], ["B", UBE2W_CAT_CYS, "SG"]],
    ]
    return d, {"LIG-1": {LIG_LINK: 1}, "UBGG": {UBGG_N: 1, UBGG_C: 1}}


def job_tetrahedral(pos, seeds, name):
    """The tetrahedral intermediate, declared DELIBERATELY.

    Same as `charged` plus an explicit bond from the XisoK amine to the thioester
    carbon. That carbon then legitimately has four substituents (=O becomes an
    alkoxide in the real intermediate; here it stays a double bond, which is the
    approximation) -- so this job is EXPECTED to fail the valence check unless the
    ligand is rebuilt with an sp3 carbon. Reported honestly rather than forced.
    """
    d = base(name, seeds)
    add_prot(d, "A", SUMO2_CONSTRUCT)
    add_prot(d, "B", UBE2W)
    add_prot(d, "U", UBIQUITIN[:74])
    add_lig(d, "L", "LIG-1")
    add_lig(d, "G", "UBGT")     # sp3 carbon, alkoxide O -- see make_ccd_v2.py
    d["bondedAtomPairs"] = [
        [["A", pos, "NZ"], ["L", 1, LIG_LINK]],
        [["U", 74, "C"], ["G", 1, UBGG_N]],
        [["G", 1, UBGG_C], ["B", UBE2W_CAT_CYS, "SG"]],
        [["L", 1, LIG_NUC], ["G", 1, UBGG_C]],
    ]
    return d, {"LisoK_openN": {LIG_LINK: 1, LIG_NUC: 1},
               "UBGT": {UBGG_N: 1, UBGG_C: 2}}


def job_product(pos, seeds, name):
    """Ub C-terminus joined to the XisoK amine: the reaction has happened.

    Full ubiquitin 1-76 as a protein chain -- the bond Ub Gly76 C -> LisoK N01 is
    polymer->ligand, which AF3 accepts, so no split is needed here.
    """
    d = base(name, seeds)
    add_prot(d, "A", SUMO2_CONSTRUCT)
    add_prot(d, "B", UBE2W)
    add_prot(d, "U", UBIQUITIN)
    add_lig(d, "L", "LIG-1")
    d["bondedAtomPairs"] = [
        [["A", pos, "NZ"], ["L", 1, LIG_LINK]],
        [["U", 76, "C"], ["L", 1, LIG_NUC]],
    ]
    # LisoK_openN: the alpha-amine has one H removed so it can accept Ub's
    # C-terminus. In the product that nitrogen genuinely is a secondary amide.
    return d, {"LisoK_openN": {LIG_LINK: 1, LIG_NUC: 1}}


def ccd_filename(code):
    """Map an ext-map key to its CCD file. All XisoK files declare comp id LIG-1,
    so the key names the FILE (which variant / which valence pattern), not the
    comp id -- that is how one bond spec serves every variant unchanged."""
    if code == "LIG-1":
        return "LisoK_userCCD.cif"
    return f"{code}_userCCD.cif"


def check_job(d, ext, ccddir, label):
    """Verify sequences, bond targets and ligand valences before writing."""
    ch = {}
    for s in d["sequences"]:
        k = list(s)[0]
        ch[s[k]["id"]] = (k, s[k])
    # every bond target must exist and be the right residue type
    for p1, p2 in d.get("bondedAtomPairs", []):
        for cid, resi, atom in (p1, p2):
            if cid not in ch:
                raise ValueError(f"{label}: bond references missing chain {cid}")
            kind, body = ch[cid]
            if kind == "protein":
                aa = body["sequence"][resi - 1]
                if atom == "NZ" and aa != "K":
                    raise ValueError(f"{label}: {cid}:{resi} is {aa}, not Lys "
                                     f"(NZ bond target)")
                if atom == "SG" and aa != "C":
                    raise ValueError(f"{label}: {cid}:{resi} is {aa}, not Cys "
                                     f"(SG bond target)")
    # no polymer-polymer bonds: AF3 3.0.1 silently discards them
    for p1, p2 in d.get("bondedAtomPairs", []):
        if ch[p1[0]][0] == "protein" and ch[p2[0]][0] == "protein":
            raise ValueError(f"{label}: polymer-polymer bond {p1}-{p2} would be "
                             f"SILENTLY DISCARDED by AF3 3.0.1")
    # ligand valences, counting the external bonds this job adds
    problems = []
    for code, extmap in ext.items():
        path = os.path.join(ccddir, ccd_filename(code))
        if not os.path.exists(path):
            raise ValueError(f"{label}: missing CCD {path}")
        atoms, v = valences(open(path).read(), extmap)
        from ccd_valence import MAX_VALENCE
        for a, tot in v.items():
            lim = MAX_VALENCE.get(atoms[a])
            if lim is not None and tot > lim:
                problems.append(f"{code}/{a} ({atoms[a]}) order {tot:g} > {lim}")
    return problems


def embed_ccd(d, ext, ccddir):
    """Concatenate every needed CCD block into the job's userCCD field."""
    blocks = []
    for code in ext:
        fn = ccd_filename(code)
        blocks.append(open(os.path.join(ccddir, fn)).read().rstrip() + "\n")
    if blocks:
        d["userCCD"] = "\n".join(blocks)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--ccd", required=True)
    ap.add_argument("--seeds", type=int, nargs="+", default=list(range(1, 26)))
    ap.add_argument("--sweep", action="store_true",
                    help="also emit an xisok job for every SUMO2 lysine")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    k11_native, k11_pos = 11, 11 + OFFSET
    assert SUMO2_CONSTRUCT[k11_pos - 1] == "K", "K11 construct position is not Lys"

    specs = [
        ("sumo2_k11_xisok_ube2w", job_xisok(k11_pos, a.seeds,
                                            "sumo2_k11_xisok_ube2w")),
        ("sumo2_k11_lyscontrol_ube2w", job_lyscontrol(a.seeds,
                                                      "sumo2_k11_lyscontrol_ube2w")),
        ("sumo2_k11_charged_ube2w", job_charged(k11_pos, a.seeds,
                                                "sumo2_k11_charged_ube2w")),
        ("sumo2_k11_tetrahedral_ube2w", job_tetrahedral(k11_pos, a.seeds,
                                                        "sumo2_k11_tetrahedral_ube2w")),
        ("sumo2_k11_product_ube2w", job_product(k11_pos, a.seeds,
                                                "sumo2_k11_product_ube2w")),
    ]
    if a.sweep:
        for nat, pos in sumo2_lysines():
            if nat == k11_native:
                continue
            n = f"sumo2_k{nat}_xisok_ube2w"
            specs.append((n, job_xisok(pos, a.seeds, n)))

    manifest, failed = [], []
    for name, (d, ext) in specs:
        problems = check_job(d, ext, a.ccd, name)
        embed_ccd(d, ext, a.ccd)
        nb = len(d.get("bondedAtomPairs", []))
        if problems:
            failed.append((name, problems))
            print(f"  {name:32s} {nb} bonds  VALENCE PROBLEM: {problems}")
            continue
        with open(os.path.join(a.outdir, f"{name}.json"), "w") as fh:
            json.dump(d, fh, indent=1)
        manifest.append(dict(name=name, bonds=nb, seeds=len(a.seeds),
                             chains=len(d["sequences"])))
        print(f"  {name:32s} {nb} bonds  {len(d['sequences'])} chains  valence OK")

    with open(os.path.join(a.outdir, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=1)
    print(f"\nwrote {len(manifest)} jobs to {a.outdir}/")
    if failed:
        print(f"{len(failed)} job(s) NOT written because their bonds would "
              f"over-fill a ligand atom:")
        for n, p in failed:
            print(f"  {n}: {'; '.join(p)}")


if __name__ == "__main__":
    main()
