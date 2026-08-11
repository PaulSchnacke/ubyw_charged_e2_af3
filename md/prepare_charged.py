#!/usr/bin/env python3
"""Convert the AF3 charged-state co-fold into an Amber-ready PDB.

Input : the AF3 charged model (SUMO2 + UBE2W + Ub 1-74 + UBGG dipeptide + XisoK)
Output: a PDB where
          * chain A residue 12 is LYQ  (lysine + XisoK acyl fused, validated residue)
          * chain U is a CONTINUOUS Ub 1-76: the UBGG ligand's two glycine units are
            promoted to standard GLY 75 and GLY 76 of the ubiquitin chain
          * chain B residue 91 is CYX  (cysteine with no SG hydrogen)
        and a tleap script that creates the thioester as an explicit inter-chain
        bond Gly76.C -> Cys91.SG.

WHY THE UBGG LIGAND EXISTS AT ALL. AF3 3.0.1 silently discards covalent bonds
between two POLYMER chains (structure_cleaning.py: "Reducing number of bonds ...
N are polymer-polymer"), exits 0, and leaves the atoms tens of angstroms apart.
The workaround was to end the ubiquitin protein chain at Arg74 and carry Gly75 and
Gly76 as a two-residue LIGAND (UBGG), because protein->ligand bonds ARE kept.
That constraint belongs to AF3 only. Amber has no such restriction, so here the
split is REVERSED: the glycines rejoin their own chain and the only remaining
inter-chain bond is the thioester itself.

UBGG atom mapping (asserted below, not assumed):
    N1 CA1 C1 O1  ->  GLY 75  N CA C O
    N2 CA2 C2 O2  ->  GLY 76  N CA C O
C2 is the thioester carbonyl carbon; after the merge it is simply GLY76.C.

Ub(1-75) VARIANT (--ub75). For the engineering question -- does UBE2W tolerate a
ubiquitin one residue shorter -- the thioester sits on Gly75 instead of Gly76. Pass
--ub75 and the second glycine is dropped, leaving Ub 1-75 whose GLY75.C carries the
thioester. Everything else is identical, so the two systems differ by exactly one
glycine.

Usage:
    python prepare_charged.py IN.cif OUT.pdb [--ub75] [--site 12]
"""
import sys

import gemmi

# XisoK ligand -> LYQ, the isopeptide residue already charge-validated against ALY.
LIG2LYQ = {"C01": "CH", "O01": "OH", "C02": "CI", "N01": "NI",
           "C03": "CJ", "C04": "CK", "C05": "CL", "C06": "CM"}

# UBGG dipeptide ligand -> two standard glycines.
UBGG2GLY = {75: {"N1": "N", "CA1": "CA", "C1": "C", "O1": "O"},
            76: {"N2": "N", "CA2": "CA", "C2": "C", "O2": "O"}}

SUB_CHAIN, ENZ_CHAIN, LIG_CHAIN, UB_CHAIN, GLY_CHAIN = "A", "B", "L", "U", "G"
CATALYTIC_CYS = 91
UB_LAST_PROTEIN_RES = 74          # Arg74; 75 and 76 arrive via the ligand


def main(cif, out, *flags):
    ub75 = "--ub75" in flags
    site = 12
    if "--site" in flags:
        site = int(flags[list(flags).index("--site") + 1])

    st = gemmi.read_structure(cif)
    st.setup_entities()
    model = st[0]
    chains = {ch.name: ch for ch in model}
    for need in (SUB_CHAIN, ENZ_CHAIN, LIG_CHAIN, UB_CHAIN, GLY_CHAIN):
        if need not in chains:
            sys.exit(f"missing chain {need}; found {sorted(chains)}")

    # ---- 1. fuse the XisoK ligand into the substrate lysine as LYQ -------------
    lig = [r for r in chains[LIG_CHAIN] if r.name == "LIG-1"]
    if len(lig) != 1:
        sys.exit(f"expected exactly one LIG-1, found {len(lig)}")
    lig = lig[0]
    lig_atoms = {a.name: a for a in lig}
    if set(LIG2LYQ) - set(lig_atoms):
        sys.exit(f"ligand missing atoms: {sorted(set(LIG2LYQ) - set(lig_atoms))}")

    target = [r for r in chains[SUB_CHAIN] if r.seqid.num == site]
    if len(target) != 1 or target[0].name != "LYS":
        sys.exit(f"chain {SUB_CHAIN} residue {site} is not a unique LYS")
    target = target[0]
    nz = target.find_atom("NZ", "*")
    if nz is None:
        sys.exit("target lysine has no NZ")
    d = nz.pos.dist(lig_atoms["C01"].pos)
    if d > 2.0:
        sys.exit(f"NZ-C01 = {d:.2f} A -- isopeptide bond not formed in this model")
    print(f"  isopeptide NZ-C01 = {d:.2f} A")
    target.name = "LYQ"
    for old, new in LIG2LYQ.items():
        a = lig_atoms[old]
        a.name = new
        target.add_atom(a)
    print(f"  LYQ built with {len(target)} heavy atoms")

    # ---- 2. promote the UBGG glycines into the ubiquitin chain ------------------
    ubgg = [r for r in chains[GLY_CHAIN] if r.name == "UBGG"]
    if len(ubgg) != 1:
        sys.exit(f"expected exactly one UBGG, found {len(ubgg)}")
    ubgg = ubgg[0]
    g_atoms = {a.name: a for a in ubgg}
    expected = {n for m in UBGG2GLY.values() for n in m}
    if set(g_atoms) != expected:
        sys.exit(f"UBGG atoms {sorted(g_atoms)} != expected {sorted(expected)}")

    ub = chains[UB_CHAIN]
    last = max(r.seqid.num for r in ub)
    if last != UB_LAST_PROTEIN_RES:
        sys.exit(f"ubiquitin chain ends at {last}, expected {UB_LAST_PROTEIN_RES}")
    # Arg74 is about to become internal, so its C-terminal OXT must go or the
    # carbonyl carbon would carry both OXT and the Gly75 amide nitrogen.
    r74 = [r for r in ub if r.seqid.num == UB_LAST_PROTEIN_RES][0]
    for nm in ("OXT", "OT2"):
        a = r74.find_atom(nm, "*")
        if a is not None:
            r74.remove_atom(nm, "\0", gemmi.Element("X"))
            print(f"  removed {nm} from Arg74 (now internal)")

    keep = [75] if ub75 else [75, 76]
    for num in keep:
        res = gemmi.Residue()
        res.name = "GLY"
        res.seqid = gemmi.SeqId(num, " ")
        for old, new in UBGG2GLY[num].items():
            a = gemmi.Atom()
            a.name = new
            src = g_atoms[old]
            a.pos, a.element, a.b_iso, a.occ = src.pos, src.element, src.b_iso, src.occ
            res.add_atom(a)
        ub.add_residue(res)
    thio_res = keep[-1]
    print(f"  ubiquitin now 1-{thio_res} "
          f"({'Ub(1-75) VARIANT' if ub75 else 'wild-type Ub(1-76)'})")

    # Verify the new peptide bonds are geometrically real rather than assumed.
    prev = r74
    for num in keep:
        cur = [r for r in ub if r.seqid.num == num][0]
        dd = prev.find_atom("C", "*").pos.dist(cur.find_atom("N", "*").pos)
        if not 1.0 < dd < 1.8:
            sys.exit(f"peptide bond {prev.seqid.num}C-{num}N is {dd:.2f} A")
        print(f"  peptide bond {prev.seqid.num}.C-{num}.N = {dd:.2f} A")
        prev = cur

    # ---- 3. the thioester: measure it before writing anything ------------------
    cys = [r for r in chains[ENZ_CHAIN] if r.seqid.num == CATALYTIC_CYS]
    if len(cys) != 1 or cys[0].name != "CYS":
        sys.exit(f"chain {ENZ_CHAIN} residue {CATALYTIC_CYS} is not CYS")
    cys = cys[0]
    sg = cys.find_atom("SG", "*")
    thio_c = [r for r in ub if r.seqid.num == thio_res][0].find_atom("C", "*")
    dt = thio_c.pos.dist(sg.pos)
    print(f"  thioester GLY{thio_res}.C - CYS91.SG = {dt:.2f} A "
          f"({'as modelled' if dt < 2.2 else 'NOT BONDED in the input model'})")
    # CYX = cysteine with no SG hydrogen, so tleap will not cap the sulfur.
    cys.name = "CYX"

    # Drop the now-empty ligand chains.
    for name in (LIG_CHAIN, GLY_CHAIN):
        model.remove_chain(name)
    st.setup_entities()

    st.write_pdb(out)
    remaining = [ch.name for ch in st[0]]
    print(f"  chains written: {remaining}")
    print(f"  wrote {out}")
    return dict(thio_res=thio_res, site=site, thio_dist=dt, chains=remaining)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], *sys.argv[3:])
