#!/usr/bin/env python3
"""Convert the AF3 co-fold into an Amber-ready PDB with LYQ as one residue.

Input:  the SUMO2(K11LisoK) + UBE2W co-fold (chain A substrate, B enzyme,
        L the LisoK acyl ligand).
Output: a PDB where chain A residue 12 is LYQ -- lysine and acyl fused into a
        single residue -- so tleap builds the isopeptide bond from the LYQ
        template instead of needing a manual bond command with hand-guessed
        parameters.

Atom mapping (ligand -> LYQ), fixed here and asserted:
    C01 -> CH    isopeptide carbonyl carbon (was bonded to lysine NZ)
    O01 -> OH    isopeptide carbonyl oxygen
    C02 -> CI    leucyl alpha-carbon
    N01 -> NI    leucyl free alpha-amine (THE NUCLEOPHILE)
    C03 -> CJ    leucyl C-beta
    C04 -> CK    leucyl C-gamma
    C05 -> CM1   leucyl C-delta1
    C06 -> CM2   leucyl C-delta2

WHAT "NON-CHARGED UBE2W" MEANS HERE. UBE2W's catalytic Cys151... no: the
catalytic cysteine is Cys91. The C-terminal Cys151 is a separate residue, and
the C-terminal ~10 residues are the flexible tail (AF3 pLDDT falls from 91 for
residues 1-140 to 41 at Cys151). Uncharged means no ubiquitin is attached at
Cys91 -- both cysteines are free thiols (CYS, not CYX/CYM), which is the
uncharged E2 resting state. That is what makes this a "test" run: the amine is
attacking an empty active site, which rounds 2 and 3 showed is permissive.

Usage: python prepare_system.py COFOLD.cif OUT.pdb
"""
import sys
import gemmi

LIG2LYQ = {"C01": "CH", "O01": "OH", "C02": "CI", "N01": "NI",
           "C03": "CJ", "C04": "CK", "C05": "CM1", "C06": "CM2"}
SUB_CHAIN, ENZ_CHAIN, LIG_CHAIN = "A", "B", "L"
# Modified-site residue number in CONSTRUCT numbering (the construct carries an
# N-terminal Pro, so native K11 -> 12 and native K21 -> 22). Overridable on the
# command line so the same script serves every site in the sweep.
MOD_RESNUM = 12
CATALYTIC_CYS = 91


def main(cif, out, mod_resnum=None):
    global MOD_RESNUM
    if mod_resnum is not None:
        MOD_RESNUM = int(mod_resnum)
    st = gemmi.read_structure(cif)
    st.setup_entities()
    model = st[0]

    chains = {ch.name: ch for ch in model}
    for need in (SUB_CHAIN, ENZ_CHAIN, LIG_CHAIN):
        if need not in chains:
            sys.exit(f"chain {need} missing -- chains present: {sorted(chains)}")

    lig = [r for r in chains[LIG_CHAIN] if r.name == "LIG-1"]
    if len(lig) != 1:
        sys.exit(f"expected exactly one LIG-1, found {len(lig)}")
    lig = lig[0]
    lig_atoms = {a.name: a for a in lig}
    missing = set(LIG2LYQ) - set(lig_atoms)
    if missing:
        sys.exit(f"ligand lacks expected atoms: {sorted(missing)}")
    extra = set(lig_atoms) - set(LIG2LYQ)
    if extra:
        sys.exit(f"ligand has unmapped atoms {sorted(extra)} -- refusing to "
                 "silently drop them")

    target = [r for r in chains[SUB_CHAIN] if r.seqid.num == MOD_RESNUM]
    if len(target) != 1 or target[0].name != "LYS":
        sys.exit(f"chain {SUB_CHAIN} residue {MOD_RESNUM} is "
                 f"{target[0].name if target else 'absent'}, expected LYS")
    target = target[0]

    cys = [r for r in chains[ENZ_CHAIN] if r.seqid.num == CATALYTIC_CYS]
    if len(cys) != 1 or cys[0].name != "CYS":
        sys.exit(f"chain {ENZ_CHAIN} residue {CATALYTIC_CYS} is not CYS")

    # sanity: the isopeptide bond must actually be formed in this model, or we
    # would be building a topology the coordinates do not support
    nz = target.find_atom("NZ", "*")
    if nz is None:
        sys.exit("target lysine has no NZ")
    # AF3 writes covalent bonds SHORT: across all 350 models in this project the
    # C-N isopeptide bond is 1.07 +- 0.09 A (range 0.77-1.36) against an ideal
    # amide 1.33 A. So the acceptance window is wide on the low side. This is
    # harmless as an MD starting point -- ff19SB restores the correct bond length
    # during minimisation -- but it must be checked, because a value outside this
    # range would mean the acyl group is not attached to the lysine at all and we
    # would be fusing two residues that the coordinates keep apart.
    d = nz.pos.dist(lig_atoms["C01"].pos)
    if not 0.6 < d < 2.0:
        sys.exit(f"NZ-C01 distance is {d:.2f} A -- isopeptide bond not formed")
    print(f"  isopeptide NZ-C01 = {d:.2f} A "
          f"(bond present; AF3 runs short, ff19SB will relax it to ~1.33 A)")

    # fuse: rename the lysine to LYQ and append the renamed ligand atoms
    target.name = "LYQ"
    for old, new in LIG2LYQ.items():
        a = gemmi.Atom()
        a.name = new
        a.element = lig_atoms[old].element
        a.pos = lig_atoms[old].pos
        a.b_iso = lig_atoms[old].b_iso
        a.occ = 1.0
        target.add_atom(a)
    print(f"  LYQ now has {len(target)} heavy atoms "
          f"(9 from Lys + {len(LIG2LYQ)} from the acyl group)")

    # drop the now-empty ligand chain
    model.remove_chain(LIG_CHAIN)

    # strip hydrogens: tleap adds them with ff19SB geometry, and AF3's are
    # absent anyway; also drop anything not protein (this test system has none)
    st.remove_hydrogens()
    st.remove_ligands_and_waters()   # no-op here, guards against surprises

    remaining = {ch.name: len(ch) for ch in st[0]}
    if set(remaining) != {SUB_CHAIN, ENZ_CHAIN}:
        sys.exit(f"unexpected chains after cleanup: {remaining}")
    tgt = [r for r in st[0][SUB_CHAIN] if r.seqid.num == MOD_RESNUM][0]
    if tgt.name != "LYQ" or tgt.find_atom("NI", "*") is None:
        sys.exit("LYQ fusion lost during cleanup")

    st.setup_entities()
    doc = st.make_pdb_string()
    with open(out, "w") as fh:
        fh.write(doc)
    print(f"  chains: {remaining}")
    print(f"  wrote {out}")


if __name__ == "__main__":
    # usage: prepare_system.py IN.cif OUT.pdb [MOD_RESNUM]
    main(*sys.argv[1:4])
