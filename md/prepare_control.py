#!/usr/bin/env python3
"""Convert the unmodified-lysine control co-fold to an Amber-ready PDB.

No custom residue: the substrate lysine stays a standard LYS, so tleap builds
this from ff19SB alone. prepare_system.py deliberately refuses this input (it
requires the ligand chain), which is why the control gets its own script rather
than a flag that weakens that assertion.

IMPORTANT about the starting geometry. The other systems start from near-attack
poses. This one cannot: across 200 AF3 models of the control, the unmodified
construct residue 12 (native K11) never comes closer than 12.8 A to Cys91 SG.
That is a finding, not a limitation to engineer around -- the modification is
what recruits the site. So this trajectory asks a different question: starting
from where AF3 actually puts an unmodified lysine, does it approach at all?

Usage: python prepare_control.py IN.cif OUT.pdb [MOD_RESNUM]
"""
import sys

import gemmi

SUB_CHAIN, ENZ_CHAIN = "A", "B"
CATALYTIC_CYS = 91


def main(cif, out, mod_resnum=12):
    mod_resnum = int(mod_resnum)
    st = gemmi.read_structure(cif)
    st.setup_entities()
    chains = {ch.name: ch for ch in st[0]}
    for need in (SUB_CHAIN, ENZ_CHAIN):
        if need not in chains:
            sys.exit(f"chain {need} missing -- chains present: {sorted(chains)}")
    if "L" in chains:
        sys.exit("a ligand chain L is present -- this is not the unmodified "
                 "control; use prepare_system.py")

    tgt = [r for r in chains[SUB_CHAIN] if r.seqid.num == mod_resnum]
    if not tgt or tgt[0].name != "LYS":
        sys.exit(f"chain {SUB_CHAIN} residue {mod_resnum} is "
                 f"{tgt[0].name if tgt else 'absent'}, expected LYS")
    cys = [r for r in chains[ENZ_CHAIN] if r.seqid.num == CATALYTIC_CYS]
    if not cys or cys[0].name != "CYS":
        sys.exit(f"chain {ENZ_CHAIN} residue {CATALYTIC_CYS} is not CYS")

    nz = tgt[0].find_atom("NZ", "*")
    sg = cys[0].find_atom("SG", "*")
    if nz is None or sg is None:
        sys.exit("missing NZ or SG")
    d = nz.pos.dist(sg.pos)

    st.remove_ligands_and_waters()
    st.remove_hydrogens()
    st.setup_entities()
    remaining = {ch.name: len(ch) for ch in st[0]}
    if set(remaining) != {SUB_CHAIN, ENZ_CHAIN}:
        sys.exit(f"unexpected chains after cleanup: {remaining}")

    with open(out, "w") as fh:
        fh.write(st.make_pdb_string())
    print(f"  LYS{mod_resnum} NZ to Cys{CATALYTIC_CYS} SG = {d:.2f} A "
          f"(NOT near-attack; see the module docstring)")
    print(f"  chains: {remaining}")
    print(f"  wrote {out}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
