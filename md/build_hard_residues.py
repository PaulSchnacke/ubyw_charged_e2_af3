#!/usr/bin/env python3
"""Build model compounds for the three residues that need new ff19SB parameters.

Priority order, and an honest statement of how defensible each one is:

  LYP (product)     -- LYQ whose leucyl alpha-amine is ACYLATED by ubiquitin
                       Gly76. Chemically this is just a second amide, and amides
                       are the best-covered motif in any protein force field, so
                       AM1-BCC charges on a capped model compound are as good
                       here as they are for LYQ itself. DEFENSIBLE.

  LYT (thioester)   -- Gly-Cys thioester at the enzyme active site. A thioester
                       is NOT standard ff19SB chemistry: the C-S bond, the
                       reduced C=O bond order and the near-zero rotation barrier
                       about C-S are all poorly described by GAFF-style
                       transferable parameters. AM1-BCC will produce *a* charge
                       set; whether the resulting dynamics are quantitative is
                       unestablished. APPROXIMATE -- report as such.

  LYX (tetrahedral) -- the sp3 tetrahedral intermediate, with an alkoxide oxygen
                       on a carbon bearing both S and N. This is a transition-
                       state-like species. A classical force field has no
                       business describing it: it cannot break or form the bonds
                       that define the intermediate, so the trajectory can only
                       show whether the geometry is *mechanically stable* under
                       ff19SB, not whether it is chemically accessible.
                       NOT DEFENSIBLE as kinetics; run only as a stability probe.

Every compound is written as ACE-X-NME so prepgen can excise the caps and force
an integer residue charge, exactly as for LYQ. Cap atoms are named CY/OY/NN/HN/
CAY/CAT so they cannot collide with the residue's own backbone names -- the
collision that cost several submissions when LYQ was built.

Usage: python build_hard_residues.py OUTDIR
"""
import argparse
import os
import sys

from rdkit import Chem
from rdkit.Chem import AllChem

# ACE-<residue>-NME model compounds.
#   [C@@H] at the lysine alpha carbon and at the leucyl alpha carbon = L in both
#   cases (verified by reconstruction below, NOT by CIP letter -- the CIP code is
#   not a reliable L/D proxy and this project has been bitten by that twice).
COMPOUNDS = {
    # product: Ub Gly76 C(=O) acylates the leucyl alpha-amine
    "LYP": dict(
        smiles=("CC(=O)N[C@@H](CCCCNC(=O)[C@H](CC(C)C)NC(=O)CN)C(=O)NC"),
        note="isopeptide + Ub Gly76 amide on the leucyl amine (product)",
        confidence="defensible: both linkages are ordinary amides",
    ),
    # thioester: Gly-C(=O)-S-CH2- (the cysteine side chain)
    "LYT": dict(
        smiles=("CC(=O)N[C@@H](CSC(=O)CN)C(=O)NC"),
        note="Cys-S-C(=O)-Gly thioester at the active site",
        confidence="approximate: thioester is not standard ff19SB chemistry",
    ),
    # tetrahedral intermediate: sp3 C bearing S, N, O(-) and CH2
    "LYX": dict(
        smiles=("CC(=O)N[C@@H](CSC(O)(N)CN)C(=O)NC"),
        note="sp3 tetrahedral intermediate: C bonded to S, N, OH and CH2",
        confidence="NOT defensible as kinetics; stability probe only",
    ),
}

CAP_RENAME = {"CY", "OY", "NN", "HN", "CAY", "CAT"}

# L-configuration reference: the signed volume of the (N, C=O, side chain, H)
# tetrahedron at an alpha carbon is NEGATIVE for L-Lys, L-Cys and L-Leu alike.
L_SIGN = -1.0


def alpha_signs(mol):
    """[(atom_idx, signed_volume)] for every backbone alpha carbon."""
    import numpy as np
    conf = mol.GetConformer()
    out = []
    for a in mol.GetAtoms():
        if a.GetSymbol() != "C":
            continue
        nb = list(a.GetNeighbors())
        if sum(1 for n in nb if n.GetSymbol() == "H") != 1 or len(nb) != 4:
            continue
        has_n = any(n.GetSymbol() == "N" for n in nb)
        has_co = any(n.GetSymbol() == "C" and
                     any(x.GetSymbol() == "O" and
                         mol.GetBondBetweenAtoms(n.GetIdx(), x.GetIdx()
                                                 ).GetBondTypeAsDouble() == 2
                         for x in n.GetNeighbors()) for n in nb)
        if not (has_n and has_co):
            continue

        def rank(x):
            if x.GetSymbol() == "N":
                return 0
            if x.GetSymbol() == "H":
                return 3
            if any(y.GetSymbol() == "O" and
                   mol.GetBondBetweenAtoms(x.GetIdx(), y.GetIdx()
                                           ).GetBondTypeAsDouble() == 2
                   for y in x.GetNeighbors()):
                return 1
            return 2

        q = sorted(nb, key=rank)
        p = [np.array(conf.GetAtomPosition(x.GetIdx())) for x in q]
        out.append((a.GetIdx(),
                    float(np.dot(np.cross(p[1] - p[0], p[2] - p[0]), p[3] - p[0]))))
    return out


def build(code, spec, outdir):
    mol = Chem.MolFromSmiles(spec["smiles"])
    if mol is None:
        raise ValueError(f"{code}: SMILES failed to parse")
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE) != 0:
        raise ValueError(f"{code}: 3D embedding failed")
    AllChem.MMFFOptimizeMolecule(mol, maxIters=4000)

    # verify stereochemistry structurally: every declared centre must be assigned,
    # and the count must match what the SMILES asked for
    Chem.AssignStereochemistryFrom3D(mol)
    # STRUCTURAL L-check at every alpha carbon: signed volume of the
    # (N, carbonyl-C, side-chain-C, H) tetrahedron must match the L reference.
    # The CIP letter is NOT the test -- L-Cys is (R) while L-Lys is (S), same
    # geometry, different letter. An earlier version of this file shipped a
    # D-leucine in LYP that the letter check did not catch.
    bad = [i for i, s in alpha_signs(mol) if s > 0]
    if bad:
        raise ValueError(f"{code}: alpha carbon(s) {bad} have D configuration "
                         f"(signed volume > 0; L reference is negative)")
    centres = Chem.FindMolChiralCenters(mol, includeUnassigned=True,
                                        useLegacyImplementation=False)
    unassigned = [c for c in centres if c[1] in ("?", None)]
    if unassigned:
        raise ValueError(f"{code}: unassigned stereocentres {unassigned}")

    # Name atoms and split into ACE / <code> / NME, then write via gemmi.
    #
    # WHY: RDKit's PDBWriter emits generic names (C1, C2, N1...) in residue UNL.
    # prepgen's mainchain file refers to atoms BY NAME, so "MAIN_CHAIN CA" matched
    # nothing, prepgen built a garbage tree, and tleap died with
    #   "Atom C2: Illegal chain specifier [X] in PREP file"
    # after antechamber and parmchk2 had both reported success. Names are
    # load-bearing, exactly as in build_lyq_model.py.
    info = assign_names(mol, code)
    out = os.path.join(outdir, f"{code}_model.pdb")
    write_pdb(mol, info, code, out)
    heavy = sum(1 for a in mol.GetAtoms() if a.GetSymbol() != "H")
    return out, heavy, centres


def assign_names(mol, code):
    """{atom_idx: (pdb_name, residue_name)} with ACE/NME caps named by convention.

    The two acetyl caps are identified structurally: each is a methyl carbon bonded
    to a carbonyl carbon (ACE) or to the amide nitrogen (NME) at the chain termini.
    """
    info = {}
    # backbone: the alpha carbon is the C with one H, an N, and a carbonyl C
    alphas = [i for i, _ in alpha_signs(mol)]
    if not alphas:
        raise ValueError(f"{code}: no alpha carbon found")
    ca = alphas[0]                      # the LYSINE/CYS alpha carbon
    a = mol.GetAtomWithIdx(ca)
    n_bb = next(x for x in a.GetNeighbors() if x.GetSymbol() == "N")
    c_bb = next(x for x in a.GetNeighbors()
                if x.GetSymbol() == "C" and
                any(y.GetSymbol() == "O" and
                    mol.GetBondBetweenAtoms(x.GetIdx(), y.GetIdx()
                                            ).GetBondTypeAsDouble() == 2
                    for y in x.GetNeighbors()))
    o_bb = next(y for y in c_bb.GetNeighbors()
                if y.GetSymbol() == "O" and
                mol.GetBondBetweenAtoms(c_bb.GetIdx(), y.GetIdx()
                                        ).GetBondTypeAsDouble() == 2)
    info[ca] = ("CA", code)
    info[n_bb.GetIdx()] = ("N", code)
    info[c_bb.GetIdx()] = ("C", code)
    info[o_bb.GetIdx()] = ("O", code)

    # ACE cap: the carbonyl carbon on the far side of the backbone N, plus its
    # methyl and oxygen
    ace_c = next((x for x in n_bb.GetNeighbors()
                  if x.GetIdx() != ca and x.GetSymbol() == "C"), None)
    if ace_c is None:
        raise ValueError(f"{code}: no ACE cap found on the backbone N")
    info[ace_c.GetIdx()] = ("CY", "ACE")
    for y in ace_c.GetNeighbors():
        if y.GetSymbol() == "O":
            info[y.GetIdx()] = ("OY", "ACE")
        elif y.GetSymbol() == "C" and y.GetIdx() != n_bb.GetIdx():
            info[y.GetIdx()] = ("CAY", "ACE")

    # NME cap: the N on the far side of the backbone carbonyl, plus its methyl
    nme_n = next((x for x in c_bb.GetNeighbors()
                  if x.GetSymbol() == "N"), None)
    if nme_n is None:
        raise ValueError(f"{code}: no NME cap found on the backbone C")
    info[nme_n.GetIdx()] = ("NN", "NME")
    for y in nme_n.GetNeighbors():
        if y.GetSymbol() == "C" and y.GetIdx() != c_bb.GetIdx():
            info[y.GetIdx()] = ("CAT", "NME")

    # everything else belongs to the residue; name by element with a counter
    counters = {}
    for at in mol.GetAtoms():
        i = at.GetIdx()
        if i in info:
            continue
        s = at.GetSymbol()
        if s == "H":
            continue                    # hydrogens named after their heavy parent
        counters[s] = counters.get(s, 0) + 1
        info[i] = (f"{s}{counters[s]}", code)

    # hydrogens: H<parent name>, numbered when a parent carries several. Cap
    # hydrogens inherit the cap residue so prepgen's OMIT_NAME list catches them.
    # Number hydrogens PER PARENT ATOM INDEX, not per parent name, and always
    # number when the parent carries more than one H. Keying on the name collided
    # whenever two parents shared a name stem (HCAY/HCAT methyls), and truncating
    # to 4 characters then merged distinct names -- the duplicate check below
    # caught both.
    hcount = {}
    for at in mol.GetAtoms():
        if at.GetSymbol() != "H":
            continue
        par = at.GetNeighbors()[0]
        pname, prn = info[par.GetIdx()]
        nh = sum(1 for x in par.GetNeighbors() if x.GetSymbol() == "H")
        stem = "H" + pname.lstrip("CNOS") if pname not in ("N", "C") else "H"
        if pname == "N":
            nm = "H"
        elif nh == 1:
            nm = f"H{pname}"
        else:
            k = par.GetIdx()
            hcount[k] = hcount.get(k, 0) + 1
            nm = f"H{pname}{hcount[k]}"
        if len(nm) > 4:
            # PDB atom names are 4 columns; keep the parent stem and the index
            nm = (f"H{pname[-2:]}{hcount.get(par.GetIdx(), '')}")[:4]
        info[at.GetIdx()] = (nm, prn)

    names = [v[0] for v in info.values()]
    dupes = {n for n in names if names.count(n) > 1}
    if dupes:
        raise ValueError(f"{code}: duplicate atom names {sorted(dupes)} -- prepgen "
                         f"and parmchk2 both resolve atoms by name")
    return info


def write_pdb(mol, info, code, out):
    """Write ACE-<code>-NME via gemmi so the fixed PDB columns are exact.

    antechamber parses PDB by column position and rejects a file whose coordinate
    fields are off by even one space; a hand-rolled format string caused exactly
    that failure when LYQ was built.
    """
    import gemmi
    conf = mol.GetConformer()
    st = gemmi.Structure()
    st.add_model(gemmi.Model("1"))
    ch = gemmi.Chain("A")
    for seq, resname in enumerate(("ACE", code, "NME"), start=1):
        res = gemmi.Residue()
        res.name = resname
        res.seqid = gemmi.SeqId(seq, " ")
        res.het_flag = "A"
        for idx, (aname, rn) in sorted(info.items(), key=lambda kv: kv[0]):
            if rn != resname:
                continue
            p = conf.GetAtomPosition(idx)
            at = gemmi.Atom()
            at.name = aname
            at.element = gemmi.Element(mol.GetAtomWithIdx(idx).GetSymbol())
            at.pos = gemmi.Position(p.x, p.y, p.z)
            at.occ = 1.0
            at.b_iso = 0.0
            res.add_atom(at)
        ch.add_residue(res)
    st[0].add_chain(ch)
    st.setup_entities()
    with open(out, "w") as fh:
        fh.write(st.make_pdb_string())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--only", nargs="*")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    for code, spec in COMPOUNDS.items():
        if a.only and code not in a.only:
            continue
        try:
            path, heavy, centres = build(code, spec, a.outdir)
            print(f"{code}: {heavy} heavy atoms, {len(centres)} stereocentre(s) "
                  f"{[c[1] for c in centres]} -> {os.path.basename(path)}")
            print(f"      {spec['note']}")
            print(f"      confidence: {spec['confidence']}")
        except Exception as e:
            print(f"{code}: BUILD FAILED -- {type(e).__name__}: {e}")


if __name__ == "__main__":
    main()
