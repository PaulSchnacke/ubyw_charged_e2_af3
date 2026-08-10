#!/usr/bin/env python3
"""Build valence-correct CCDs for the XisoK acyl group and the Ub Gly-Gly tail.

Fixes the two bugs found in Julian's run:

  1. The acyl carbon must NOT carry a hydrogen. It was built from an aldehyde
     SMILES (O=C[...]), so RDKit added an aldehyde H to complete the valence and
     that H was written into the CCD. With the external isopeptide bond to Lys NZ
     the carbon then had 5 connections, and NZ-C01 collapsed to ~0.9 A.
     -> the H is removed and the valence is left open for the external bond.

  2. The Ub C-terminal Gly must NOT carry OXT. The standard GLY CCD keeps it, so
     the thioester carbon had C, O, OXT, CA before the Cys91 SG bond -- 5 again,
     and the carbon came out sp3 (angle sum 328 deg) where a thioester carbonyl
     must be planar sp2 (360 deg).
     -> a custom Gly-Gly CCD whose terminal carbon has only =O, CA and the
        incoming S.

Every file is checked with ccd_valence.check_ccd before it is written, with the
external bonds declared, so a valence error fails here rather than after the job
reaches the cluster.

Usage: python make_ccd_v2.py OUTDIR
"""
import argparse
import os
import sys

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from ccd_valence import check_ccd, report

# variant -> (SMILES written carbonyl-first, expected CIP at C-alpha, description)
# CIP letter is NOT L/D: L-Cys is (R) and beta-chloro-L-Ala is (R) because S and
# Cl outrank the carboxyl. Validated structurally elsewhere, not by the letter.
VARIANTS = {
    "LisoK":    ("O=C[C@@H](N)CC(C)C",      "S", "L-leucine"),
    "AisoK":    ("O=C[C@@H](N)C",           "S", "L-alanine"),
    "SisoK":    ("O=C[C@@H](N)CO",          "S", "L-serine"),
    "TisoK":    ("O=C[C@@H](N)[C@@H](C)O",  "S", "L-threonine"),
    "VisoK":    ("O=C[C@@H](N)C(C)C",       "S", "L-valine"),
    "CisoK":    ("O=C[C@@H](N)CS",          "R", "L-cysteine"),
    "PisoK":    ("O=C[C@@H]1CCCN1",         "S", "L-proline (2-deg amine)"),
    "aisoK":    ("O=C[C@H](N)C",            "R", "D-alanine"),
    "sisoK":    ("O=C[C@H](N)CO",           "R", "D-serine"),
    "pLisoK":   ("O=C[C@@H](N)CC1(C)N=N1",  "S", "photo-leucine"),
    "PrgGisoK": ("O=C[C@@H](N)CC#C",        "S", "propargylglycine"),
    "ClAisoK":  ("O=C[C@@H](N)CCl",         "R", "beta-chloro-L-alanine"),
}

AF3_KEYS = ["_chem_comp.id", "_chem_comp.name", "_chem_comp.type",
            "_chem_comp.formula", "_chem_comp.mon_nstd_parent_comp_id",
            "_chem_comp.pdbx_synonyms", "_chem_comp.formula_weight"]
ORDER = {Chem.BondType.SINGLE: "SING", Chem.BondType.DOUBLE: "DOUB",
         Chem.BondType.TRIPLE: "TRIP", Chem.BondType.AROMATIC: "AROM"}


def strip_acyl_hydrogen(mol):
    """Remove the H on the carbonyl carbon, leaving its valence open.

    THE FIX FOR BUG 1. The carbon keeps =O and the side chain; the fourth bond is
    supplied externally by AF3's bondedAtomPairs to the substrate lysine NZ.
    """
    # Select the ALDEHYDE carbon specifically: C=O bearing a hydrogen. Selecting
    # on "C=O with one carbon neighbour" is ambiguous in a peptide -- in Gly-Gly
    # the internal amide carbonyl also matches that pattern and has no H, so the
    # check fired on the wrong atom. The aldehyde H is the unambiguous marker,
    # since it is exactly the atom that has to go.
    cands = []
    for a in mol.GetAtoms():
        if a.GetSymbol() != "C":
            continue
        dbl_o = any(n.GetSymbol() == "O" and
                    mol.GetBondBetweenAtoms(a.GetIdx(), n.GetIdx()
                                            ).GetBondTypeAsDouble() == 2
                    for n in a.GetNeighbors())
        hs = [n.GetIdx() for n in a.GetNeighbors() if n.GetSymbol() == "H"]
        if dbl_o and len(hs) == 1:
            cands.append((a, hs))
    if len(cands) != 1:
        raise ValueError(f"expected exactly 1 aldehyde carbon (C=O with one H), "
                         f"found {len(cands)} -- the SMILES must present the "
                         f"reactive carbon as an aldehyde and nothing else")
    acyl, hs = cands[0]
    em = Chem.RWMol(mol)
    em.RemoveAtom(hs[0])
    out = em.GetMol()
    Chem.SanitizeMol(out, Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
    return out, acyl.GetIdx()


def name_atoms(mol, acyl_idx):
    """Deterministic names: acyl C=C01, its O=O01, alpha C=C02, its N=N01, rest C03+."""
    names, ci, oi, ni, hi = {}, 3, 2, 2, 1
    names[acyl_idx] = "C01"
    for n in mol.GetAtomWithIdx(acyl_idx).GetNeighbors():
        if n.GetSymbol() == "O":
            names[n.GetIdx()] = "O01"
    alpha = [n.GetIdx() for n in mol.GetAtomWithIdx(acyl_idx).GetNeighbors()
             if n.GetSymbol() == "C"]
    if len(alpha) != 1:
        raise ValueError("acyl carbon must have exactly one carbon neighbour")
    names[alpha[0]] = "C02"
    for n in mol.GetAtomWithIdx(alpha[0]).GetNeighbors():
        if n.GetSymbol() == "N":
            names[n.GetIdx()] = "N01"
    for a in mol.GetAtoms():
        if a.GetIdx() in names:
            continue
        s = a.GetSymbol()
        if s == "H":
            names[a.GetIdx()] = f"H{hi:02d}"; hi += 1
        elif s == "C":
            names[a.GetIdx()] = f"C{ci:02d}"; ci += 1
        elif s == "O":
            names[a.GetIdx()] = f"O{oi:02d}"; oi += 1
        elif s == "N":
            names[a.GetIdx()] = f"N{ni:02d}"; ni += 1
        else:
            names[a.GetIdx()] = f"{s.upper()}01"
    for a in mol.GetAtoms():
        a.SetProp("name", names[a.GetIdx()])
    return names


def to_cif(mol, code, name):
    conf = mol.GetConformer()
    L = [f"data_{code}", "#"]
    vals = {"_chem_comp.id": code, "_chem_comp.name": f"'{name}'",
            "_chem_comp.type": "non-polymer", "_chem_comp.formula": "?",
            "_chem_comp.mon_nstd_parent_comp_id": "?",
            "_chem_comp.pdbx_synonyms": "?", "_chem_comp.formula_weight": "?"}
    for k in AF3_KEYS:
        L.append(f"{k} {vals[k]}")
    L += ["#", "loop_", "_chem_comp_atom.comp_id", "_chem_comp_atom.atom_id",
          "_chem_comp_atom.type_symbol", "_chem_comp_atom.charge",
          "_chem_comp_atom.pdbx_leaving_atom_flag",
          "_chem_comp_atom.pdbx_model_Cartn_x_ideal",
          "_chem_comp_atom.pdbx_model_Cartn_y_ideal",
          "_chem_comp_atom.pdbx_model_Cartn_z_ideal"]
    for a in mol.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        L.append(f"{code} {a.GetProp('name')} {a.GetSymbol()} "
                 f"{a.GetFormalCharge()} N {p.x:.3f} {p.y:.3f} {p.z:.3f}")
    L += ["#", "loop_", "_chem_comp_bond.atom_id_1", "_chem_comp_bond.atom_id_2",
          "_chem_comp_bond.value_order", "_chem_comp_bond.pdbx_aromatic_flag"]
    for b in mol.GetBonds():
        L.append(f"{b.GetBeginAtom().GetProp('name')} "
                 f"{b.GetEndAtom().GetProp('name')} "
                 f"{ORDER[b.GetBondType()]} "
                 f"{'Y' if b.GetIsAromatic() else 'N'}")
    L.append("#")
    return "\n".join(L) + "\n"


def open_amine_valence(mol):
    """Drop one H from the alpha-amine so it can receive an external bond.

    Needed for the PRODUCT and TETRAHEDRAL species, where ubiquitin's C-terminus
    (or the thioester carbon) bonds to this nitrogen. As built the amine is a free
    -NH2 with 3 declared bonds, so an external bond would make nitrogen 4-valent.
    This is the same error class as bugs 1 and 2 -- caught by check_ccd here rather
    than by eye in a rendered model.

    Note the chemistry: in the product the nitrogen genuinely is a secondary amide
    -NH-, so removing one H is correct, not a fudge.
    """
    n = [a for a in mol.GetAtoms()
         if a.GetSymbol() == "N" and
         sum(1 for x in a.GetNeighbors() if x.GetSymbol() == "H") == 2]
    if len(n) != 1:
        raise ValueError(f"expected one -NH2 alpha-amine, found {len(n)}")
    h = [x.GetIdx() for x in n[0].GetNeighbors() if x.GetSymbol() == "H"][0]
    em = Chem.RWMol(mol)
    em.RemoveAtom(h)
    out = em.GetMol()
    Chem.SanitizeMol(out, Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
    return out


def build_xisok(code, smiles, expect_cip, desc, outdir, open_amine=False):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE)
    AllChem.MMFFOptimizeMolecule(mol)
    mol, acyl = strip_acyl_hydrogen(mol)
    ext = {"C01": 1}
    if open_amine:
        mol = open_amine_valence(mol)
        # re-locate the acyl carbon: indices shifted when the H was removed
        acyl = [a.GetIdx() for a in mol.GetAtoms()
                if a.GetSymbol() == "C" and
                any(n.GetSymbol() == "O" and
                    mol.GetBondBetweenAtoms(a.GetIdx(), n.GetIdx()
                                            ).GetBondTypeAsDouble() == 2
                    for n in a.GetNeighbors()) and
                sum(1 for n in a.GetNeighbors() if n.GetSymbol() == "C") == 1][0]
        ext["N01"] = 1
    # re-embed after the edit so the ideal coordinates are consistent
    AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    names = name_atoms(mol, acyl)
    suffix = " (amine open for the incoming Ub bond)" if open_amine else ""
    txt = to_cif(mol, "LIG-1", f"{desc} isopeptide moiety ({code}){suffix}")
    # THE GATE: C01 always receives the isopeptide bond to Lys NZ; N01 also
    # receives one in the product/tetrahedral species.
    check_ccd(txt, external=ext, label=code)
    fn = f"{code}_openN_userCCD.cif" if open_amine else f"{code}_userCCD.cif"
    path = os.path.join(outdir, fn)
    with open(path, "w") as fh:
        fh.write(txt)
    return path, mol.GetNumAtoms()


GLYGLY_NAME = "Ub Gly75-Gly76 with an open thioester valence"
GLYGLY_TET_NAME = ("Ub Gly75-Gly76, tetrahedral intermediate: C2 is sp3 with an "
                   "alkoxide O, accepting both SG and the incoming amine")


def build_glygly(outdir):
    """Gly-Gly dipeptide whose C-terminal carbon has NO OXT.

    THE FIX FOR BUG 2. Built from a SMILES whose C-terminus is an aldehyde, then
    the aldehyde H is stripped, so the terminal carbon has exactly =O, CA and one
    open valence for the Cys91 SG bond. Its N-terminal N takes the external bond
    back to Ub Arg74 C.
    """
    # NCC(=O)NCC=O : Gly75 N-CA-C(=O) - Gly76 N-CA-C(=O)H  (H stripped below)
    mol = Chem.AddHs(Chem.MolFromSmiles("NCC(=O)NCC=O"))
    AllChem.EmbedMolecule(mol, randomSeed=0xBEEF)
    AllChem.MMFFOptimizeMolecule(mol)
    mol, acyl = strip_acyl_hydrogen(mol)
    # The N-terminal amine also needs an open valence: it receives the peptide
    # bond back to Ub Arg74 C. As built it is a free -NH2 (3 declared bonds), so
    # the external bond would make it 4. Drop one H, leaving -NH- as in a real
    # backbone amide. check_ccd caught this; it is the same class of error as the
    # two bugs in Julian's run, found here instead of on the cluster.
    freeN_idx = [a.GetIdx() for a in mol.GetAtoms()
                 if a.GetSymbol() == "N" and
                 sum(1 for n in a.GetNeighbors() if n.GetSymbol() == "H") == 2]
    if len(freeN_idx) != 1:
        raise ValueError(f"expected one free -NH2 nitrogen, found {len(freeN_idx)}")
    nh = [n.GetIdx() for n in mol.GetAtomWithIdx(freeN_idx[0]).GetNeighbors()
          if n.GetSymbol() == "H"][0]
    em = Chem.RWMol(mol)
    em.RemoveAtom(nh)
    mol = em.GetMol()
    Chem.SanitizeMol(mol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
    acyl = [a.GetIdx() for a in mol.GetAtoms()
            if a.GetSymbol() == "C" and
            sum(1 for n in a.GetNeighbors()
                if n.GetSymbol() == "O" and
                mol.GetBondBetweenAtoms(a.GetIdx(), n.GetIdx()
                                        ).GetBondTypeAsDouble() == 2) == 1 and
            sum(1 for n in a.GetNeighbors() if n.GetSymbol() in ("C", "N")) == 1][0]
    AllChem.EmbedMolecule(mol, randomSeed=0xBEEF)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    # name by PDB convention so the geometry is readable: N/CA/C/O per residue
    conf = mol.GetConformer()
    # walk the backbone from the free N
    freeN = [a.GetIdx() for a in mol.GetAtoms()
             if a.GetSymbol() == "N" and
             sum(1 for n in a.GetNeighbors() if n.GetSymbol() == "C") == 1][0]
    order, seen, cur = [], set(), freeN
    while cur is not None and cur not in seen:
        seen.add(cur); order.append(cur)
        nxt = None
        for n in mol.GetAtomWithIdx(cur).GetNeighbors():
            if n.GetIdx() in seen or n.GetSymbol() == "H":
                continue
            if n.GetSymbol() == "O":
                continue
            nxt = n.GetIdx(); break
        cur = nxt
    labels = ["N1", "CA1", "C1", "N2", "CA2", "C2"]
    names = {}
    for idx, lab in zip(order, labels):
        names[idx] = lab
    # carbonyl oxygens
    for idx, lab in list(names.items()):
        if lab.startswith("C") and not lab.startswith("CA"):
            for n in mol.GetAtomWithIdx(idx).GetNeighbors():
                if n.GetSymbol() == "O":
                    names[n.GetIdx()] = "O" + lab[1:]
    hi = 1
    for a in mol.GetAtoms():
        if a.GetIdx() not in names:
            names[a.GetIdx()] = f"H{hi:02d}"; hi += 1
        a.SetProp("name", names[a.GetIdx()])
    txt = to_cif(mol, "UBGG", GLYGLY_NAME)
    # N1 takes the bond back to Ub Arg74 C; C2 takes the thioester to Cys91 SG
    check_ccd(txt, external={"N1": 1, "C2": 1}, label="UBGG")
    path = os.path.join(outdir, "UBGG_userCCD.cif")
    with open(path, "w") as fh:
        fh.write(txt)
    return path, mol.GetNumAtoms()


def build_glygly_tetrahedral(outdir):
    """Gly-Gly whose C-terminal carbon is sp3, for the TETRAHEDRAL INTERMEDIATE.

    In the real intermediate the amine has added across the C=O: the carbon is
    sp3 with four sigma bonds (CA, O-, S, N) and the oxygen carries the negative
    charge as an alkoxide before it collapses. So the C=O must become a C-O
    SINGLE bond, otherwise the carbon cannot accept both the sulfur and the
    incoming nitrogen without exceeding four.

    This is a chemistry decision, not a formatting workaround: modelling the
    intermediate REQUIRES declaring the carbon sp3. AF3 produced this species by
    accident last time on a carbon that was simply over-valent; here it is
    declared on purpose with the correct hybridisation.
    """
    # C-terminal carbon as an sp3 alcohol carbon: -C(O)(H) with the H stripped,
    # leaving two open valences (for SG and for the incoming amine N).
    mol = Chem.AddHs(Chem.MolFromSmiles("NCC(=O)NCC(O)"))
    AllChem.EmbedMolecule(mol, randomSeed=0xBEEF)
    AllChem.MMFFOptimizeMolecule(mol)
    # open the N-terminal amine for the bond back to Ub Arg74 C
    mol = open_amine_valence(mol)
    # strip BOTH hydrogens from the terminal sp3 carbon so it can take SG and N
    term = None
    for a in mol.GetAtoms():
        if a.GetSymbol() != "C":
            continue
        ohs = [n for n in a.GetNeighbors() if n.GetSymbol() == "O" and
               mol.GetBondBetweenAtoms(a.GetIdx(), n.GetIdx()
                                       ).GetBondTypeAsDouble() == 1]
        heavy = [n for n in a.GetNeighbors() if n.GetSymbol() != "H"]
        if ohs and len(heavy) == 2:
            term = a
            break
    if term is None:
        raise ValueError("no terminal sp3 C-OH carbon found")
    hs = [n.GetIdx() for n in term.GetNeighbors() if n.GetSymbol() == "H"]
    if len(hs) != 2:
        raise ValueError(f"terminal sp3 carbon has {len(hs)} H, expected 2")
    em = Chem.RWMol(mol)
    for h in sorted(hs, reverse=True):
        em.RemoveAtom(h)
    mol = em.GetMol()
    Chem.SanitizeMol(mol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
    AllChem.EmbedMolecule(mol, randomSeed=0xBEEF)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)

    # name by backbone walk from the free N, same convention as build_glygly
    freeN = [a.GetIdx() for a in mol.GetAtoms()
             if a.GetSymbol() == "N" and
             sum(1 for n in a.GetNeighbors() if n.GetSymbol() == "H") == 1][0]
    order, seen, cur = [], set(), freeN
    while cur is not None and cur not in seen:
        seen.add(cur); order.append(cur)
        nxt = None
        for n in mol.GetAtomWithIdx(cur).GetNeighbors():
            if n.GetIdx() in seen or n.GetSymbol() in ("H", "O"):
                continue
            nxt = n.GetIdx(); break
        cur = nxt
    names = {}
    for idx, lab in zip(order, ["N1", "CA1", "C1", "N2", "CA2", "C2"]):
        names[idx] = lab
    for idx, lab in list(names.items()):
        if lab.startswith("C") and not lab.startswith("CA"):
            for n in mol.GetAtomWithIdx(idx).GetNeighbors():
                if n.GetSymbol() == "O":
                    names[n.GetIdx()] = "O" + lab[1:]
    hi = 1
    for a in mol.GetAtoms():
        if a.GetIdx() not in names:
            names[a.GetIdx()] = f"H{hi:02d}"; hi += 1
        a.SetProp("name", names[a.GetIdx()])
    txt = to_cif(mol, "UBGT", GLYGLY_TET_NAME)
    # N1 -> Ub Arg74 C; C2 -> BOTH Cys91 SG and the XisoK amine
    check_ccd(txt, external={"N1": 1, "C2": 2}, label="UBGT")
    path = os.path.join(outdir, "UBGT_userCCD.cif")
    with open(path, "w") as fh:
        fh.write(txt)
    return path, mol.GetNumAtoms()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--only", nargs="*", help="subset of variant codes")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    print("=== XisoK acyl ligands (C01 receives 1 external bond) ===")
    for code, (smi, cip, desc) in VARIANTS.items():
        if a.only and code not in a.only:
            continue
        path, n = build_xisok(code, smi, cip, desc, a.outdir)
        print(f"  {code:9s} {n:2d} atoms  valence OK  -> {os.path.basename(path)}")

    # LisoK with the alpha-amine valence left open, for the product and
    # tetrahedral species where ubiquitin bonds to that nitrogen.
    print("\n=== XisoK with open amine (C01 and N01 each take 1 external bond) ===")
    path, n = build_xisok("LisoK", VARIANTS["LisoK"][0], VARIANTS["LisoK"][1],
                          VARIANTS["LisoK"][2], a.outdir, open_amine=True)
    print(f"  LisoK-openN {n:2d} atoms  valence OK  -> {os.path.basename(path)}")

    print("\n=== Ub Gly-Gly (N1 and C2 each receive 1 external bond) ===")
    path, n = build_glygly(a.outdir)
    print(f"  UBGG      {n:2d} atoms  valence OK  -> {os.path.basename(path)}")

    print("\n=== Ub Gly-Gly, sp3 for the tetrahedral intermediate "
          "(N1 +1, C2 +2 external) ===")
    path, n = build_glygly_tetrahedral(a.outdir)
    print(f"  UBGT      {n:2d} atoms  valence OK  -> {os.path.basename(path)}")
    print(f"\nall files in {a.outdir}/ passed check_ccd with external bonds declared")


if __name__ == "__main__":
    main()
