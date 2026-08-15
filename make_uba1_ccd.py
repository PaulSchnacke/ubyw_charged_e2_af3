#!/usr/bin/env python3
"""Build the SINGLE-glycine CCD (UBG1) needed for a Ub(1-75) thioester.

WHY THIS FILE EXISTS. The Ub(1-76) thioester already has a valence-correct CCD:
UBGG, the Gly75-Gly76 dipeptide with no OXT (make_ccd_v2.build_glygly). Ub(1-75)
needs the same thing one residue shorter -- Gly75 alone, carrying the thioester on
its own carbonyl.

The trap this avoids is documented in results/JULIAN_RUN_ANALYSIS.md as Bug 2:
the standard GLY CCD keeps OXT, so the carbon that becomes the thioester already
has C, O, OXT and CA. Adding the external bond to the catalytic SG makes it
five-coordinate, and AF3 then builds an sp3 carbon (bond-angle sum 328 deg) where
a thioester carbonyl must be planar sp2 (360 deg). It exits 0 and writes a
plausible model.

That bug was fixed for the sweep (jobs_v2/, jobs_sweep/ use UBGG) but jobs_ub75/
was written afterwards with bare GLY ligands again, so the Ub(1-75) comparison
carries the regression. Hence: one custom CCD per tail length, and no bare GLY on
any reactive C-terminus, ever.

  UBG1 atoms   N1 CA1 C1 O1 + 3 H
  open valences  N1 <- Ub Arg74 C (extends the chain)
                 C1 <- catalytic Cys SG (the thioester)

Usage: python make_uba1_ccd.py OUTDIR
"""
import argparse
import os

from rdkit import Chem
from rdkit.Chem import AllChem

from ccd_valence import check_ccd, report
from make_ccd_v2 import strip_acyl_hydrogen, to_cif

UBG1_NAME = "Ub Gly75 with an open thioester valence (no OXT)"


def build_gly1(outdir):
    """One glycine: N-CA-C(=O), no OXT, both terminal valences left open.

    Built the same way as UBGG so the two tail lengths differ only in length:
    present the reactive carbon as an ALDEHYDE in the SMILES, then strip the
    aldehyde hydrogen. Leaving the valence open is what lets AF3's external bond
    complete it at exactly four.
    """
    mol = Chem.AddHs(Chem.MolFromSmiles("NCC=O"))      # N-CA-C(=O)H
    AllChem.EmbedMolecule(mol, randomSeed=0xBEEF)
    AllChem.MMFFOptimizeMolecule(mol)
    mol, _ = strip_acyl_hydrogen(mol)                  # -> C has =O, CA, open

    # The N-terminal amine receives the peptide bond back to Ub Arg74 C. As built
    # it is a free -NH2 with three declared bonds, so the external bond would make
    # nitrogen 4-valent -- the same error class as the two bugs in Julian's run.
    # Drop one H, leaving -NH- as in a real backbone amide.
    freeN = [a.GetIdx() for a in mol.GetAtoms()
             if a.GetSymbol() == "N" and
             sum(1 for n in a.GetNeighbors() if n.GetSymbol() == "H") == 2]
    if len(freeN) != 1:
        raise ValueError(f"expected one free -NH2, found {len(freeN)}")
    nh = [n.GetIdx() for n in mol.GetAtomWithIdx(freeN[0]).GetNeighbors()
          if n.GetSymbol() == "H"][0]
    em = Chem.RWMol(mol)
    em.RemoveAtom(nh)
    mol = em.GetMol()
    Chem.SanitizeMol(mol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
    AllChem.EmbedMolecule(mol, randomSeed=0xBEEF)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)

    # Name explicitly rather than by parse order: antechamber and prepgen both
    # infer meaning from atom NAMES downstream, so a generic C1/N1/UNL set
    # silently produces a garbage residue tree (silent failure #4).
    acyl = [a.GetIdx() for a in mol.GetAtoms()
            if a.GetSymbol() == "C" and
            any(n.GetSymbol() == "O" and
                mol.GetBondBetweenAtoms(a.GetIdx(), n.GetIdx()
                                        ).GetBondTypeAsDouble() == 2
                for n in a.GetNeighbors())][0]
    n_idx = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "N"][0]
    ca_idx = [n.GetIdx() for n in mol.GetAtomWithIdx(n_idx).GetNeighbors()
              if n.GetSymbol() == "C"][0]
    o_idx = [n.GetIdx() for n in mol.GetAtomWithIdx(acyl).GetNeighbors()
             if n.GetSymbol() == "O"][0]
    names = {n_idx: "N1", ca_idx: "CA1", acyl: "C1", o_idx: "O1"}
    hi = 1
    for a in mol.GetAtoms():
        if a.GetIdx() not in names:
            names[a.GetIdx()] = f"H{hi:02d}"
            hi += 1
        a.SetProp("name", names[a.GetIdx()])

    txt = to_cif(mol, "UBG1", UBG1_NAME)
    # THE GATE: N1 takes the bond back to Ub Arg74 C, C1 takes the thioester to
    # the catalytic SG. Declaring both is the whole point -- each half is fine
    # alone and only the sum is over-valent.
    check_ccd(txt, external={"N1": 1, "C1": 1}, label="UBG1")
    path = os.path.join(outdir, "UBG1_userCCD.cif")
    with open(path, "w") as fh:
        fh.write(txt)
    return path, mol.GetNumAtoms(), txt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    path, n, txt = build_gly1(a.outdir)
    print(f"UBG1: {n} atoms -> {path}")
    report(txt, external={"N1": 1, "C1": 1}, label="UBG1")


if __name__ == "__main__":
    main()
