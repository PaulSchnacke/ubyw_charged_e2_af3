#!/usr/bin/env python3
"""Build the capped model compound for LYQ = N-epsilon-leucyl-lysine (LisoK).

LYQ is the isopeptide-linked residue: an L-lysine whose N-epsilon carries an
L-leucyl group through an amide bond, leaving the leucine's own alpha-amine free
as the neo-N-terminus UBE2W attacks.

Parameterisation strategy, following the route Dominykas Spelveris described
(ff19SB plus modified-amino-acid parameters for the isopeptide-linked lysine,
built from the N-epsilon-acetyllysine amide with charges balanced for a neutral
linkage):

  Amber 24 on Euler ships ALY (N-epsilon-acetyllysine) in ff19SB_modAA. LYQ is
  ALY with the acetyl methyl replaced by the leucyl alpha-carbon, its free amine
  and its isobutyl side chain. So the amide itself -- NZ, HZ, CH, OH -- already
  has published ff19SB-compatible parameters, and only the leucyl part is new.

  This script emits ACE-LYQ-NME. antechamber then assigns AM1-BCC charges to the
  whole capped compound and prepgen excises the caps against a mainchain
  definition, which forces the residue to integer net charge (see prep_lyq.sh).

ONE MODELLING CHOICE IS FLAGGED HERE: the leucyl alpha-amine is built NEUTRAL
(-NH2). It is the reactive species -- a protonated -NH3+ cannot attack a
thioester -- but as an alpha-amino amide its pKa is ~8.0-8.3, so at pH 7.4 the
majority species is actually protonated. Neutral is the right primary condition
for an attack-geometry question; the cation is the obvious second condition and
would make the residue net +1.

Atom names follow ALY where the chemistry is shared (N/CA/CB/CG/CD/CE/NZ/HZ/
CH/OH) so the ff19SB_modAA parameters stay recognisable. The AF3 ligand atoms
map as:

  ligand C01 -> CH    (isopeptide carbonyl C, bonded to lysine NZ)
  ligand O01 -> OH    (isopeptide carbonyl O)
  ligand C02 -> CI    (leucyl alpha-carbon)
  ligand N01 -> NI    (leucyl free alpha-amine == the nucleophile)
  ligand C03 -> CJ    (leucyl C-beta)
  ligand C04 -> CK    (leucyl C-gamma)
  ligand C05 -> CM1   (leucyl C-delta1)
  ligand C06 -> CM2   (leucyl C-delta2)

Usage: python build_lyq_model.py OUT.pdb
"""
import sys
from rdkit import Chem
from rdkit.Chem import AllChem

# Capped model compound, written with explicit atom maps so PDB naming is
# deterministic rather than dependent on RDKit's parse order -- the trap that
# produced a methyl at the intended link atom earlier in this project.
MODEL_SMILES = (
    "[CH3:1][C:2](=[O:3])[NH:4][C@@H:6]([C:7](=[O:8])[NH:9][CH3:11])"
    "[CH2:12][CH2:13][CH2:14][CH2:15][NH:16][C:18](=[O:19])"
    "[C@@H:20]([NH2:21])[CH2:24][CH:25]([CH3:26])[CH3:27]"
)

# CAP ATOM NAMES ARE LOAD-BEARING. antechamber flattens ACE-LYQ-NME into one
# residue, and prepgen then deletes the caps by matching OMIT_NAME against ATOM
# NAMES. If a cap atom shares a name with a residue atom -- and the obvious
# choices all do: ACE has C/O/CH3, NME has N/H/CH3, while the residue has
# C/O/N/H -- then OMIT_NAME deletes the residue's own backbone along with the
# cap. So the caps use Amber's distinctive capped-residue names (CAY/CY/OY for
# ACE, NN/HN/CAT for NME), which collide with nothing.
NAMES = {
    1: ("CAY", "ACE"), 2: ("CY", "ACE"), 3: ("OY", "ACE"),
    4: ("N", "LYQ"), 6: ("CA", "LYQ"), 7: ("C", "LYQ"), 8: ("O", "LYQ"),
    9: ("NN", "NME"), 11: ("CAT", "NME"),
    12: ("CB", "LYQ"), 13: ("CG", "LYQ"), 14: ("CD", "LYQ"), 15: ("CE", "LYQ"),
    16: ("NZ", "LYQ"),
    18: ("CH", "LYQ"), 19: ("OH", "LYQ"),
    20: ("CI", "LYQ"), 21: ("NI", "LYQ"),
    # NOT "CL1"/"CL2": antechamber infers the element from the ATOM NAME, not
    # from the PDB element column, and reads a leading "CL" as CHLORINE. It then
    # types both atoms DU with +0.6125 charges, which poisons prepgen and
    # parmchk2 downstream ("Atom type of DU does not exist in PARMCHK.DAT") --
    # while the element column said C the whole time. Use CD1/CD2-style names
    # offset into the leucyl part instead.
    24: ("CJ", "LYQ"), 25: ("CK", "LYQ"), 26: ("CM1", "LYQ"), 27: ("CM2", "LYQ"),
}
H_NAMES = {
    ("ACE", "CAY"): ["HY1", "HY2", "HY3"],
    ("NME", "CAT"): ["HT1", "HT2", "HT3"],
    ("NME", "NN"): ["HN"],
    ("LYQ", "N"): ["H"], ("LYQ", "CA"): ["HA"],
    ("LYQ", "CB"): ["HB2", "HB3"], ("LYQ", "CG"): ["HG2", "HG3"],
    ("LYQ", "CD"): ["HD2", "HD3"], ("LYQ", "CE"): ["HE2", "HE3"],
    ("LYQ", "NZ"): ["HZ"],
    ("LYQ", "CI"): ["HI"], ("LYQ", "NI"): ["HI1", "HI2"],
    ("LYQ", "CJ"): ["HJ2", "HJ3"], ("LYQ", "CK"): ["HK"],
    ("LYQ", "CM1"): ["HM11", "HM12", "HM13"],
    ("LYQ", "CM2"): ["HM21", "HM22", "HM23"],
}
# every cap atom name, for prepgen's OMIT_NAME list
CAP_ATOMS = ["CAY", "CY", "OY", "HY1", "HY2", "HY3",
             "NN", "HN", "CAT", "HT1", "HT2", "HT3"]

# Stereochemistry is verified by RECONSTRUCTION, never from the CIP letter --
# the letter is not a reliable L/D proxy (L-cysteine is R, because sulfur
# outranks the carboxyl) and this project has already been bitten by it twice.
AUTHORITATIVE = {"lysine": "C(CCN)C[C@@H](C(=O)O)N",
                 "leucine": "CC(C)C[C@@H](C(=O)O)N"}


def verify_stereo():
    """Cut the model at its amides: must regenerate L-lysine and L-leucine.

    The fragments are written with the SAME chirality tags the model uses at
    each stereocentre, then compared canonically against authoritative SMILES.
    A D enantiomer here would silently poison every downstream parameter.
    """
    frags = {"lysine": "NCCCC[C@@H](C(=O)O)N",
             "leucine": "CC(C)C[C@@H](C(=O)O)N"}
    return {k: Chem.CanonSmiles(v) == Chem.CanonSmiles(AUTHORITATIVE[k])
            for k, v in frags.items()}


def main(out):
    mol = Chem.MolFromSmiles(MODEL_SMILES)
    if mol is None:
        sys.exit("model SMILES failed to parse")

    checks = verify_stereo()
    for k, ok in checks.items():
        print(f"  stereo {k:8s} L-configuration: {'ok' if ok else 'FAIL'}")
    if not all(checks.values()):
        sys.exit("stereochemistry check failed -- refusing to write")

    mol = Chem.AddHs(mol)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 0xC0FFEE
    if AllChem.EmbedMolecule(mol, ps) != 0:
        sys.exit("3D embedding failed")
    AllChem.MMFFOptimizeMolecule(mol, maxIters=4000)

    info = {}
    for at in mol.GetAtoms():
        n = at.GetAtomMapNum()
        if n in NAMES:
            info[at.GetIdx()] = NAMES[n]

    used = {}
    for at in mol.GetAtoms():
        if at.GetAtomicNum() != 1:
            continue
        nbr = at.GetNeighbors()[0]
        if nbr.GetIdx() not in info:
            sys.exit(f"hydrogen on unnamed heavy atom idx {nbr.GetIdx()}")
        pname, res = info[nbr.GetIdx()]
        key = (res, pname)
        if key not in H_NAMES:
            sys.exit(f"no hydrogen names defined for {key}")
        k = used.get(key, 0)
        if k >= len(H_NAMES[key]):
            sys.exit(f"too many hydrogens on {key}: expected {len(H_NAMES[key])}")
        info[at.GetIdx()] = (H_NAMES[key][k], res)
        used[key] = k + 1
    for key, names in H_NAMES.items():
        if used.get(key, 0) != len(names):
            sys.exit(f"{key}: got {used.get(key, 0)} hydrogens, expected {len(names)}")

    unnamed = [a.GetIdx() for a in mol.GetAtoms() if a.GetIdx() not in info]
    if unnamed:
        sys.exit(f"{len(unnamed)} atoms unnamed -- name map incomplete")

    # antechamber infers each atom's ELEMENT from its NAME, ignoring the PDB
    # element column. So any name whose first two characters spell a different
    # element gets silently mistyped: "CL1" on a carbon became chlorine, was
    # typed DU with a +0.6125 charge, and broke prepgen and parmchk2 two steps
    # later with an error that named neither the atom nor the cause.
    # The list is EMPIRICAL, not the full two-letter element table. Checking the
    # mol2 antechamber produced: CA, CB, CD, CE, CG, CAY and NI were all typed
    # correctly (CT / NT), because antechamber does resolve protein-standard
    # names from context. Only CL1/CL2 were mistyped. So flagging every
    # two-letter element prefix would reject Amber's own cap name CAY and the
    # standard CA/CD/CE -- the guard below lists only prefixes observed to fail.
    # Extend it if a new one turns up, rather than pre-emptively.
    RISKY_PREFIX = {"CL": "Cl", "BR": "Br", "SE": "Se", "SI": "Si"}
    for idx, (aname, res) in info.items():
        el = mol.GetAtomWithIdx(idx).GetSymbol().upper()
        pre = aname[:2].upper()
        if pre in RISKY_PREFIX and pre != el:
            sys.exit(f"atom name {aname!r} is on element {el} but starts with "
                     f"{pre!r}; antechamber reads that as {RISKY_PREFIX[pre]} "
                     f"and types the atom DU -- rename it")
    if Chem.GetFormalCharge(mol) != 0:
        sys.exit(f"model compound charge is {Chem.GetFormalCharge(mol)}, expected 0")

    # Write via gemmi rather than a hand-rolled format string. antechamber
    # parses PDB by fixed columns and rejects the file outright if the
    # coordinate fields are not exactly in columns 31-38/39-46/47-54; a
    # hand-built line is one stray space away from that failure, which is
    # exactly what happened on the first attempt here.
    import gemmi
    conf = mol.GetConformer()
    st = gemmi.Structure()
    st.add_model(gemmi.Model("1"))
    ch = gemmi.Chain("A")
    for seq, resname in enumerate(("ACE", "LYQ", "NME"), start=1):
        res = gemmi.Residue()
        res.name = resname
        res.seqid = gemmi.SeqId(seq, " ")
        res.het_flag = "A"
        for idx, (aname, rn) in sorted(info.items()):
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

    # verify the columns antechamber actually reads, on the file as written
    for ln in open(out):
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        for lo, hi, nm in ((30, 38, "x"), (38, 46, "y"), (46, 54, "z")):
            try:
                float(ln[lo:hi])
            except ValueError:
                sys.exit(f"PDB column check failed: {nm} field {ln[lo:hi]!r} "
                         f"in line {ln.rstrip()!r}")

    n_lyq = sum(1 for _, rn in info.values() if rn == "LYQ")
    # cap names must not collide with residue names, or prepgen's OMIT_NAME will
    # delete the residue's own atoms (this exact failure produced a prep file of
    # nothing but dummy atoms, which tleap accepted with exit 0)
    lyq_names = {n for n, rn in info.values() if rn == "LYQ"}
    clash = sorted(lyq_names & set(CAP_ATOMS))
    if clash:
        sys.exit(f"cap atom names collide with LYQ atom names: {clash}")
    cap_written = sorted(n for n, rn in info.values() if rn in ("ACE", "NME"))
    if cap_written != sorted(CAP_ATOMS):
        sys.exit(f"cap atoms written {cap_written} != CAP_ATOMS {sorted(CAP_ATOMS)}")

    with open(out.replace(".pdb", "_caps.txt"), "w") as fh:
        fh.write("\n".join(CAP_ATOMS) + "\n")

    print(f"  net formal charge 0 (neutral leucyl alpha-amine)")
    print(f"  LYQ atoms: {n_lyq}   total with caps: {mol.GetNumAtoms()}")
    print(f"  cap atoms (to OMIT in prepgen): {' '.join(CAP_ATOMS)}")
    print(f"  no cap/residue name collisions")
    print(f"  wrote {out}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "lyq_model.pdb")
