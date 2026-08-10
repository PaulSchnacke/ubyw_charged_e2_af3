#!/usr/bin/env python3
"""Valence check for AF3 user-CCD ligands, INCLUDING external bondedAtomPairs.

This is the check that was missing. Both bugs in the first charged run were
over-valent carbons, and every other assertion in the pipeline passed:
stereochemistry, the CCD key set, filename collisions, and which residue each
bond lands on were all verified. Nothing summed the bonds on an atom.

  Bug 1: LisoK C01 declared O01(DOUB) + C02(SING) + H01(SING) = 4, then AF3 added
         the isopeptide bond to Lys NZ -> 5 bonds on a carbon. Result: NZ-C01
         collapsed to 0.9 A and the O-C-N cluster rendered as a cyclopropene.
  Bug 2: the standard GLY CCD keeps OXT, so the thioester carbon had C, O, OXT,
         CA before the Cys91 SG bond -> 5 again. Result: sp3 (angle sum 328 deg)
         where a thioester carbonyl must be planar sp2 (360 deg).

Usage:
    from ccd_valence import check_ccd
    check_ccd(ccd_text, external={"C01": 1})     # raises on violation
"""
import re
import sys

# max total bond order per element, for the neutral species we build here
MAX_VALENCE = {"C": 4, "N": 3, "O": 2, "S": 2, "H": 1, "CL": 1, "BR": 1}
ORDER = {"SING": 1, "DOUB": 2, "TRIP": 3, "AROM": 1.5}


def parse_ccd(text):
    """Return (atoms: name->element, bonds: [(a, b, order_str)]) from a CCD block."""
    atoms, bonds = {}, []
    lines = text.split("\n")
    mode = None
    for ln in lines:
        s = ln.strip()
        if s.startswith("_chem_comp_atom."):
            mode = "atom"
            continue
        if s.startswith("_chem_comp_bond."):
            mode = "bond"
            continue
        if s in ("#", "loop_", "") or s.startswith("_") or s.startswith("data_"):
            continue
        f = s.split()
        # atom rows: <comp_id> <atom_id> <element> <charge> <flag> x y z
        if mode == "atom" and len(f) >= 8 and f[3].lstrip("-+").isdigit():
            atoms[f[1]] = f[2].upper()
        # bond rows: <atom1> <atom2> <order> <aromatic_flag>
        elif mode == "bond" and len(f) == 4 and f[2] in ORDER:
            bonds.append((f[0], f[1], f[2]))
    return atoms, bonds


def valences(text, external=None):
    """Total bond order per atom, counting declared bonds plus external ones.

    `external` maps atom name -> number of external single bonds AF3 will add
    from bondedAtomPairs (e.g. {"C01": 1} for the isopeptide bond to Lys NZ).
    """
    atoms, bonds = parse_ccd(text)
    v = {a: 0.0 for a in atoms}
    for a, b, o in bonds:
        for x in (a, b):
            if x not in v:
                raise ValueError(f"bond references undeclared atom {x!r}")
            v[x] += ORDER[o]
    for a, n in (external or {}).items():
        if a not in v:
            raise ValueError(f"external bond names undeclared atom {a!r}")
        v[a] += n
    return atoms, v


def check_ccd(text, external=None, label="ccd", strict=True):
    """Raise (or return problems) if any atom exceeds its element's valence."""
    atoms, v = valences(text, external)
    problems = []
    for a, tot in sorted(v.items()):
        el = atoms[a]
        lim = MAX_VALENCE.get(el)
        if lim is None:
            continue
        if tot > lim:
            ext = (external or {}).get(a, 0)
            problems.append(
                f"{a} ({el}): total bond order {tot:g} exceeds {lim} "
                f"({tot - ext:g} declared + {ext} external)")
        # an acyl carbon that will receive an external bond must leave room:
        # flag the specific pattern that caused bug 1
        if el == "C" and (external or {}).get(a, 0) and tot == lim:
            pass  # exactly full is correct
    if problems and strict:
        raise ValueError(f"{label}: valence violation\n  " + "\n  ".join(problems))
    return problems


def report(text, external=None, label="ccd"):
    atoms, v = valences(text, external)
    print(f"{label}: {len(atoms)} atoms")
    for a in sorted(v, key=lambda x: (atoms[x] == "H", x)):
        if atoms[a] == "H":
            continue
        ext = (external or {}).get(a, 0)
        lim = MAX_VALENCE.get(atoms[a], "?")
        flag = "  <-- OVER" if isinstance(lim, int) and v[a] > lim else ""
        note = f" (+{ext} external)" if ext else ""
        print(f"  {a:6s} {atoms[a]:2s}  order {v[a]:>4g} / {lim}{note}{flag}")
    return v


if __name__ == "__main__":
    txt = open(sys.argv[1]).read()
    ext = {}
    for kv in sys.argv[2:]:
        k, n = kv.split("=")
        ext[k] = int(n)
    report(txt, ext, label=sys.argv[1])
    p = check_ccd(txt, ext, label=sys.argv[1], strict=False)
    sys.exit(1 if p else 0)
