#!/usr/bin/env python3
"""Check LYQ's amide-half charges against ALY's published ff19SB_modAA values.

Why this check exists. LYQ's charges come from AM1-BCC on a capped model
compound -- a reasonable but unvalidated route. However the reactive half of
LYQ (the isopeptide amide: NZ, HZ, CH, OH) is chemically identical to ALY
(N-epsilon-acetyllysine), which Amber 24 ships with published, peer-reviewed
ff19SB-compatible charges. So that half HAS a reference and should be compared
against it rather than trusted.

A large disagreement here would mean AM1-BCC has put the amide dipole somewhere
different from the published parameterisation -- which matters, because the
amide dipole is what orients the arm and the amine is the atom we measure.

Also checks the thing that actually breaks simulations: integer net charge.

Usage: python check_lyq_charges.py lyq.prep [ALY_REFERENCE.lib]
"""
import sys
import re

# ALY charges from $AMBER_EULER_ROOT/dat/leap/lib/mod_amino19.lib (Amber 24).
# Quoted here so the check runs without the Amber tree; verified by reading that
# file on Euler.
ALY_REF = {"NZ": -0.703093, "HZ": 0.378922, "CH": 0.857267, "OH": -0.609159,
           "CE": -0.213091, "CD": 0.065973, "CG": 0.155883, "CB": -0.232522,
           "N": -0.415700, "H": 0.271900, "C": 0.597200, "O": -0.567900}
# atoms where LYQ and ALY are the SAME chemistry -> charges should be close
AMIDE_ATOMS = ["NZ", "HZ", "CH", "OH"]
BACKBONE_ATOMS = ["N", "H", "C", "O"]
TOL_AMIDE = 0.25      # e; AM1-BCC vs RESP will not agree exactly
TOL_TOTAL = 0.001     # e; net charge must be integer to this precision


def parse_prep(path):
    """Atom name -> charge from an Amber prep (prepi) file."""
    q = {}
    for ln in open(path):
        f = ln.split()
        # prep atom lines: index name type topo_type ... charge
        if len(f) >= 11 and f[0].isdigit() and re.match(r"^[A-Za-z]", f[1]):
            try:
                q[f[1]] = float(f[-1])
            except ValueError:
                continue
    return q


EXPECTED_HEAVY = {"N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ",
                  "CH", "OH", "CI", "NI", "CJ", "CK", "CM1", "CM2"}


def main(prep):
    q = parse_prep(prep)
    if not q:
        sys.exit(f"no charges parsed from {prep}")

    # FIRST CHECK, and the one that catches the failure mode that fooled tleap:
    # a wrong prepgen mainchain definition deletes the whole residue and emits a
    # prep file of nothing but DUMM atoms. tleap then builds a peptide with exit
    # code 0, noting only "LYQ: no atoms" among its output. Count real atoms.
    real = {k: v for k, v in q.items() if not k.startswith("DUMM")}
    missing = EXPECTED_HEAVY - set(real)
    if missing:
        sys.exit(f"LYQ has {len(real)} real atoms; MISSING heavy atoms "
                 f"{sorted(missing)}.\nThis is almost always a bad prepgen "
                 f"mainchain definition (MAIN_CHAIN/OMIT_NAME) deleting the "
                 f"residue -- check lyq.mc and prepgen.log for "
                 f"'Number of omited atoms'.")
    q = real
    total = sum(q.values())
    print(f"LYQ: {len(q)} atoms (all {len(EXPECTED_HEAVY)} expected heavy atoms "
          f"present), net charge {total:+.6f} e")

    fails = []
    if abs(total - round(total)) > TOL_TOTAL:
        fails.append(f"net charge {total:+.6f} is not integer "
                     f"(off by {total - round(total):+.6f} e) -- this accumulates "
                     f"as a system-wide neutrality error")
    if round(total) != 0:
        print(f"  NOTE net charge rounds to {round(total):+d}. Expected 0 for the "
              f"neutral alpha-amine form; +1 would mean the cation was built.")

    print(f"\namide half vs ALY (published ff19SB_modAA, tolerance {TOL_AMIDE} e):")
    for a in AMIDE_ATOMS:
        if a not in q:
            fails.append(f"amide atom {a} absent from LYQ")
            continue
        d = q[a] - ALY_REF[a]
        flag = "ok" if abs(d) <= TOL_AMIDE else "DEVIATES"
        print(f"  {a:4s} LYQ {q[a]:+.4f}   ALY {ALY_REF[a]:+.4f}   "
              f"delta {d:+.4f}  {flag}")
        if abs(d) > TOL_AMIDE:
            fails.append(f"{a} charge deviates from ALY by {d:+.4f} e")

    print(f"\nbackbone vs ALY (informational -- prepgen redistributes cap charge "
          f"here, so deviation is expected):")
    for a in BACKBONE_ATOMS:
        if a in q:
            print(f"  {a:4s} LYQ {q[a]:+.4f}   ALY {ALY_REF[a]:+.4f}   "
                  f"delta {q[a] - ALY_REF[a]:+.4f}")

    # the nucleophile must be a negatively polarised amine nitrogen, or the
    # electrostatics of the very interaction we are studying are wrong
    if "NI" in q:
        print(f"\nnucleophile NI charge {q['NI']:+.4f} e")
        if q["NI"] > -0.3:
            fails.append(f"nucleophile NI is only {q['NI']:+.4f} e -- too "
                         f"positive for a neutral amine nitrogen")
    else:
        fails.append("nucleophile atom NI absent")

    if fails:
        print("\nFAILED:")
        for f in fails:
            print(f"  - {f}")
        sys.exit(1)
    print("\nPASS: net charge integer, amide half consistent with ALY, "
          "nucleophile correctly polarised")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "lyq.prep")
