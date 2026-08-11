#!/usr/bin/env python3
"""Refuse a residue whose GUESSED force-field terms sit at the reactive centre.

Why this exists. LYT (thioester) and LYX (tetrahedral intermediate) both passed
every check we had: antechamber assigned AM1-BCC charges, prepgen produced a
valid prep with integer charge, and tleap built a peptide from them and reported
`Total unperturbed charge: 0.000000`. All three signals said "usable".

They are not usable. parmchk2 flagged 6 terms for LYT and 3 for LYX as
`ATTN, need revision`, and set every one of them to a force constant of **zero**:

    LYT:  C -S bond, C -S -CT angle, O -C -S angle, CT-C -S angle,
          O -C -S -CT torsion, CT-C -S -CT torsion
    LYX:  OH-CT-S, NT-CT-S, NT-CT-OH angles

Those are not peripheral terms. They are the thioester bond itself and the angles
that define the tetrahedral centre. With zero force constants there is no
restoring force: the C-S bond has no equilibrium length, rotation about it is
free, and the sp3 centre can flatten or invert. A trajectory would run happily and
report numbers that describe nothing.

So a passing tleap is not evidence of a usable residue. This script is the missing
check: it fails when any zero-force-constant term involves an atom of the reactive
centre, so such a residue cannot silently reach production MD.

Usage: python check_frcmod_attn.py RESIDUE.frcmod --reactive S C  [--max-attn 0]
"""
import argparse
import re
import sys


def parse(frcmod):
    """[(section, term, tokens, is_attn)] for every parameter line."""
    rows, section = [], None
    for ln in open(frcmod):
        s = ln.rstrip("\n")
        if re.match(r"^(MASS|BOND|ANGLE|DIHE|IMPROPER|NONBON)", s.strip()):
            section = s.strip().split()[0]
            continue
        if not s.strip() or s.strip().startswith("#"):
            continue
        if section is None:
            continue
        attn = "ATTN" in s
        term = s[:11].strip()
        nums = re.findall(r"-?\d+\.\d+", s.split("ATTN")[0])
        rows.append((section, term, [float(x) for x in nums], attn))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("frcmod")
    ap.add_argument("--reactive", nargs="*", default=[],
                    help="atom TYPES defining the reactive centre, e.g. S C")
    ap.add_argument("--max-attn", type=int, default=None,
                    help="fail if more than this many ATTN terms (any location)")
    a = ap.parse_args()

    rows = parse(a.frcmod)
    attn = [r for r in rows if r[3]]
    reactive = {t.strip().upper() for t in a.reactive}

    zero_at_centre = []
    for section, term, nums, is_attn in attn:
        types = {t.strip().upper() for t in term.split("-") if t.strip()}
        # the first numeric field of a BOND/ANGLE/DIHE line is its force constant
        k = nums[0] if nums else 0.0
        if k == 0.0 and (not reactive or types & reactive):
            zero_at_centre.append((section, term, k))

    print(f"{a.frcmod}: {len(rows)} parameter lines, {len(attn)} flagged ATTN")
    for section, term, nums, _ in attn:
        k = nums[0] if nums else float("nan")
        mark = "  <-- ZERO force constant" if k == 0.0 else ""
        print(f"    [{section:8s}] {term:14s} k={k:8.3f}{mark}")

    fail = False
    if zero_at_centre:
        print(f"\nFAIL: {len(zero_at_centre)} guessed term(s) with ZERO force "
              f"constant involve the reactive centre "
              f"({', '.join(sorted(reactive)) or 'any atom'}):")
        for section, term, k in zero_at_centre:
            print(f"    {section:8s} {term}")
        print("  These terms have no restoring force, so the geometry they define "
              "is unconstrained.\n  This residue is NOT usable for production MD "
              "without QM-derived parameters.")
        fail = True
    if a.max_attn is not None and len(attn) > a.max_attn:
        print(f"\nFAIL: {len(attn)} ATTN terms exceeds --max-attn {a.max_attn}")
        fail = True
    if not fail:
        print("\nPASS: no zero-force-constant guessed terms at the reactive centre")
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
