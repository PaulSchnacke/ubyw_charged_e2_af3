#!/usr/bin/env python3
"""Write the handful of C-S terms a protein-protein thioester needs, taken from GAFF2.

THE ACTUAL PROBLEM, stated precisely. Modelling a thioester is solved -- GAFF2 has
every term. Bonding two protein chains in tleap is also solved -- an explicit
`bond` command does it. What is NOT covered is the intersection: a bond between
two chains whose atoms carry *protein* atom types, because `parm19.dat` defines no
`C-S` bond at all. No standard amino acid has a thioester, so the protein force
field never needed one.

This is why the lysine-linked ubiquitin-chain case is easier than ours. An
isopeptide bond needs `C-N` amide terms, and those are standard protein
parameters, present in parm19 already. Our thioester needs `C-S`, which is not.
Same topology, one extra ingredient.

So rather than deriving anything, map GAFF2's thioester values onto the protein
types they correspond to:

    GAFF2 c  (sp2 carbonyl carbon)  -> protein C   (backbone/acyl carbonyl)
    GAFF2 ss (thioether sulfur)     -> protein S   (Met/Cys sulfur)
    GAFF2 c3 (sp3 carbon)           -> protein CT
    GAFF2 o  (carbonyl oxygen)      -> protein O

That mapping is chemically like-for-like, and it is a transfer of published
bonded terms, not a derivation. Values are extracted from gaff2.dat at runtime so
there is no transcription error, and each term is checked against parm19.dat first
so only genuinely missing ones are written.

Usage: python make_thioester_link_frcmod.py AMBERHOME OUT.frcmod
"""
import os
import re
import sys

# GAFF2 term -> (protein term, kind). Order matters only for the report.
# ff19SB does NOT use CT for every sp3 carbon: the test showed tleap asking for
# XC-C-S, because ff19SB gives the alpha carbon type XC. Any sp3 carbon that can
# neighbour the linkage needs its own angle, so emit the GAFF2 c3 value for each.
SP3 = ("CT", "XC", "2C", "3C")
WANT = [("c -ss", "C -S", "BOND"),
        ("o -c -ss", "O -C -S", "ANGLE")]
WANT += [("c3-c -ss", f"{x:<2s}-C -S", "ANGLE") for x in SP3]
WANT += [("c -ss-c3", f"C -S -{x:<2s}", "ANGLE") for x in SP3]
WANT += [("X -c -ss-X", "X -C -S -X", "DIHE")]


def read_terms(path, keys):
    """{key: rest-of-line} for lines whose leading term matches a key exactly."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path, errors="ignore") as fh:
        for ln in fh:
            s = ln.rstrip("\n")
            for k in keys:
                if s.startswith(k) and k not in out:
                    # guard against 'c -ss' also matching 'c -ss-c3'
                    tail = s[len(k):]
                    if tail[:1] in ("", " ", "\t"):
                        out[k] = s
    return out


def numbers(line, term):
    return [float(x) for x in re.findall(r"-?\d+\.?\d*", line[len(term):])]


def main(amberhome, out):
    gaff2 = os.path.join(amberhome, "dat/leap/parm/gaff2.dat")
    parm19 = os.path.join(amberhome, "dat/leap/parm/parm19.dat")
    if not os.path.exists(gaff2):
        sys.exit(f"no gaff2.dat at {gaff2}")

    src = read_terms(gaff2, sorted({g for g, _, _ in WANT}))
    have = read_terms(parm19, [p for _, p, _ in WANT])

    missing = [w for w in WANT if w[1] not in have]
    print(f"gaff2.dat: found {len(src)}/{len({g for g,_,_ in WANT})} distinct source terms")
    for g, p, kind in WANT:
        status = "already in parm19" if p in have else "MISSING from parm19"
        got = "yes" if g in src else "NOT FOUND in gaff2"
        print(f"  {kind:6s} {g:12s} -> {p:12s}  gaff2={got:18s} {status}")

    absent = [g for g, _, _ in WANT if g not in src]
    if absent:
        sys.exit(f"cannot proceed: {absent} not found in gaff2.dat")
    if not missing:
        print("\nNothing to write: parm19 already covers every term.")
        return

    # frcmod format: ONE title line, then sections. Amber rejects the entire file
    # if any line is malformed, so keep the title to a single line and write the
    # canonical field layout for each section.
    # frcmod structure, verified against files tleap accepts (our lyq.frcmod and
    # the shipped frcmod.ff19SB): ONE title line, then MASS, BOND, ANGLE, DIHE,
    # IMPROPER, NONBON -- each keyword on its own line, each section terminated by a
    # blank line. MASS must be present even when empty; omitting it made Amber
    # reject the ENTIRE file with "Could not load parameter set", which looks like a
    # missing-parameter problem rather than a formatting one.
    lines = ["Thioester C-S terms transferred from GAFF2 for a protein-protein acyl link",
             "MASS", ""]
    for kind in ("BOND", "ANGLE", "DIHE"):
        rows = [(g, p) for g, p, k in missing if k == kind]
        lines.append(kind)
        if not rows:
            lines.append("")
            continue
        for g, p in rows:
            n = numbers(src[g], g)
            if kind == "BOND":
                lines.append(f"{p}  {n[0]:8.2f}  {n[1]:8.4f}")
            elif kind == "ANGLE":
                lines.append(f"{p}  {n[0]:8.3f}  {n[1]:9.3f}")
            else:
                # GAFF2 dihedral line: IDIVF, barrier/2, phase, periodicity
                lines.append(f"{p}  {int(n[0]):2d}  {n[1]:8.3f}  {n[2]:9.3f}  "
                             f"{n[3]:5.1f}")
        lines.append("")
    lines += ["IMPROPER", "", "NONBON", ""]

    with open(out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwrote {out} with {len(missing)} term(s):")
    for g, p, kind in missing:
        print(f"  {kind:6s} {p}   <- gaff2 {g}: {src[g][len(g):].strip()[:40]}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
