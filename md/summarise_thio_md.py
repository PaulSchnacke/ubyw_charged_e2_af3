#!/usr/bin/env python3
"""Summarise a charged-state trajectory analysis directory.

Usage: python summarise_thio_md.py ANALYSIS_DIR PREFIX

Reads <PREFIX>_thio.dat, _attack.dat, _planar.dat, _rmsd_sub.dat, _rmsd_ub.dat,
_contacts_subenz.dat, _contacts_ubenz.dat as written by cpptraj and prints one line
per observable.

Exists as a FILE rather than an embedded heredoc: when this was inlined in an
f-string job command, the outer f-string evaluated the inner `{...}` expressions in
the submitting kernel instead of passing them through. A staged file cannot be
interpolated by accident.

Stdlib only -- the cluster's bare python3 has neither numpy nor gemmi.
"""
import os
import statistics as st
import sys

FRAME_NS = 0.05          # ntwx=25000 x dt=0.002 ps, verified against the mdin
ATTACK = 4.0


def col(path, c=1):
    out = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if len(f) > c:
            try:
                out.append(float(f[c]))
            except ValueError:
                pass
    return out


def main(d, prefix):
    series = [("thioester C-S", "thio", "A"),
              ("attack NI->acyl", "attack", "A"),
              ("O=C-S angle", "planar", "deg"),
              ("substrate CA RMSD", "rmsd_sub", "A"),
              ("ubiquitin CA RMSD", "rmsd_ub", "A")]
    for label, key, unit in series:
        path = os.path.join(d, f"{prefix}_{key}.dat")
        if not os.path.exists(path):
            print(f"  {label:20s} MISSING {os.path.basename(path)}")
            continue
        v = col(path)
        if not v:
            print(f"  {label:20s} NO DATA")
            continue
        q = v[int(len(v) * 0.75):]
        extra = ""
        if key == "attack":
            n_in = sum(1 for x in v if x <= ATTACK)
            extra = f"   <={ATTACK:.0f}A in {n_in}/{len(v)} ({n_in / len(v) * 100:.0f}%)"
        print(f"  {label:20s} n={len(v):4d} ({len(v) * FRAME_NS:5.1f} ns)  "
              f"first {v[0]:7.2f}  mean {st.mean(v):7.2f}  "
              f"final-quarter {st.mean(q):7.2f} {unit}{extra}")

    for label, key in (("substrate-enzyme", "contacts_subenz"),
                       ("ubiquitin-enzyme", "contacts_ubenz")):
        path = os.path.join(d, f"{prefix}_{key}.dat")
        if not os.path.exists(path):
            continue
        v = col(path)
        if v:
            q = v[int(len(v) * 0.75):]
            print(f"  {label:20s} native contacts: first {v[0]:.0f}  "
                  f"mean {st.mean(v):.0f}  final-quarter {st.mean(q):.0f}  "
                  f"({st.mean(q) / v[0] * 100:.0f}% retained)")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
