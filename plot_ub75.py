#!/usr/bin/env python3
"""Figure: does UBE2W refuse a ubiquitin one residue shorter?

Usage: python plot_ub75.py CONTACT_FREQ_JSON OUT.png

Panel a  thioester geometry per variant. The charged jobs had the bond DECLARED, so
         "formed" there is compliance, not evidence -- the non-covalent bar is the
         only one where AF3 chose the distance freely.
Panel b  UBE2W residues whose contact frequency changes, split by WHICH ubiquitin
         residue makes the contact. Contacts made only by Gly76 vanish trivially in
         the shorter variant; only contacts made by residues present in BOTH
         (Leu73) are informative.
"""
import json
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt

WT, SHORT, FREE = "#0b5394", "#e08214", "#8c8c8c"
ALARM = "#c1121f"
ATTACK = 4.0
BOND = 2.2

# Which ubiquitin residue makes each contact, from the donor probe over 30 models per
# variant. This split is the whole point of panel b: without it a 100% -> 0% drop
# looks like a structural finding when it is the absence of the contacting atom.
DONOR = {"TYR85": "gly76", "PRO82": "gly76", "TRP144": "gly76",
         "ILE94": "leu73", "SER118": "leu73"}

GEOM = [("Ub(1–76)\ncharged", 1.59, 1.73, 126, WT,    "bond imposed"),
        ("Ub(1–75)\ncharged", 1.60, 1.68, 126, SHORT, "bond imposed"),
        ("Ub(1–76)\nnon-covalent", 3.23, 3.48, 0, FREE, "free")]


def main(freq_json, out, n=126):
    freq = json.load(open(freq_json))
    A, B = "ube2w_ub76_charged", "ube2w_ub75_charged"

    rows = []
    for x in sorted(set(freq[A]) | set(freq[B])):
        fa, fb = freq[A].get(x, 0) / n * 100, freq[B].get(x, 0) / n * 100
        if abs(fb - fa) >= 15:
            rows.append((x, fa, fb, DONOR.get(x, "mixed")))
    rows.sort(key=lambda r: r[1] - r[2], reverse=True)

    fig = plt.figure(figsize=(7.2, 3.1))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.5], wspace=0.42)
    ax, ax2 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

    # ---- panel a: closest approach of the acyl carbon to Cys91 SG ---------------
    xs = range(len(GEOM))
    for i, (lab, tmin, tmed, formed, col, kind) in enumerate(GEOM):
        # lollipop: a single summary per variant, so a dot with a stem rather than a
        # bar (§6.2) -- and the value axis is a distance, where a bar's length would
        # imply a ratio to an arbitrary floor.
        ax.plot([i, i], [0, tmed], color=col, lw=1.0, alpha=0.5, zorder=1)
        ax.plot(i, tmed, "o", color=col, ms=7, zorder=3,
                mec="white", mew=0.6)
        ax.plot(i, tmin, "v", color=col, ms=4.5, zorder=3, alpha=0.85)
        ax.annotate(f"{tmed:.2f}", (i, tmed), textcoords="offset points",
                    xytext=(9, 1), fontsize=6, color=col, va="center")
    ax.axhline(BOND, color=ALARM, lw=0.8, ls="--", zorder=0)
    ax.text(-0.42, BOND, "bonded 2.2 Å", fontsize=6, color=ALARM,
            va="bottom", ha="left")
    ax.set_xticks(list(xs))
    ax.set_xticklabels([g[0] for g in GEOM], fontsize=6)
    ax.set_ylabel("acyl C to Cys91 SG (Å)", fontsize=7)
    ax.set_ylim(0, 4.3)
    ax.set_title("Shorter donor bonds just as readily —\nbut only when told to",
                 loc="left")
    # Mark which bars had the bond imposed: without this the panel overstates.
    for i, g in enumerate(GEOM):
        if g[5] == "bond imposed":
            ax.annotate("imposed", (i, 0.12), ha="center", fontsize=5.5,
                        color=META_GREY)
    ax.margins(x=0.14)

    # ---- panel b: contact frequency change, split by contacting residue --------
    y = range(len(rows))
    for i, (res, fa, fb, donor) in enumerate(rows):
        ax2.plot([fa, fb], [i, i], color=META_GREY, lw=0.8, zorder=1)
        ax2.plot(fa, i, "o", color=WT, ms=5, zorder=3, mec="white", mew=0.5)
        ax2.plot(fb, i, "o", color=SHORT, ms=5, zorder=3, mec="white", mew=0.5)
    ax2.set_yticks(list(y))
    labels = []
    for res, fa, fb, donor in rows:
        mark = " *" if donor == "gly76" else ""
        labels.append(f"{res.capitalize()}{mark}")
    ax2.set_yticklabels(labels, fontsize=6)
    ax2.invert_yaxis()
    ax2.set_xlabel("models with contact (%)", labelpad=3, fontsize=7)
    ax2.set_xlim(-6, 112)
    ax2.set_title("Most lost contacts are made by Gly76 itself;\nonly Ile94 and Ser118 reflect a real shift",
                  loc="left")
    from matplotlib.lines import Line2D
    ax2.legend(handles=[Line2D([], [], marker="o", ls="", color=WT, ms=5,
                               label="Ub(1–76)"),
                        Line2D([], [], marker="o", ls="", color=SHORT, ms=5,
                               label="Ub(1–75)")],
               frameon=False, fontsize=6, loc="lower left", handletextpad=0.3,
               bbox_to_anchor=(0.02, 0.10))
    ax2.margins(y=0.05)

    for a, L in ((ax, "a"), (ax2, "b")):
        panel_letter(a, L)

    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"wrote {out}")
    print(f"  panel b rows: {len(rows)}; "
          f"{sum(1 for r in rows if r[3] == 'gly76')} are Gly76-only (trivial)")
    return fig


if __name__ == "__main__":
    main(*sys.argv[1:3])
