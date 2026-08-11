#!/usr/bin/env python3
"""Figure: the charged (thioester) state under a force field.

Usage: python plot_charged_md.py SERIES_JSON OUT.png

Panel a  the attack coordinate, against the apo run that had to be stopped. This is
         the comparison that matters: with ubiquitin absent the pose relaxed to
         7.4 +- 0.3 A and NEVER came inside 4 A; with the thioester present it
         tightens instead.
Panel b  thioester bond length and the acyl carbon's planarity -- the chemistry AF3
         could not hold (it pyramidalised the same carbon to 338.6 deg).
Panel c  interface native contacts. Included because backbone RMSD cannot tell
         "came off" from "rearranged while staying bound", and because the honest
         reading of this panel is more cautious than panels a and b.
"""
import json
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

MD, APO, ALARM = "#0b5394", "#8c8c8c", "#c1121f"
UB_COL, SUB_COL = "#e08214", "#0b5394"
ATTACK = 4.0
FRAME_NS = 0.05          # ntwx=25000 x dt=0.002 ps, verified against the mdin
APO_MEAN, APO_SD = 7.4, 0.3   # apo K11 final quarter, from the stopped runs


def main(series_json, out):
    S = json.load(open(series_json))
    n = len(S["attack"])
    t = np.arange(n) * FRAME_NS

    fig = plt.figure(figsize=(7.2, 4.9))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.30, 0.52, 0.52],
                          hspace=0.75, wspace=0.34)
    ax = fig.add_subplot(gs[0, :])
    ax2 = fig.add_subplot(gs[1, 0])                 # thioester bond length
    ax2b = fig.add_subplot(gs[2, 0], sharex=ax2)    # planarity, its own strip
    ax3 = fig.add_subplot(gs[1:, 1])                # contacts, full height

    # ---- panel a: the attack coordinate vs the apo baseline --------------------
    a = np.array(S["attack"])
    ax.axhspan(APO_MEAN - APO_SD, APO_MEAN + APO_SD, color=APO, alpha=0.22, zorder=0)
    ax.axhline(APO_MEAN, color=APO, lw=1.0, ls="-", zorder=1)
    ax.annotate("apo (no ubiquitin): 7.4 ± 0.3 Å, 0/1200 frames ≤4 Å",
                (t[-1] * 0.99, APO_MEAN), fontsize=6, color=APO,
                ha="right", va="bottom")
    ax.axhline(ATTACK, color=ALARM, lw=0.9, ls="--", zorder=1)
    ax.annotate("attack distance 4 Å", (t[1], ATTACK), fontsize=6, color=ALARM,
                ha="left", va="bottom")
    ax.plot(t, a, color=MD, lw=0.7, alpha=0.85, zorder=3)
    # running mean, so the trend is visible without smoothing away the spread
    w = 21
    if n > w:
        rm = np.convolve(a, np.ones(w) / w, mode="valid")
        ax.plot(t[w // 2: w // 2 + len(rm)], rm, color=MD, lw=1.8, zorder=4)
    frac = (a <= ATTACK).mean() * 100
    q = a[int(n * 0.75):]
    ax.set_ylabel("XisoK amine to acyl C (Å)")
    ax.set_xlabel("production time (ns)")
    ax.set_ylim(0, max(8.4, a.max() * 1.05))
    ax.set_title(f"With ubiquitin bonded, the nucleophile closes on the acyl carbon: "
                 f"{q.mean():.2f} Å at the end, {frac:.0f}% of frames ≤4 Å",
                 loc="left")
    ax.margins(x=0.01)

    # ---- panel b: the chemistry AF3 could not hold -----------------------------
    th = np.array(S["thio"])
    pl = np.array(S["planar"])
    ax2.plot(t, th, color=MD, lw=0.7)
    ax2.axhline(1.8104, color=META_GREY, lw=0.8, ls=":")
    ax2.set_ylabel("C–S (Å)", fontsize=7)
    ax2.set_ylim(1.70, 1.95)
    ax2.set_yticks([1.75, 1.85])
    ax2.tick_params(labelbottom=False)
    # Short enough for one line: the AF3 contrast moves to panel a's caption text
    # rather than competing with panel a's x-axis for space.
    ax2.set_title(f"Thioester stable: {th.mean():.2f} Å, {pl.mean():.0f}° planar",
                  loc="left")

    ax2b.plot(t, pl, color=MD, lw=0.7)
    ax2b.axhline(123.0, color=META_GREY, lw=0.8, ls=":")
    ax2b.set_ylabel("O=C–S (°)", fontsize=7)
    ax2b.set_xlabel("production time (ns)", fontsize=7)
    ax2b.set_ylim(112, 136)
    ax2b.set_yticks([115, 125, 135])
    # One shared reference note for both strips, in the lower strip's whitespace.
    ax2b.annotate("dotted: ideal sp$^2$ values", (t[-1], 113.5), fontsize=5.5,
                  color=META_GREY, ha="right", va="bottom")

    # ---- panel c: did the interface hold? --------------------------------------
    cs = np.array(S["c_subenz"], dtype=float)
    cu = np.array(S["c_ubenz"], dtype=float)
    ax3.plot(t, cs / cs[0] * 100, color=SUB_COL, lw=1.0, label="substrate–enzyme")
    ax3.plot(t, cu / cu[0] * 100, color=UB_COL, lw=1.0, label="ubiquitin–enzyme")
    ax3.set_ylabel("native contacts retained (%)", fontsize=7)
    ax3.set_xlabel("production time (ns)", fontsize=7)
    ax3.set_ylim(0, 105)
    ax3.legend(frameon=False, fontsize=6, loc="lower left", handlelength=1.2,
               borderaxespad=0.4)
    ax3.set_title("Interfaces settle from the AF3 pose", loc="left")

    for a_, L in ((ax, "a"), (ax2, "b"), (ax3, "c")):
        panel_letter(a_, L)

    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"wrote {out}")
    print(f"  {n} frames = {n * FRAME_NS:.1f} ns")
    print(f"  attack: first {a[0]:.2f}  final-quarter {q.mean():.2f} A  "
          f"{(a <= ATTACK).sum()}/{n} frames <=4 A")
    print(f"  thioester {th.mean():.3f} +- {th.std():.3f} A; planarity {pl.mean():.1f} deg")
    print(f"  contacts retained at end: substrate {cs[-1]/cs[0]*100:.0f}%, "
          f"ubiquitin {cu[-1]/cu[0]*100:.0f}%")
    return fig


if __name__ == "__main__":
    main(*sys.argv[1:3])
