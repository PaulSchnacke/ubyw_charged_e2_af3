#!/usr/bin/env python3
"""Figure: the corrected-chemistry sweep. Reach still does not predict reactivity.

Usage: python plot_sweep.py QC_CSV OUT.png
"""
import sys

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

WORKS, FAILS, UNTESTED, ALARM = "#0b5394", "#404040", "#9e9e9e", "#c1121f"
ATTACK = 4.0
KNOWN = {11: "works", 21: "fails"}


def load(csv):
    Q = pd.read_csv(csv)
    Q["job_base"] = Q.job.str.replace(r"_\d{8}_\d{6}$", "", regex=True)
    Q = Q[Q.model.str.contains(r"seed-\d+_sample-\d+")].copy()
    Q["sample_id"] = Q.model.str.extract(r"(seed-\d+_sample-\d+)")[0]
    return Q.drop_duplicates(subset=["job_base", "sample_id"])


def main(csv, out):
    try:
        from kernel import apply_figure_style, panel_letter
        apply_figure_style(sizes=(8, 7, 6))
    except Exception:
        mpl.rcParams.update({"font.size": 8, "axes.titlesize": 8,
                             "axes.labelsize": 8, "xtick.labelsize": 6,
                             "ytick.labelsize": 6, "legend.fontsize": 7})

    Q = load(csv)
    sites = Q[Q.job_base.str.contains("_xisok_")].copy()
    sites["site"] = sites.job_base.str.extract(r"_k(\d+)_")[0].astype(int)
    spec = Q[Q.job_base.str.startswith("sumo2_k11_")].copy()
    spec["species"] = spec.job_base.str.replace("sumo2_k11_", "").str.replace("_ube2w", "")

    fig = plt.figure(figsize=(9.6, 3.2))
    gs = fig.add_gridspec(1, 3, wspace=0.42, left=0.07, right=0.985,
                          top=0.84, bottom=0.20)

    # ---- panel a: all 8 lysines, ranked by closest approach --------------------
    ax = fig.add_subplot(gs[0, 0])
    order = sorted(sites.site.unique(),
                   key=lambda s: sites[sites.site == s].n01_catcys_sg.min())
    rng = np.random.default_rng(0)
    for i, s in enumerate(order):
        v = sites[sites.site == s].n01_catcys_sg.dropna().values
        col = WORKS if KNOWN.get(s) == "works" else (FAILS if KNOWN.get(s) == "fails"
                                                     else UNTESTED)
        ax.scatter(np.full(len(v), i) + rng.normal(0, 0.10, len(v)), v,
                   s=3, color=col, alpha=0.45, lw=0)
        ax.hlines(v.min(), i - 0.34, i + 0.34, color=col, lw=1.6, zorder=3)
    ax.axhline(ATTACK, color=ALARM, ls="--", lw=0.9, zorder=1)
    # place the threshold label in the upper-left whitespace, not on the line:
    # every site column has points near 4 A, so a label at the line overlaps data
    ax.text(0.02, 0.995, "– – attack distance 4 Å", transform=ax.transAxes,
            fontsize=6, color=ALARM, ha="left", va="top")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"K{s}" for s in order], fontsize=6)
    for t, s in zip(ax.get_xticklabels(), order):
        t.set_color(WORKS if KNOWN.get(s) == "works" else "#444444")
    ax.set_ylabel("XisoK amine to catalytic Cys91 (Å)")
    rank = order.index(11) + 1
    ax.set_title(f"The one working site ranks {rank} of {len(order)}", loc="left")
    ax.scatter([], [], s=14, color=WORKS, lw=0, label="UbyW works (K11)")
    ax.scatter([], [], s=14, color=FAILS, lw=0, label="UbyW fails (K21)")
    ax.scatter([], [], s=14, color=UNTESTED, lw=0, label="untested")
    ax.legend(frameon=False, fontsize=6, loc="upper right", ncol=1,
              columnspacing=0.9, handletextpad=0.3, borderaxespad=0.3)
    ax.margins(y=0.10)
    ax.spines[["top", "right"]].set_visible(False)

    # ---- panel b: the modification drives docking (control vs modified) -------
    ax2 = fig.add_subplot(gs[0, 1])
    lc = spec[spec.species == "lyscontrol"].nearest_lys_nz_d.dropna().values
    xi = spec[spec.species == "xisok"].n01_catcys_sg.dropna().values
    for i, (v, col, lab) in enumerate([(lc, UNTESTED, "unmodified Lys"),
                                       (xi, WORKS, "XisoK")]):
        ax2.scatter(np.full(len(v), i) + rng.normal(0, 0.09, len(v)), v,
                    s=3, color=col, alpha=0.45, lw=0)
        ax2.hlines(np.median(v), i - 0.30, i + 0.30, color=col, lw=1.8, zorder=3)
    ax2.axhline(ATTACK, color=ALARM, ls="--", lw=0.9)
    ax2.set_xticks([0, 1])
    ax2.set_xticklabels(["unmodified\nLys", "XisoK"], fontsize=6)
    ax2.set_ylabel("nucleophile to Cys91 (Å)")
    ax2.set_title(f"Modification drives docking\n"
                  f"median {np.median(lc):.0f} Å \u2192 {np.median(xi):.1f} Å", loc="left")
    ax2.margins(x=0.30, y=0.08)
    ax2.spines[["top", "right"]].set_visible(False)

    # ---- panel c: thioester hybridisation across species ---------------------
    ax3 = fig.add_subplot(gs[0, 2])
    got = [(s, spec[spec.species == s].thio_planarity.dropna().values)
           for s in ("charged", "tetrahedral")]
    got = [(s, v) for s, v in got if len(v)]
    for i, (s, v) in enumerate(got):
        col = WORKS if s == "charged" else UNTESTED
        ax3.scatter(np.full(len(v), i) + rng.normal(0, 0.08, len(v)), v,
                    s=3, color=col, alpha=0.5, lw=0)
        ax3.hlines(np.median(v), i - 0.30, i + 0.30, color=col, lw=1.8, zorder=3)
    ax3.axhline(360, color="#666666", ls=":", lw=0.9)
    ax3.text(0.02, 360, "planar sp$^2$", fontsize=6, color="#666666",
             ha="left", va="bottom", transform=ax3.get_yaxis_transform())
    ax3.axhline(328.5, color=ALARM, ls="--", lw=0.9)
    ax3.text(0.02, 328.5, "tetrahedral sp$^3$", fontsize=6, color=ALARM,
             ha="left", va="top", transform=ax3.get_yaxis_transform())
    ax3.set_xticks(range(len(got)))
    ax3.set_xticklabels([s for s, _ in got], fontsize=6)
    ax3.set_ylabel("bond-angle sum at thioester C (°)")
    ax3.set_title("Charged state pyramidalises anyway", loc="left")
    ax3.margins(x=0.34, y=0.16)
    ax3.spines[["top", "right"]].set_visible(False)

    try:
        from kernel import panel_letter
        for a, L in ((ax, "a"), (ax2, "b"), (ax3, "c")):
            panel_letter(a, L)
    except Exception:
        pass

    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"wrote {out}")
    print(f"  K11 ranks {rank}/{len(order)} on closest approach")
    return fig


if __name__ == "__main__":
    main(*sys.argv[1:3])
