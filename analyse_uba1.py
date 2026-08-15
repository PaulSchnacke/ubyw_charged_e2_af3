#!/usr/bin/env python3
"""Turn qc_uba1.py output into the answer: does Ub(1-75) fail at UBA1, and if so where?

Run AFTER qc_uba1.py. Needs pandas/numpy/matplotlib (run locally, not on the cluster).

THE STATISTICS, and why they are what they are. This project's earlier AF3 rounds
reported a metric as "anti-predictive" on the strength of an AUC computed with the score
as +distance while the hypothesis was that SHORTER distance means more reactive -- the
sign convention inverted the biology. And the supporting p-values counted AF3 models as
independent replicates when the experimental unit is the site. So:

  * EVERY AUC HERE STATES ITS DIRECTION IN THE OUTPUT. The score is -distance
    ("closer = more reactive"), printed alongside the number. An AUC without a stated
    sign convention is exactly as ambiguous as a bare CIP letter.

  * MODELS ARE NOT REPLICATES. 25 seeds of one construct are 25 pictures of one system.
    The comparison here is 2 tail lengths x 2 UBE2W states x 2 sites = 8 conditions, so
    the honest unit is the CONDITION, and a condition-level test with this design has a
    hard floor on attainable p. That floor is computed and printed, so a small p from
    model-counting can never be mistaken for evidence.

  * NO MINIMA. "Closest approach over N models" is an extreme-value statistic: it keeps
    sliding down as sampling deepens, so it ranks tail behaviour (how mobile the arm is)
    rather than central tendency (how well positioned it is). In the sweep data K21's
    minimum fell from 3.84 A at 5 models to 1.77 A at 100 while K11's barely moved.
    Reported instead: median with a bootstrap CI, and FRACTION WITHIN A CUTOFF, both of
    which converge.

  * AF3 CONFIDENCE IS MODEL QUALITY, NOT REACTIVITY. Reported for QC only, never as
    evidence about the engineering -- it has been uninformative across five rounds.

Usage: python analyse_uba1.py QC.csv OUTDIR
"""
import os
from math import comb
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ATTACK = 4.0          # nucleophilic attack distance on a carbonyl carbon / phosphorus
BOOT = 10000
RNG = np.random.default_rng(0)


def parse_job(j):
    """uba1_<site>_ub<len>_<ube2w> -> (site, tail, with_ube2w)."""
    p = j.split("_")
    return p[1], int(p[2].replace("ub", "")), ("with" in j)


def auc_closer(pos, neg):
    """AUC with score = -distance, i.e. 'closer predicts the positive class'.

    Returned value is P(a random positive is CLOSER than a random negative). 0.5 is
    chance. Ties count 0.5. The DIRECTION IS IN THE NAME so it cannot be misread.
    """
    pos, neg = np.asarray(pos), np.asarray(neg)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    c = sum((1.0 if p < n else 0.5 if p == n else 0.0) for p in pos for n in neg)
    return c / (len(pos) * len(neg))


def boot_ci(x, stat=np.median, n=BOOT):
    x = np.asarray(x)
    if len(x) == 0:
        return (float("nan"),) * 2
    vals = [stat(RNG.choice(x, len(x), replace=True)) for _ in range(n)]
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def frac_within(x, cut=ATTACK):
    x = np.asarray(x)
    return float((x <= cut).mean()) if len(x) else float("nan")


def main(qc_csv, outdir):
    os.makedirs(outdir, exist_ok=True)
    df = pd.read_csv(qc_csv)
    # assign each column explicitly: df[[...]] = df.job.apply(pd.Series) silently
    # misaligns when the returned Series carries its own index labels
    parsed = [parse_job(j) for j in df.job]
    df["site"] = [p[0] for p in parsed]
    df["tail_len"] = [p[1] for p in parsed]
    df["with_ube2w"] = [p[2] for p in parsed]

    # ---- integrity gates before any interpretation -------------------------
    print("=" * 74)
    print("INTEGRITY")
    print(f"  models: {len(df)}   jobs: {df.job.nunique()}")
    bad = df[~df.cys632_is_cys.astype(str).str.lower().isin(["true"])]
    print(f"  Cys632 confirmed a cysteine in {len(df)-len(bad)}/{len(df)} models"
          + (f"   ** {len(bad)} FAILURES **" if len(bad) else ""))
    print("  acyl_source distribution (must match job types, no silent fallbacks):")
    for k, v in df.acyl_source.value_counts().items():
        print(f"     {k:24s} {v}")
    thio = df[df.site == "thio"]
    if len(thio):
        kept = thio.thio_formed.astype(str).str.lower().eq("true").mean()
        print(f"  declared thioester bond present in {kept*100:.0f}% of covalent models")
        print(f"  planarity at acyl C: median {thio.planarity.median():.1f} deg "
              f"(sp2=360, sp3=328.5)")

    # ---- the reach question, per site --------------------------------------
    rows = []
    for site, col in (("aden", "reach_aden"), ("cys", "reach_cys")):
        sub = df[(df.site == site) & df[col].notna()]
        if not len(sub):
            continue
        for e2 in (False, True):
            s = sub[sub.with_ube2w == e2]
            if not len(s):
                continue
            d76 = s[s.tail_len == 76][col].dropna().values
            d75 = s[s.tail_len == 75][col].dropna().values
            if not (len(d76) and len(d75)):
                continue
            lo76, hi76 = boot_ci(d76)
            lo75, hi75 = boot_ci(d75)
            rows.append(dict(
                site=site, ube2w="with" if e2 else "without",
                n76=len(d76), n75=len(d75),
                med76=np.median(d76), ci76=f"[{lo76:.1f},{hi76:.1f}]",
                med75=np.median(d75), ci75=f"[{lo75:.1f},{hi75:.1f}]",
                frac76=frac_within(d76), frac75=frac_within(d75),
                # AUC direction: positive class = Ub(1-76), the one that WORKS
                auc_76_closer=auc_closer(d76, d75)))
    res = pd.DataFrame(rows)

    print("\n" + "=" * 74)
    print(f"REACH: is Ub(1-75) further from the catalytic centre than Ub(1-76)?")
    print(f"  AUC direction: score = -distance, positive class = Ub(1-76).")
    print(f"  So AUC > 0.5 means the WORKING (1-76) construct sits CLOSER, as predicted.")
    print(f"  frac = fraction of models within {ATTACK} A (a converging statistic;")
    print(f"  minima are deliberately NOT used -- they slide with sampling depth).\n")
    if len(res):
        with pd.option_context("display.width", 200):
            print(res.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
        res.to_csv(os.path.join(outdir, "uba1_reach_summary.csv"), index=False)

    # ---- the design floor: what p is even attainable? ----------------------
    n_cond = df.groupby(["site", "tail_len", "with_ube2w"]).ngroups
    print("\n" + "=" * 74)
    print("WHAT THIS DESIGN CAN AND CANNOT SHOW")
    print(f"  {len(df)} models, but the experimental unit is the CONDITION: "
          f"{n_cond} conditions.")
    print("  A p-value computed over models treats 25 seeds of one construct as 25")
    print("  independent observations. They are 25 pictures of one system.")
    # intra-condition clustering -> effective n
    for site, col in (("aden", "reach_aden"), ("cys", "reach_cys")):
        sub = df[(df.site == site) & df[col].notna()]
        if len(sub) < 4:
            continue
        g = sub.groupby(["tail_len", "with_ube2w"])[col]
        msw = g.var(ddof=1).mean()
        msb = g.mean().var(ddof=1)
        m = g.size().mean()
        if msw + msb > 0:
            icc = msb / (msb + msw)
            deff = 1 + (m - 1) * icc
            print(f"  {site}: ICC = {icc:.2f}, design effect = {deff:.1f} -> "
                  f"{len(sub)} models carry ~{len(sub)/deff:.0f} observations' worth")

    # THE DESIGN FLOOR. State this BEFORE quoting any p-value. It is a property of the
    # experiment, not of the data: with k conditions of which p are 'positive', a
    # permutation test at the condition level has only C(k,p) distinct label assignments,
    # so the smallest attainable one-sided p is 1/C(k,p). A perfect, flawless separation
    # cannot beat it. Round 2 of this project reported p = 0.004 from a design whose
    # floor was 0.5, which is what makes this worth printing every time.
    print("\n  DESIGN FLOOR -- the smallest p this comparison could possibly reach:")
    for site in sorted(df.site.unique()):
        sub = df[df.site == site]
        for e2 in (False, True):
            s = sub[sub.with_ube2w == e2]
            k = s.groupby("tail_len").ngroups
            if k < 2:
                continue
            # one tail length is the 'positive' (the working construct)
            floor = 1.0 / comb(k, 1)
            tag = "with UBE2W" if e2 else "no UBE2W"
            print(f"    {site:5s} {tag:11s} {k} tail lengths, 1 positive -> "
                  f"C({k},1) = {comb(k,1)} assignments -> floor p = {floor:.2f}")
    print("    So a condition-level test on 2 tail lengths CANNOT reach significance at")
    print("    all: the comparison is descriptive (effect size and its CI), not inferential.")
    print("    Any small p you see elsewhere for this design came from counting models.")

    # ---- the differential contact map: the mechanistic answer --------------
    print("\n" + "=" * 74)
    print("DIFFERENTIAL CONTACTS: UBA1 residues contacted by the SHORT tail but not")
    print("the long one (>=25% of models), i.e. what the shorter monomer runs into.\n")
    diff_rows = []
    for site in df.site.unique():
        for e2 in (False, True):
            s = df[(df.site == site) & (df.with_ube2w == e2)]
            if not len(s):
                continue
            def freq(tail):
                t = s[s.tail_len == tail]
                if not len(t):
                    return {}
                cnt = {}
                for c in t.contacts.fillna(""):
                    for r in str(c).split(";"):
                        if r:
                            cnt[int(r)] = cnt.get(int(r), 0) + 1
                return {k: v / len(t) for k, v in cnt.items()}
            f76, f75 = freq(76), freq(75)
            if not (f76 and f75):
                continue
            for r in sorted(set(f75) | set(f76)):
                a, b = f75.get(r, 0.0), f76.get(r, 0.0)
                if abs(a - b) >= 0.25:
                    diff_rows.append(dict(site=site,
                                          ube2w="with" if e2 else "without",
                                          uba1_residue=r, freq_ub75=a, freq_ub76=b,
                                          delta=a - b))
    if diff_rows:
        d = pd.DataFrame(diff_rows).sort_values("delta", ascending=False)
        d.to_csv(os.path.join(outdir, "uba1_differential_contacts.csv"), index=False)
        print("  top residues MORE contacted by Ub(1-75):")
        print(d.head(12).to_string(index=False, float_format=lambda x: f"{x:.2f}"))
        print("\n  top residues LESS contacted by Ub(1-75):")
        print(d.tail(8).to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    else:
        print("  none at the 25% threshold.")

    # ---- figure ------------------------------------------------------------
    sites = [s for s in ("aden", "cys", "thio") if (df.site == s).any()]
    fig, axes = plt.subplots(1, max(2, len(sites)), figsize=(4.6 * max(2, len(sites)), 4.3))
    axes = np.atleast_1d(axes)
    for ax, site in zip(axes, sites):
        col = {"aden": "reach_aden", "cys": "reach_cys", "thio": "reach_cys"}[site]
        sub = df[(df.site == site) & df[col].notna()]
        data, labs, cols = [], [], []
        for e2 in (False, True):
            for tail, c in ((76, "#1b6ca8"), (75, "#e0a458")):
                v = sub[(sub.tail_len == tail) & (sub.with_ube2w == e2)][col].dropna().values
                if len(v):
                    data.append(v)
                    labs.append(f"Ub 1-{tail}\n{'+UBE2W' if e2 else 'no E2'}")
                    cols.append(c)
        if not data:
            continue
        bp = ax.boxplot(data, tick_labels=labs, widths=.55, patch_artist=True,
                        showfliers=False, medianprops=dict(color="black", lw=1.5))
        for b, c in zip(bp["boxes"], cols):
            b.set_facecolor(c)
            b.set_alpha(.75)
        for i, v in enumerate(data, 1):
            ax.plot(RNG.normal(i, .06, len(v)), v, ".", ms=2.6, color="#33383f", alpha=.35)
        ax.axhline(ATTACK, color="#2e7d32", lw=1.3, ls="--")
        ax.text(len(data) + .45, ATTACK, f"{ATTACK:g} A", fontsize=8, color="#2e7d32",
                ha="right", va="bottom")
        title = {"aden": "Adenylation site\nUb C-term C -> ATP P-alpha",
                 "cys": "Catalytic Cys632\nUb C-term C -> SG",
                 "thio": "Covalent thioester\ndeclared bond"}[site]
        ax.set_title(title, fontsize=10, loc="left")
        ax.set_ylabel("distance (A)")
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "uba1_reach.png"), dpi=200)
    print(f"\nwrote {outdir}/uba1_reach.png and the CSVs")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
