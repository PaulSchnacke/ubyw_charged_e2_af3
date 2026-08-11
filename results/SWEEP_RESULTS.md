# Corrected-chemistry sweep: 12 jobs, 1200 models

All ran overnight on Euler with the valence-corrected ligands. Every declared
bond formed in **all 1200 models**, no `Reducing number of bonds` line anywhere.

Note on counting: AF3 writes 100 seed×sample models plus a top-ranked duplicate
at the job root (101 files), and it appends `_YYYYmmdd_HHMMSS` when an output
directory already exists, so six jobs appear twice. Analysis uses the 100
`seed-N_sample-M` models per job after de-duplicating on that key.

## Result 1: reach still does not predict reactivity — now on all 8 lysines

Ranked by closest approach across 100 models per site:

| site | outcome | closest | median | ≤4 Å | interface |
|---|---|---|---|---|---|
| K21 | **fails** | **1.77** | 4.84 | 15/100 | 15.6 |
| K45 | untested | 1.87 | 5.68 | 9/100 | 11.7 |
| K33 | untested | 1.89 | 5.15 | 21/100 | 13.8 |
| K35 | untested | 1.97 | 4.59 | 36/100 | 14.3 |
| K42 | untested | 2.59 | 4.70 | 26/100 | 14.3 |
| K7 | untested | 2.75 | 4.96 | 2/100 | 13.6 |
| **K11** | **WORKS** | 3.78 | 4.52 | 7/100 | 9.2 |
| K5 | untested | 4.18 | 4.96 | 0/100 | 11.7 |

**K11, the only site known to work, ranks 7th of 8.** The known-failing K21 ranks
first. On the two sites with known outcomes the sample-level AUC is **0.391**
(p = 0.004) — significantly *worse* than chance, i.e. the metric is
anti-predictive rather than merely uninformative.

This is the fourth consistent negative, and the best-controlled: same protein,
same construct, same enzyme, **same MSA**, 100 models per site, and the six
untested sites were scored blind.

### The one metric that "works" is noise

Median distance ranks K11 first — by 0.07 Å over K35. Bootstrapping the median
(2000 resamples): K11 95% CI [4.43, 4.63], K35 [4.25, 5.07],
**P(K11 < K35) = 0.58**. A coin flip. It is also the metric round 2 rejected, for
the same reason: the median of a flexible arm reports where the arm idles, not
where it can reach.

Interface size runs the wrong way too — K11 has the *smallest* interface (9.2
residues vs K21's 15.6), so it ranks 8th of 8 there.

## Result 2: the modification drives docking — replicated at n=100

| | nucleophile to Cys91 |
|---|---|
| unmodified Lys (control) | median **24.86 Å** (min 4.15) |
| XisoK at the same site | median **4.52 Å** (min 3.78) |

p = 2×10⁻³³. In the control the nearest lysine to the active site is usually a
*different* one (K46 in 60/100 models, K22 in 22/100) — the unmodified K11 is not
even the closest candidate. The modification is what recruits the site.

This is the project's most robust positive finding, and it is a **recognition**
signal, not a reactivity predictor: it holds regardless of whether the site
actually reacts.

## Result 3: the charged state pyramidalises on its own

This was the open question from last night's stub run, and the answer inverted.

| species | thioester C–S | bond-angle sum at C | amine → thioester C |
|---|---|---|---|
| charged | 1.71 Å | **338.6 ± 1.7°** | median **3.12 Å**, 97/100 ≤4 Å |
| tetrahedral (declared) | 1.80 Å | 328.9 ± 1.1° | 1.40–1.55 Å |

The stub run (one seed) gave 356.8° — nearly planar. Across 100 models the
charged state sits at **338.6°**, between planar sp² (360°) and tetrahedral sp³
(328.5°), only ~10° from the job where the bond is *declared*. And the free amine
is 2.46–3.12 Å from that carbon in 97/100 models, with no bond specified.

Our input geometry was correct and unambiguous — the ideal O–C–CA angle is
124.9° in `UBGG` (sp²) versus 110.1° in `UBGT` (sp³). AF3 received the right
hybridisation and partially collapsed it.

**A hypothesis I tested and could not support:** that the approaching amine causes
the pyramidalisation. Spearman(amine distance, planarity) = **+0.049, p = 0.63**
across the 100 charged models — no relationship. So the pyramidalisation is not
demonstrably driven by the amine; it may simply be AF3's prior for that local
environment. Stated as a limitation, not a mechanism.

What this means practically: **the "charged" and "tetrahedral" jobs are not cleanly
separable species in AF3's hands.** Asking for a pre-transfer state with a free
amine and getting something 10° from a tetrahedral intermediate means the
distinction we wanted to draw is not one AF3 maintains. A force field would be
needed to hold the thioester planar — the MD pipeline already has a validated
isopeptide residue and could be extended.

## Result 4: AF3 does not refine bond lengths

`NZ–C01` is 0.93–1.09 Å across every job (ideal amide C–N 1.33 Å), and this is
not our CCD: our ideal carbonyl is 1.221 Å and AF3 wrote 1.343 Å. The same
looseness appeared in rounds 2–5 (1.053 Å mean over 350 models), before the
valence bug existed. AF3 3.0.1 treats `bondedAtomPairs` as *connectivity*,
satisfied loosely.

Usable for **where things are**; not for bond-level geometry.

## Recommendation

Stop pursuing attack geometry as a UbyW predictor. Four rounds, four negatives,
and this one is anti-predictive with 100 models per site and blind untested sites.
The consistent picture across all of them: AF3 models **reachability**, and
reachability is not reactivity.

Worth doing instead:
* **MD on the charged state** with a force field that holds the thioester planar,
  which is the only way to ask this structurally without AF3 improvising the
  chemistry.
* **The six untested sites are now predictions.** K35 (36/100 ≤4 Å) and K33
  (21/100) score high, K5 (0/100) lowest. If any get assayed, that is a real test
  of whether the metric means anything — worth recording now rather than after.
