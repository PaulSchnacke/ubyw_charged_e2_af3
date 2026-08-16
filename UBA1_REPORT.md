# Does Ub(1-75) fail at UBA1? — AF3 results, both catalytic sites

**Status: 7 of 12 variants complete (875 models). The 4 adenylation variants are
re-running after a bug; 1 remains queued.** Everything below is measured from the harvested
models, recomputed here rather than copied from any log. The geometry statistics use all
875; the QC integrity table and the contact map were computed on the first 750 and are
unchanged by the seventh variant.

---

## 1. The answer, in one paragraph

**The UBA1 models do not show a steric block.** Ub(1-75) engages the adenylation site in
every model; it simply sits **2–3.5 Å further out** than Ub(1-76), which is what removing
one glycine from a flexible tail predicts on geometry alone. There is **no single UBA1
residue that admits the long tail and refuses the short one** — what changes is a diffuse
loss of contact across the whole ATP pocket (Fig. 2). And **UBE2W's presence is
chemically irrelevant** — it shifts the C-terminus by 0.03–0.14 Å, roughly **25× less than
tail length does** (Fig. 1c). See §4 for why "irrelevant" and not "not significant": on the
full 125-model set those tiny CIs *exclude* zero, which is a clean demonstration of why
this design reports effect sizes rather than p-values.

Read against the UBE2W half of the project, that localises the defect: if UBA1 accommodates
the shorter donor with only a proportionate loss of reach, the failure is more likely to
sit downstream. **But that comparison is currently unsafe** — the existing UBE2W Ub(1-75)
jobs carry two valence bugs (§5), so their geometry cannot be compared to these numbers
until `jobs_ub75_v2/` runs.

**What this design cannot do is give you a p-value.** With two tail lengths there are two
possible label assignments, so a condition-level permutation test has floor *p* = 0.50 —
no amount of AF3 sampling can beat it. The result is an effect size with a bootstrap CI
and a contact list, and that is the honest ceiling for this experiment.

---

## 2. What was run

| | |
|---|---|
| UBA1 | human, UniProt **P22314**, 1058 aa; Cys632 verified as the annotated glycyl-thioester cysteine |
| sites | **adenylation** (ATP·Mg, non-covalent) and **catalytic Cys632** (non-covalent *and* covalent thioester) |
| tails | Ub(1-76) and Ub(1-75) |
| UBE2W | present and absent |
| sampling | 25 seeds × 5 samples = **125 models per variant** (+1 top-ranked duplicate = 126 files) |
| MSAs | UBA1 16,643 unpaired / 50,000 paired / 4 templates; ubiquitin 7,905–7,964; UBE2W 14,993 |

**Integrity checks, all passing on all 750 models:** Cys632 confirmed a cysteine in
750/750; `acyl_source` resolving to exactly the four expected atom layouts with **no silent
fallbacks** (`UBG1:C1` 250, `protein_Cterm:U75` 250, `UBGG:C2` 125, `protein_Cterm:U76`
125); the declared thioester present in **100%** of covalent models; and **no "Reducing
number of bonds"** in any production log — AF3 kept every bond it was given.

---

## 3. Results

![Reach, Cys632 separation, and thioester planarity](figures/uba1_results.png)

**Fig. 1 — (a)** Distance from the Ub C-terminal carbonyl carbon to the centroid of the
six annotated ATP-binding residues (478, 504, 515, 528, 576, 577). Dots are individual
models, bars are medians. Ub(1-75) sits **+3.43 Å** further out in the non-covalent
treatment (bootstrap 95% CI [+3.41, +3.45]) and **+1.94 Å** in the thioester treatment
(CI [+1.82, +2.05]). Both are consistent with the ~3.5 Å a glycine contributes; the
smaller thioester shift reflects the tail being anchored at Cys632. **(b)** Distance to
Cys632 Sγ in the non-covalent jobs, the honest reach measure (in the thioester jobs this
is 1.5 Å by construction). Every model sits **28–36 Å** away — far outside the reactive
range — for both tail lengths, with and without UBE2W. **(c)** All four measured effects on
a common scale, medians with bootstrap 95% CIs. Tail length moves the C-terminus 1.9–3.4 Å;
UBE2W moves it 0.03–0.14 Å. The shaded band marks the sub-0.5 Å region, below which a shift
carries no chemical meaning for a nucleophilic attack — every UBE2W effect lies inside it.
n = 125 models per variant throughout.

Three things follow.

**The tails differ by exactly the length of a glycine, and nothing more.** A Gly residue
contributes ~3.5 Å of reach, and the measured non-covalent shift is +3.43 Å. There is no
extra penalty beyond the arithmetic — no sign that the shorter tail is excluded, misdocked,
or forced into a different pose.

**Neither tail reaches Cys632, and that is expected.** The two catalytic centres are
**32–36 Å apart** in every model. UBA1's transthiolation requires a large domain rotation
that AF3 does not model, so a static prediction cannot place one Ub C-terminus at both
sites. This is a limitation of the method, **not** evidence that the shorter tail fails —
and it is why the covalent thioester jobs were run: they impose the geometry AF3 will not
find on its own.

**The valence fix holds at production scale.** 358.7° across 250 covalent models, both
tail lengths behaving identically, so the new single-glycine component `UBG1` is as sound
as the established dipeptide `UBGG`.

![Differential contacts between the two tail lengths](figures/uba1_contacts.png)

**Fig. 2 — What the shorter tail loses.** UBA1 residues whose contact frequency differs
by ≥25% between the two tail lengths (contact = any heavy atom within 4.5 Å of the Ub
C-terminal carbonyl carbon; n = 125 models per variant). Asterisked, bold residues are
annotated ATP-binding sites in UniProt P22314. Ub(1-75) loses contact with **477, 478,
479, 480, 515, 574, 577, 598, 626** and **893** — a contiguous sweep of the ATP pocket
rather than one gatekeeper — and gains contact with **624, 625, 631**, which sit beside
Cys632. In other words the shorter tail slides *away* from the adenylation pocket and
*toward* the catalytic cysteine, exactly as a shorter tether would.

**No single residue explains the experimental failure.** If one UBA1 side chain were
rejecting the shorter donor, you would expect one residue contacted at ~100% by Ub(1-76)
and ~0% by Ub(1-75) while its neighbours were unchanged. Instead the loss is spread across
ten residues spanning the pocket, each losing 30–100%. That is the signature of the whole
C-terminus withdrawing by a couple of ångström, not of a steric clash.

---

## 4. Why there is no p-value here, and why that is the right answer

This is worth stating explicitly, because the project's earlier AF3 rounds got it wrong in
the other direction and reported *p* = 0.004 from a design whose floor was 0.5.

**The experimental unit is the construct, not the model.** 875 models sound like a lot,
but they are 7 conditions sampled 125 times each. AF3 seeds of the same complex are
repeated measurements of one system, in the same sense that 125 photographs of one mouse
are not 125 mice. Measured here on all 875: the intraclass correlation is **0.93** — 93% of
the variance sits *between* conditions — giving a design effect of 117, so all 875 models
carry roughly **8 independent observations'** worth of information.

**The design floor.** With 2 tail lengths and one "working" construct there are C(2,1) = 2
distinct label assignments, so the smallest attainable one-sided permutation *p* is
**0.50**. A perfect, flawless separation could not produce a significant result. The
analysis script now prints this before any statistic, so the number can never be quoted
without its ceiling.

**So the deliverable is an effect size with a CI plus the contact list.** The bootstrap
intervals in Fig. 1a are tight ([+3.41, +3.45]) because AF3's within-condition spread is
small, and that tightness is a statement about AF3's reproducibility, not about biology.
Read them as "this is what the model consistently predicts", not "this is significant".

**A worked example of exactly that trap, from this dataset.** When I had 25 models per
variant, UBE2W's effect on reach came out as −0.04 Å with a CI spanning zero, and I wrote
that UBE2W "changes nothing". With the full 125 models the same effect is **+0.03 Å with a
CI of [+0.01, +0.06] — which excludes zero.** Nothing about the biology changed; I just
sampled the same two structures more times, and a bootstrap CI narrows as √n regardless of
whether the underlying quantity means anything.

If this project's convention were to read "CI excludes zero" as a result, that would now
read as *UBE2W significantly repositions the ubiquitin C-terminus* — on a shift of
**0.03 Å**, about 1/30th of a C–C bond and far below the coordinate error of the method
that produced it. The correct reading is the one Fig. 1c is drawn to make unavoidable:
**tail length matters, UBE2W does not, and the ratio is ~25×.** Effect size and a
pre-declared relevance threshold decide that; the interval does not.

This is the same failure as round 2's *p* = 0.004 wearing different clothes. There the
model count inflated a *p*-value; here it narrowed a CI. In both cases the fix is to fix
the *unit* — 125 seeds of one construct are 125 pictures of one system — and to state a
threshold for what size of effect would matter *before* looking at the interval.

**Minima are not used anywhere.** "Closest approach over N models" keeps sliding downward
as sampling deepens — it ranks tail behaviour (arm mobility) rather than positioning. All
statistics here are medians, bootstrap CIs, and fraction-within-cutoff, all of which
converge.

---

## 5. A finding that affects your existing data

**The `jobs_ub75/` UBE2W comparison was run on chemistry known to be wrong**, and it is the
direct comparator for everything above. `jobs_ub75/` was written *after* the valence fix but
reverted to the pre-fix components. Verified with the repo's own checker under the exact
external bonds each job declares:

| bug | effect |
|---|---|
| `LIG-1` retains the spurious `C01–H01` | isopeptide acyl carbon at **5/4** once the bond to the substrate lysine Nζ is added |
| terminal glycines are bare `GLY` (OXT retained) | **N at 4/3 and C at 5/4** — both over-valent; `ube2w_ub76_charged` carries it twice |

The consequence is measured, not inferred: 328° (sp³) with bare `GLY` versus 357–359°
(planar sp²) with the corrected components. Rebuilt jobs are in `jobs_ub75_v2/`, valence-
clean, not yet run.

**Until they run, do not compare the UBE2W Ub(1-75) geometry to the numbers above.** Reach
distances measured far from the reactive centre are less affected, but AF3's distortion
reaches one bond *beyond* the atoms it was asked to bond (§6), so anything near the
isopeptide or thioester carbon is unusable.

---

## 6. Method limitations that bound these conclusions

**AF3 treats declared bonds as connectivity, not geometry.** The thioester C–S comes out at
**1.50 Å** where 1.78–1.81 is correct. More subtly, the *pre-existing* residue is corrupted
too: the bonded cysteine's Cβ–Sγ is ~1.25 Å against ~1.81 for every untouched cysteine in
the same models. Never report a bond length or angle at or adjacent to a declared bond as
physical. The practical cost is nil — tleap rebuilds that junction when it applies real
parameters — but it means the thioester models are useful for *placement*, not for
*geometry*.

**AF3 does not model UBA1's domain rotation.** The 32–36 Å separation between the two
catalytic centres in every model is real for the static structure and uninformative about
the transthiolation step.

**Confidence metrics are QC only.** pLDDT and ipTM have been uninformative for reactivity
across all rounds of this project and are not used as evidence here.

**This is a prediction, not a measurement.** Nothing above establishes that Ub(1-75) *can*
be adenylated — only that the models show no structural reason it could not be. A negative
from a structure predictor is weak evidence; the value of this run is that it makes the
UBE2W hypothesis more attractive by comparison, and it tells you which residues to mutate
if you want to test the UBA1 side experimentally.

---

## 7. What is still running, and what to do next

| | |
|---|---|
| 4 `aden` variants | **resubmitted** (`5aa4b377`, `4238b0ff`) after the graft bug — these carry ATP·Mg and are the direct adenylation-site test |
| batch 4 delivered | `uba1_cys_ub76_with_ube2w`, completing the non-covalent 2×2 — this is the seventh variant |
| batch 4 (`309c38e2`) | still running |
| `jobs_ub75_v2/` | **not submitted** — the corrected UBE2W comparison, needs GPU time |

Read out with:

```bash
for d in ~/ubyw_uba1_models/uba1_*; do echo "$(basename $d) $(ls $d/*.cif|wc -l)/126"; done
python qc_uba1.py ~/ubyw_uba1_models qc_uba1.csv     # on the cluster, stdlib only
python analyse_uba1.py qc_uba1.csv results_uba1/      # locally
```

**The experiment I would prioritise** is `jobs_ub75_v2/`, because it is the comparator that
makes §1's localisation argument valid. The four `aden` variants will refine Fig. 1a but
are unlikely to change its direction — Ub(1-75) already engages the site in every model of
the non-covalent treatment.

**If you want to test the UBA1 side at the bench**, Fig. 2 gives the shortlist: 477–480 and
574/577 lose contact most completely, and 624/625/631 gain it. But note that a diffuse
contact loss is a weak basis for a point mutant — the models are saying *no single residue
is responsible*, which argues for testing the downstream hypothesis instead.
