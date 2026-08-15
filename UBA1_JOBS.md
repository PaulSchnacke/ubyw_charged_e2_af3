# UBA1 with Ub(1-76) vs Ub(1-75): the two catalytic sites

**Status: submitted to Euler 2026-08-15.** Stub validation `f2ea3a80` (GPU, ~30 min),
MSA search `4e124915` (CPU, 12-20 h). Production inference is NOT yet submitted — it
waits on both. See [RECOVERY.md](RECOVERY.md) to resume after any interruption.

## The question

Ub(1-75) was engineered to be accepted in place of Ub(1-76); experimentally it failed.
The failure can sit at **UBA1** (charging) or at **UBE2W** (accepting the shorter
donor). The UBE2W half is `jobs_ub75/` + `analyse_ub75.py`. This is the UBA1 half.

UBA1 has two chemically distinct catalytic centres and they are modelled differently
**on purpose**:

| site | chemistry | how modelled here |
|---|---|---|
| **adenylation** (ATP site, residues 478/504/515/528/576-577) | Ub C-terminal carboxyl + ATP → **Ub-AMP**, a mixed anhydride (acyl phosphate) | **non-covalent**: Ub + ATP + Mg, no bond declared |
| **catalytic cysteine (Cys632)** | Ub Gly76 **thioester** on Cys632 Sγ (UniProt: "Glycyl thioester intermediate") | **both** non-covalent (reach) *and* covalent (ligand trick) |

Why the asymmetry: this project's thioester machinery is validated (GAFF2 terms, the
split-ubiquitin topology, a valence gate), whereas **no acyl-phosphate CCD exists here
yet**, and an acyl phosphate is the chemistry AF3 handles worst — it has no thioester
chemistry either (acyl carbon pyramidalises to 338.6°) and it does not refine bond
lengths (isopeptide C–N came out at 1.07 Å where an amide is 1.33 Å). Declaring a
covalent Ub–AMP bond would hand AF3 a five-coordinate-prone phosphorus and a carboxyl
it wants to keep OXT on. The adenylate is therefore a **reach** measurement until MD,
where the parameters actually matter.

## The 12 jobs

`jobs_uba1/`, 25 seeds each. 2 sites × 2 tail lengths × UBE2W present/absent, plus the
4 covalent thioester variants.

| job | chains | res | ligands | bonds |
|---|---|---|---|---|
| `uba1_aden_ub76_no_ube2w` | A,U | 1134 | ATP, MG | – |
| `uba1_aden_ub76_with_ube2w` | A,U,B | 1285 | ATP, MG | – |
| `uba1_aden_ub75_no_ube2w` | A,U | 1133 | ATP, MG | – |
| `uba1_aden_ub75_with_ube2w` | A,U,B | 1284 | ATP, MG | – |
| `uba1_cys_ub76_no_ube2w` | A,U | 1134 | – | – |
| `uba1_cys_ub76_with_ube2w` | A,U,B | 1285 | – | – |
| `uba1_cys_ub75_no_ube2w` | A,U | 1133 | – | – |
| `uba1_cys_ub75_with_ube2w` | A,U,B | 1284 | – | – |
| `uba1_thio_ub76_no_ube2w` | A,U(1-74) | 1132 | **UBGG** | 2 |
| `uba1_thio_ub76_with_ube2w` | A,U(1-74),B | 1283 | **UBGG** | 2 |
| `uba1_thio_ub75_no_ube2w` | A,U(1-74) | 1132 | **UBG1** | 2 |
| `uba1_thio_ub75_with_ube2w` | A,U(1-74),B | 1283 | **UBG1** | 2 |

Chains: `A` = UBA1 (UniProt P22314, 1058 aa, Cys632 verified), `U` = ubiquitin,
`B` = UBE2W (Q96B02, 151 aa, Cys91), `G` = the glycine-tail ligand.

Covalent bonds, both thioester variants:
```
U:74:C   -> G:1:N1     extend the backbone into the ligand
G:1:C2   -> A:632:SG   the thioester   (UBGG, Ub 1-76)
G:1:C1   -> A:632:SG   the thioester   (UBG1, Ub 1-75)
```

## A regression this fixes: `jobs_ub75/` uses bare GLY

`results/JULIAN_RUN_ANALYSIS.md` Bug 2: the standard `GLY` CCD **retains OXT**, so the
carbon that becomes the thioester already has C, O, OXT and CA. AF3's external bond to
SG makes it five-coordinate, and the result is an sp3 carbon (angle sum 328°) where a
thioester carbonyl must be planar sp2 (360°). Exit 0, plausible model, wrong chemistry.

That was fixed for the sweep — `jobs_v2/` and `jobs_sweep/` use the custom `UBGG`
(no OXT) — but **`jobs_ub75/` was written afterwards and reverted to bare `GLY`**:

```
jobs_v2/sumo2_k11_charged_ube2w.json    ligands=[LIG-1, UBGG]   <- correct
jobs_sweep/sumo2_k11_charged_ube2w.json ligands=[LIG-1, UBGG]   <- correct
jobs_ub75/ube2w_ub75_charged.json       ligands=[LIG-1, GLY]    <- Bug 2 again
jobs_ub75/ube2w_ub76_charged.json       ligands=[LIG-1, GLY, GLY]  <- Bug 2 again
```

So the Ub(1-75) **UBE2W** comparison, if it was run from those JSONs, carries the same
sp3-thioester artefact. That needs a rerun with `UBGG`/`UBG1` before its numbers mean
anything. Recorded here rather than silently fixed, per the project's convention.

The rule going forward: **never use bare `GLY` on a reactive C-terminus.** One custom
CCD per tail length —

* `UBGG` (Gly75-Gly76, from `make_ccd_v2.build_glygly`) — thioester carbon `C2`
* `UBG1` (Gly75 alone, from `make_uba1_ccd.py`, **new**) — thioester carbon `C1`

Both verified with the external bonds declared:
```
UBG1: C1 order 4/4 (+1 external)   N1 order 3/3 (+1 external)   no OXT
UBGG: C2 order 4/4 (+1 external)   N1 order 3/3 (+1 external)   no OXT
```

## The gate in `build_uba1_jobs.py`

No job is written unless it passes. Each check maps to a failure that has actually
happened in this project:

1. **Valence** of every custom ligand *including the external bonds AF3 will add* —
   the pentavalent-carbon class. Each half looks fine alone; only the sum is wrong.
2. **Every bond has a ligand partner.** AF3 3.0.1 silently discards polymer–polymer
   bonds (exit 0, atoms 26 Å apart), so a protein→protein thioester would vanish.
3. **The SG target really is a cysteine** at that residue number in that chain's
   sequence — catches an off-by-one from numbering, which has bitten this project
   (the N-terminal-proline offset).
4. **Chain ids in bonds exist** in the job.

It fired for real during the build: the first run refused `uba1_thio_ub76_*` because
`UBGG_userCCD.cif` was not in `ccd_uba1/`, and wrote nothing for those two rather than
emitting jobs with a missing CCD.

## Why a stub run first

`stubs_uba1/` (via `make_stubs.py`): single-sequence MSAs, 1 seed. Inference runs in
~10 min of GPU instead of 4-6 h of database search, and it settles the only things that
must be true before production compute:

* did AF3 **keep** both declared bonds (`grep -i "Reducing number of bonds"`) — the one
  line separating a real model from a meaningless one;
* the thioester **C–S length and the bond-angle sum at the acyl carbon** — is the carbon
  sp2 (~360°) as it must be, or has it collapsed toward sp3 (328°)?
* **how AF3 renumbered the chains.** It absorbs a leading single-residue ligand into the
  protein chain (in the earlier round chain `U` came back as 75 residues, not 74, with
  the second glycine as chain `UA`), so QC must resolve chains from the `covale`
  records, never from the input naming.
* that **ATP + Mg** is accepted alongside a 1058-residue chain.

A stub says nothing about *placement* — single-sequence AF3 has no coevolution signal.
It is a chemistry and topology check only.

## MSA strategy

MSAs depend **only on protein sequence**, so three donor searches cover all 12 jobs
(verified programmatically: 0 uncovered protein chains). `msa_donors/` holds
UBA1 + Ub(76|75|74) + UBE2W; `graft_msa.py` copies alignments onto each production job
**by sequence, not by chain id**, because chain letters differ between job layouts and
matching by id would silently attach the wrong alignment.

Ubiquitin and UBE2W alignments already exist in `~/ubyw_msa_archive` and
`~/ubyw_charged_msa`. **UBA1 is new** and is the 4-6 h item.

Two things the MSA job does deliberately:

* writes to **`$HOME/ubyw_uba1_msa/`**, not scratch (scratch is purged after 15 days),
  and copies each donor's result to `$HOME` **the moment that donor finishes**, so a
  kill or requeue cannot lose a completed search;
* **resumes** — a donor whose `_data.json` is already present and non-empty is skipped,
  so re-running the same job costs nothing for work already done;
* returns only small **text summaries**, not the augmented JSONs. Those are ~34 MB
  each, and a 34 MB augmented JSON previously exceeded the transfer threshold and made
  a job that had *succeeded* report failure. They stay on Euler and are read from
  `$HOME` by the inference step.

## Sizing

These are 4× the residue count of the earlier charged jobs (321 res → up to 1285), so
the earlier `--gres=gpumem:24g` guidance does not carry over: the stub requests
**40 GB**. Confirm the stub's actual GPU memory use before sizing the production runs.

## What to measure (QC, not yet written)

For the adenylation jobs, the reach coordinate is **Ub C-terminal carbonyl carbon → ATP
Pα**, since that is the bond that forms. For the cysteine jobs it is **Ub C-terminal
carbonyl carbon → Cys632 SG**. For the covalent thioester jobs the bond is declared, so
the questions are whether AF3 kept it and what geometry it built, plus the
**differential contact list** — UBA1 residues contacted by the short tail but not the
long one at ≥25% model frequency, which is the direct answer to "are there particular
residues that do not permit the shorter monomer".

Report **ensembles, not top-ranked models**, and prefer **fraction within a cutoff** or
a low quantile with a bootstrap CI over "closest approach": a minimum is an
extreme-value statistic that keeps sliding as more models are sampled, and in the sweep
data K21's minimum fell from 3.84 Å (5 models) to 1.77 Å (100 models) while K11's barely
moved — the ranking was largely a ranking of arm mobility. Site count, not model count,
is the replicate.

AF3 confidence (pLDDT/ipTM) is **model quality only**, not evidence about the
engineering: it has been uninformative for reactivity across five rounds.
