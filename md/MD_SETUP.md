# Test MD: SUMO2(K11LisoK) + uncharged UBE2W

A feel-for-it run, not the real experiment. Purpose: find out whether MD on this
system is tractable on Euler, what it costs, and whether the AF3 pose survives
contact with a force field.

## Why Amber rather than GROMACS

Dominykas Spelveris's protocol builds topologies in Amber and converts to
GROMACS with ParmEd. On Euler that conversion is the wrong trade:

| | available |
|---|---|
| `gromacs/2021.4` | only version; needs `stack/2024-06 gcc/12.2.0 openmpi/4.1.6`; CPU/MPI build |
| `amber/24` | `pmemd.cuda` and `pmemd.cuda.MPI`, CUDA 12.1, plus tleap/antechamber/parmed/cpptraj/MMPBSA.py |

So Amber gives GPU throughput here and GROMACS does not, and staying in Amber
end-to-end removes a ParmEd round-trip that Dominykas explicitly had to verify
did not break the custom isopeptide linkage. Everything else follows his setup:
ff19SB, OPC water, 1.5 nm padding, 0.15 M NaCl, staged minimisation, restrained
NVT, stepwise NPT with restraints released, 310 K / 1 bar, 2 fs, PME, 1.0 nm
cutoffs, independent velocities per replicate.

Two deliberate substitutions, both because `pmemd` lacks the GROMACS algorithms:
Langevin thermostat instead of velocity-rescale, Monte Carlo barostat instead of
Parrinello-Rahman. Both are the standard well-tested `pmemd.cuda` choices.

## The isopeptide residue: LYQ

Dominykas built his CYL-LYI linkage "from the corresponding ff19SB cysteine and
Nε-acetyllysine amide parameters, with partial charges adjusted to retain a
neutral linkage". The same logic applies here, and Amber 24 makes it easier than
expected: **ALY (Nε-acetyllysine) ships in `ff19SB_modAA`**, so the amide half of
our residue has published ff19SB-compatible parameters:

```
ALY  NZ  N   q=-0.703093     CH  C  q= 0.857267
     HZ  H   q= 0.378922     OH  O  q=-0.609159
```

LYQ (Nε-leucyl-lysine) is ALY with the acetyl methyl replaced by the leucyl
α-carbon, its free α-amine, and its isobutyl side chain. Route:

1. `build_lyq_model.py` writes the capped model compound ACE-LYQ-NME with
   deterministic PDB atom names, and verifies by reconstruction that both
   stereocentres are L. **The CIP letter is not checked** — it is not a reliable
   L/D proxy (L-cysteine is *R*), and this project has been bitten by that twice.
2. `antechamber -c bcc -at amber` assigns AM1-BCC charges to the capped compound.
3. `prepgen` excises the ACE/NME caps against a mainchain definition, which is
   what forces **integer** residue charge. This is preferable to hand-adjusting
   charges: a hand-balanced "neutral linkage" is where a residual −0.03 e per
   copy hides.
4. `parmchk2` fills missing bonded terms; the count flagged `ATTN` is reported.
5. tleap builds `ACE-ALA-LYQ-ALA-NME` as a test, proving LYQ splices into a chain
   from both sides before any solvation is attempted.

`check_lyq_charges.py` then compares the amide-half charges against ALY's
published values — the reactive part of this residue has a reference, so it gets
checked rather than trusted.

## One modelling choice worth arguing about

The leucyl α-amine is built **neutral** (−NH2). It is the reactive species — a
protonated −NH3+ cannot attack a thioester — but as an α-amino amide its pKa is
~8.0–8.3, so at pH 7.4 the *majority* species is protonated. Neutral is the
right primary condition for an attack-geometry question; the cation is the
obvious second condition and would make LYQ net +1. Worth running both if this
goes beyond a test.

## What "uncharged UBE2W" means here

No ubiquitin attached at Cys91 — both cysteines are free thiols (CYS, not
CYX/CYM), the E2 resting state. This is deliberately the permissive case:
rounds 2 and 3 showed a *bare* catalytic cysteine accepts an amine from
reactive and unreactive sites alike. So this run cannot test reactivity. It
tests whether the pose is physical, which is a prerequisite for the
loaded-thioester experiment.

## Observables

`analyse_md.py` measures three things AF3 cannot answer:

1. **NI–Cys91 SG distance over time.** AF3 gave 4.6 Å median across 300
   independent single-shot guesses. MD asks whether that is a minimum the system
   *stays* in.
2. **Interface stability.** Round 5 showed the modification drives the docking;
   if the interface dissolves here, that docking is a prediction artefact.
3. **UBE2W C-terminal tail RMSF.** The specific thing worth looking at: AF3
   pLDDT falls from 91 (residues 1–140) to **41 at Cys151**, i.e. AF3 is
   signalling it does not know where the tail is. MD gives an amplitude, and
   whether the tail touches the substrate or the isopeptide arm.

## Scale

3 × 20 ns (60 ns total) rather than Dominykas's 3 × 100 ns. Enough to see whether
the pose holds and what the tail does; **not** enough for converged MM/GBSA, so
no interaction energies are reported from this run.

## Files

```
build_lyq_model.py    capped ACE-LYQ-NME model compound, stereochemistry verified
prep_lyq.sh           antechamber -> prepgen -> parmchk2 -> tleap test peptide
check_lyq_charges.py  compare the amide half against ALY's published charges
prepare_system.py     fuse AF3 ligand + lysine into one LYQ residue
run_md.sh             solvate, minimise, equilibrate, 3 x 20 ns production
analyse_md.py         distance trace, interface, tail RMSF
```

## LYQ parameterisation: result

Built and validated on Euler (Amber 24). `lyq.prep` + `lyq.frcmod` are the
usable residue.

```
LYQ: 40 atoms (all 17 expected heavy atoms present), net charge +0.000002 e

amide half vs ALY (published ff19SB_modAA, tolerance 0.25 e):
  NZ   LYQ -0.5593   ALY -0.7031   delta +0.1437  ok
  HZ   LYQ +0.3064   ALY +0.3789   delta -0.0725  ok
  CH   LYQ +0.6109   ALY +0.8573   delta -0.2463  ok
  OH   LYQ -0.6122   ALY -0.6092   delta -0.0031  ok

nucleophile NI charge -0.8782 e
frcmod: 112 lines, 0 parameters flagged ATTN (none guessed)
```

Net charge is integer to 2e-6, the amide dipole agrees with the published ALY
parameterisation within AM1-BCC-vs-RESP tolerance, and the nucleophile carries a
strongly negative charge as a neutral amine nitrogen should. `parmchk2` flagged
nothing as guessed, meaning every bonded term resolved against `parm19.dat`.

## Four failures worth recording, because each one failed SILENTLY

Getting here took seven submissions. None of the failures produced an error that
named its own cause, and two of them exited 0.

1. **Hand-written PDB columns.** `antechamber` parses PDB by fixed column and
   rejects the file outright: *"Coordinates must be in Columns 31-38, 39-46 and
   47-54"*. A hand-built format string was one space off. Now written via gemmi
   with a post-write column assertion.

2. **`prepgen -i` requires an AC file, not a mol2.** Fed a mol2 it parses zero
   atoms, writes a prep file containing only three DUMM placeholders, and
   **exits 0**. `tleap` then builds the test peptide successfully, noting only
   `LYQ: no atoms` in its log among hundreds of lines. Two submissions passed
   their own "PASS" message while producing an empty residue. `prep_lyq.sh` now
   counts real atoms in the prep file and `check_lyq_charges.py` asserts all 17
   heavy atoms are present.

3. **`OMIT_NAME` matches by atom name, so cap names must not collide.** The
   obvious cap names all clash with residue atoms — ACE has C/O/CH3, NME has
   N/H/CH3, and the residue has C/O/N/H — so omitting the caps deletes the
   residue's own backbone. Caps renamed to Amber's CAY/CY/OY/NN/HN/CAT
   convention, with a collision assertion in the builder.

4. **`antechamber` infers the element from the atom NAME, ignoring the PDB
   element column.** The leucyl δ-carbons named `CL1`/`CL2` were read as
   **chlorine**, typed `DU`, and given +0.6125 charges. Every downstream tool
   then failed with `Atom type of DU does not exist in PARMCHK.DAT`, which names
   neither the atom nor the cause. Renamed to `CM1`/`CM2`.

   The guard for this is deliberately narrow. An initial version flagged every
   two-letter element prefix and immediately rejected `CAY` — and would have
   rejected `CA` and `NI`, which are correctly carbon and nitrogen in this very
   residue. Checking the mol2 antechamber actually produced showed that `CA`,
   `CB`, `CD`, `CE`, `CG`, `CAY` and `NI` were all typed correctly; only `CL`
   failed. So the guard lists `CL`, `BR`, `SE`, `SI` — prefixes observed to fail
   — and should be extended only when a new one turns up.

The general lesson, and it is the same one the earlier `set -eo pipefail` and
case-collision bugs taught: in this toolchain, exit status is not evidence of
success. Assert on the CONTENT of every intermediate file.

---

## Result: 3 x 20 ns complete

**The AF3 attack pose does not survive a force field.** 60 ns total, three
independent replicates from the same equilibrated start with fresh velocities.

| replicate | ns | start | mean | final quarter | min | frames <= 4 A |
|---|---|---|---|---|---|---|
| 1 | 20.0 | 8.07 | 7.41 | 7.15 +- 0.56 | 5.15 | 0 / 400 |
| 2 | 20.0 | 8.86 | 7.92 | 7.78 +- 0.39 | 5.97 | 0 / 400 |
| 3 | 20.0 | 8.07 | 7.28 | 7.30 +- 0.59 | 4.76 | 0 / 400 |

Across replicates: **7.41 +- 0.33 A**, against 4.86 A in the AF3 model used as
input and 4.6 A median across all 300 AF3 co-fold samples. **0 of 1200 frames
reach 4 A**; the closest single frame in 60 ns is 4.76 A.

The relaxation happened during **equilibration**, not production -- by the first
production frame the distance was already 8.1 A. Production shows the pose
staying relaxed, and the tight across-replicate SD (0.33 A) says 7.4 A is where
this system sits, not a transient.

### The pose is specific, just too far

Cys91 remains the nearest cysteine to the amine throughout, by a factor of two:

| | mean | min |
|---|---|---|
| **Cys91 (catalytic)** | **7.5 A** | **4.8 A** |
| Cys135 | 14.0 A | 10.3 A |
| Cys151 | 16.5 A | 6.7 A |

So AF3 did identify the right site -- the amine stays in the catalytic
cysteine's neighbourhood rather than wandering to another thiol. What it did not
get right is the distance, and 7.4 A is not attack geometry.

### The interface rearranges without dissolving

Ca RMSD reaches 3.5-6.0 A, which alone cannot distinguish "substrate came off"
from "interface rearranged". Contact counts separate the two:

| replicate | total contacts (start -> end) | retained from the AF3 pose |
|---|---|---|
| 1 | 477 -> 476 | 32% |
| 2 | 545 -> 517 | 30% |
| 3 | 502 -> 497 | 15% |

Total interface size is **flat**, while only 15-32% of AF3's specific contacts
survive. The complex persists -- in a different pose from the one AF3 predicted.
That is a more interesting outcome than either "stable" or "fell apart": AF3 got
the partner and the site right and the interface details wrong.

### UBE2W's C-terminal tail

Confirmed as the flexible element, with RMSF rising monotonically to the very
last residue:

```
res  141  142  143  144  145  146  147  148  149  150  151
RMSF 1.7  1.6  1.6  1.4  1.5  1.6  2.3  2.9  4.7  6.5  8.4   (A, replicate 1)
```

Tail (141-151) mean 1.5-3.1 A across replicates versus 1.2-1.6 A for the core
(1-140) -- 1.3x to 2.0x. This tracks AF3's own pLDDT collapse over the same
residues (91 for 1-140, falling to 41 at Cys151), so AF3 was correctly
signalling genuine disorder rather than failing to predict.

### What this means for the AF3 rounds

Three independent lines now agree. Round 5 showed the modification drives the
docking; rounds 2-3 showed the docking happens whether or not the site is
reactive; and MD now shows the specific sub-4 A contact is not a minimum of the
force field. Taken together: **AF3 places the modification in the active-site
neighbourhood, and the sub-4 A distances that looked like attack geometry are
noise on top of that placement.**

This does not make AF3 useless here -- identifying the correct partner and the
correct site is real. It makes the DISTANCE the wrong readout.

### Caveats

* 20 ns per replicate is short. The distance distribution is unimodal and the
  across-replicate SD is 0.33 A, so 7.4 A is well determined for this starting
  structure -- but a rare excursion to attack geometry on a longer timescale is
  not excluded by 60 ns.
* **No MM/GBSA.** 201 snapshots over the final 20 ns is Dominykas's protocol for
  a 100 ns run; 20 ns does not converge endpoint energies, and quoting an
  interaction energy from this run would produce a number that looks like a
  result.
* The alpha-amine is modelled **neutral**. It is the reactive species, but at
  pH 7.4 the majority species is protonated (pKa ~8.0-8.3). The cation is the
  obvious second condition and would make LYQ net +1.
* One condition, one substrate, one site. This is a feel-for-it run.

### Cost, measured

**131 ns/day** on one A100 for 25,392 atoms (96-residue substrate + 151-residue
UBE2W + 23,482 OPC waters + 63 NaCl pairs). So Dominykas's 3 x 100 ns is about
55 GPU-hours per construct -- one day across three GPUs.

Note on scheduling: the three replicates were originally written to run
sequentially in one 12 h job. At 131 ns/day that is 11.0 h of production plus
1.2 h of setup, which exceeds the limit, so replicate 3 would have been killed
mid-run. Replicates are independent by construction (same equilibrated restart,
`ntx=1` with `ig=-1` for fresh velocities), so they were split into separate
jobs and finished in 3.7 h.
