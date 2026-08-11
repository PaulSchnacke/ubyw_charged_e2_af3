# The charged state runs: thioester stable, amine engaged

First MD of the state this project has needed throughout — SUMO2(K11 XisoK) + UBE2W +
ubiquitin 1-76 joined by a **real thioester** to Cys91. Earlier MD rounds were **apo**
(no ubiquitin at all) and were stopped as uninterpretable; that was correct.

## Through minimisation, heating and 500 ps NPT equilibration

| observable | first | last | mean | reference | verdict |
|---|---|---|---|---|---|
| thioester C–S | 1.81 Å | 1.85 Å | **1.82 Å** | GAFF2 r_eq 1.8104 Å | held |
| O=C–S angle | 122.45° | 122.07° | **123.17°** | sp² planar ~123° | **planar** |
| attack N(I)→acyl C | 3.66 Å | 3.28 Å | 3.66 Å | ≤4 Å attack window | **engaged** |

Three things this establishes.

**The GAFF2 transfer holds a real bond.** 1.82 Å mean against a 1.8104 Å equilibrium,
stable across the whole run rather than drifting. The parameters were transferred from
GAFF2's `c-ss` family onto ff19SB protein atom types, which is the step that could
have failed silently and did not.

**The force field has chemistry AF3 lacks.** The thioester carbon stays **planar at
123.2°**. AF3 pyramidalised the same carbon to 338.6 ± 1.7° — between sp² (360°) and
sp³ (328.5°) — despite receiving correct input hybridisation. This is the concrete
argument for MD as the measurement and AF3 as only a starting-structure generator.

**The nucleophile is inside the attack window and moved closer.** The XisoK α-amine
sits 3.28–3.66 Å from the acyl carbon it must attack, and it **approached** during
equilibration (3.66 → 3.28 Å) rather than relaxing away. Contrast the apo runs, which
relaxed to 7.4 ± 0.3 Å with **0/1200** frames inside 4 Å. Ubiquitin's presence is what
changes this — which is exactly why running without it answered nothing.

## What this is not

Equilibration, not production. 500 ps with the barostat still settling is not evidence
of a stable pose, and a starting structure chosen for near-attack geometry is expected
to *begin* near attack. The question production answers is whether 3.3 Å is a genuine
residence or a memory of the starting structure. Three 20 ns replicates follow.

Nor is it evidence about reactivity yet: with fixed connectivity, classical MD cannot
show the transfer occurring. What it can measure is whether the geometry required for
transfer persists — and how that differs between a site that works (K11) and one that
does not (K21).

## Cost of getting here: six submissions, four real errors

Every one was caught by a loud failure rather than producing a plausible trajectory.

| # | error | how it announced itself |
|---|---|---|
| 1 | `prep_lyq.sh` needs `build_lyq_model.py` to have RUN | "run build_lyq_model.py first" |
| 2 | that builder takes a FILENAME, not a directory | `IsADirectoryError` after the stereo check passed |
| 3 | gate tested for `frcmod.lyq`/`lyq.lib`; it writes `lyq.frcmod`/`lyq.prep` | aborted a job whose residue had validated perfectly |
| 4 | LYQ leucyl methyls are `CM1`/`CM2`, not `CL`/`CM` | "Atom .R\<LYQ 12\>.A\<CL 41\> does not have a type" |
| 5 | thioester Gly76 mapped to the CHARGED C-terminal variant | "Could not find angle parameter for atom types: S - C - O2" |
| 6 | ff19SB has several sp3/α carbon types | asked for `XC-C-S`, then `S-C-CX` |

Error 5 is the one worth remembering: removing OXT from the coordinates was **not
enough**, because tleap picks the C-terminal variant from chain **position**, not from
which atoms are present. `addPdbResMap { { 1 "GLY" "GLY" } { 0 "GLY" "GLY" } }` fixes
it, and is surgical here because chain U is the only chain ending in glycine (A ends
TYR96, B ends CYS151).

Error 4 now has a permanent guard: `prepare_charged.py` asserts the fused residue's
heavy-atom set equals what `lyq.prep` declares, so a rename in either place fails
locally instead of on the cluster.


---

# Production, 11.1 ns of replicate 1: the nucleophile stays engaged

Analysis of the running trajectory (222 frames at 50 ps, verified against the mdin).

| observable | first | mean | final quarter | reading |
|---|---|---|---|---|
| thioester C–S | 1.81 Å | **1.81 ± 0.04 Å** | 1.81 Å | rock stable at GAFF2's 1.8104 Å |
| O=C–S angle | 125.4° | **124.1°** | 124.4° | sp² planar throughout |
| attack N(I)→acyl C | 3.87 Å | 3.67 Å | **3.39 Å** | **tightens**; 176/222 frames (79%) ≤4 Å |
| SUMO2 Cα RMSD | 0 | 2.67 Å | 3.05 Å | settles, no unfolding |
| Ub Cα RMSD | 0 | 0.90 Å | 1.10 Å | ubiquitin very rigid |

## The comparison that matters

| | apo (stopped runs) | charged (this run) |
|---|---|---|
| final attack distance | **7.4 ± 0.3 Å** | **3.39 Å** |
| frames ≤4 Å | **0 / 1200** | **176 / 222 (79%)** |

The apo complex relaxed away from attack geometry and never returned. With ubiquitin
bonded, the nucleophile *closes* on the acyl carbon instead. Ubiquitin's presence is
the entire difference — removing it did not simplify the system, it removed the
interaction being measured.

## Two things this does and does not show

**Does:** the force field reproduces thioester chemistry AF3 cannot. The acyl carbon
stays planar at 124° across 11 ns where AF3 pyramidalised the same carbon to
338.6 ± 1.7° total angle. The bond neither stretches nor drifts.

**Does not:** show reactivity. With fixed connectivity, classical MD cannot model the
transfer. What it measures is whether the geometry required for transfer *persists* —
and here it does, for one site, in one replicate.

## The cautious panel

Native contacts fall to ~41% (substrate–enzyme) and ~44% (ubiquitin–enzyme) of their
starting values. Two reasons not to read that as dissociation, and one reason to watch
it:

* the reference frame is the **unrelaxed AF3 pose**, so a large early drop is expected
  as an idealised prediction settles into solvent;
* the substrate interface is **flat** after that initial settling (mean 371 vs final
  quarter 388 pairs) — settled, not leaving;
* but the **ubiquitin interface is still declining** (mean 804 vs final quarter 656).
  If that continues in the full 20 ns × 3 replicates, the engaged geometry may be
  transient. That is exactly the question production answers, and it is why the
  replicates matter rather than this single partial run.
