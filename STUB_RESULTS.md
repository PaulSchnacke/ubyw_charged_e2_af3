# Stub validation: the valence fix works, and a new AF3 limitation

Job `f2ea3a80`, 3 jobs × 1 seed × 5 samples, single-sequence MSAs, ~9 min GPU each.
MSA search `4e124915` also completed — much faster than the 4-6 h estimated.

**Verdict: the chemistry is correct, submit production.** One new finding worth
recording, which does not block anything but does constrain what these models can be
used for.

## 1. No bonds dropped

```
grep -i "Reducing number of bonds" *.log  ->  NONE
```

Both declared bonds present in every model, and the `covale` records name the intended
atoms:

```
covale1   U 74 C (ARG 74)  ->  G 1 N1        extend the backbone
covale2   G 1 C2 (UBGG)    ->  A 632 SG      the thioester   [Ub 1-76]
covale2   G 1 C1 (UBG1)    ->  A 632 SG      the thioester   [Ub 1-75]
```

**Chain renumbering did NOT happen this time.** The ligand stayed its own chain `G`
(8 atoms for UBGG, 4 for UBG1) and ubiquitin stayed 74 residues. In the earlier charged
round AF3 absorbed the first `GLY` ligand into chain `U` and left the second as chain
`UA`. The difference is that a *custom multi-atom CCD* is not absorbed the way a
single-residue standard `GLY` is. QC must still resolve chains from the covale records
rather than assume either layout.

## 2. The valence fix works — the thioester is now planar

| | C–S (Å) | **bond-angle sum at acyl C** | OXT present |
|---|---|---|---|
| Julian run (bare `GLY`) | 1.723 | **328.4°** | yes, 50/50 models |
| **stub, Ub 1-76 (`UBGG`)** | 1.612 | **356.9°** (355.2–360.0) | no |
| **stub, Ub 1-75 (`UBG1`)** | 1.514 | **359.3°** (358.2–360.0) | no |
| ideal sp² thioester | ~1.78 | 360° | — |
| ideal sp³ | — | 328.5° | — |

This is the direct confirmation that Bug 2 was the cause. Removing OXT moved the acyl
carbon from textbook tetrahedral (328°) to essentially planar (357-359°). Both tail
lengths behave identically, so `UBG1` is as sound as `UBGG`.

Note this also revises the earlier interpretation. `SWEEP_RESULTS.md` reported the
charged species at **338.6 ± 1.7°** and offered "AF3's prior for that local environment"
as a possible explanation, since the input hybridisation was correct. It was not AF3's
prior — it was residual over-valence. With the valence genuinely correct, AF3 *does*
build a planar sp² thioester. That is a point in AF3's favour and against our earlier
chemistry.

## 3. NEW: AF3 distorts the atom it bonds — the bonded Cys is the sole outlier

This was not looked for; it fell out of a control check.

`HANDOVER_UBL_CHARGING.md` already says AF3 "does not refine bond lengths" for the
*isopeptide it creates*. The new observation is that the distortion propagates **into
the pre-existing residue**, corrupting a bond AF3 was not asked to touch:

```
STUB (UBA1, thioester at Cys632)        JULIAN RUN (UBE2W, thioester at Cys91)
  A:632   1.237  <-- bonded Cys           B:91    1.279  <-- bonded Cys
  B:151   1.805                           A:49    1.805
  B:135   1.808                           B:135   1.816
  A:431   1.811                           B:119   1.816
  ...19 more, all 1.80-1.83...            ...all 1.80-1.82...
```

24 cysteines in the stub, 6 in the Julian run. **In both, the single cysteine carrying
the declared bond is the only one with a distorted CB–SG bond** — 1.24 and 1.28 Å against
an ideal 1.81 Å, while every untouched cysteine sits at 1.80-1.83 Å. Two different
proteins, two different residue numbers, two MSA regimes (single-sequence and full),
same result.

The general pattern, from bonds grouped by who created them:

| bond class | measured | ideal | verdict |
|---|---|---|---|
| ligand internal (from our CCD) | C=O 1.221, 1.277 | 1.23 | **preserved** |
| protein backbone (AF3's own) | CA–C 1.528 ± 0.01 (n=60) | 1.52 | **correct** |
| AF3-created thioester | 1.565 | 1.78-1.81 | too short |
| AF3-created peptide | 1.817 | 1.33 | too long |
| **pre-existing bond at a bonded atom** | **CB–SG 1.237** | **1.81** | **corrupted** |

So AF3's backbone geometry is fine and our CCD survives intact; the damage is local to
the atoms named in `bondedAtomPairs` and spreads one bond further than documented.

**Consequence.** Do not read any bond length or angle *at or adjacent to* a declared
bond as physical — including in the residue AF3 already knew about. Distances between
atoms that are not part of a declared bond, and everything else in the structure, remain
usable. This strengthens rather than weakens the existing plan: it is another reason the
thioester geometry has to come from a force field, and it means the MD starting
structure needs the junction geometry rebuilt (tleap does this anyway when it applies
the real `C–S` parameters, so the practical cost is nil).

It also means the `1.57 Å` thioester here should **not** be reported as "AF3 gives a
short C–S". The right statement is that AF3 treats `bondedAtomPairs` as connectivity and
does not refine the junction — already this project's position, now with a sharper
boundary on how far the unrefined region extends.

## 4. MSAs are deep

| chain | length | unpaired | paired | templates |
|---|---|---|---|---|
| UBA1 (A) | 1058 | 16,643 | 50,000 | 4 |
| ubiquitin (U) | 74/75/76 | 7,905-7,964 | 50,000 | 4 |
| UBE2W (B) | 151 | 14,993 | 50,000 | 4 |

All three donors complete in `~/ubyw_uba1_msa/`. Ubiquitin's alignment differs slightly
by tail length (7,905 / 7,934 / 7,964), which is expected — the query differs — and is
why the donors were built per length rather than grafted from one.

## 5. Sizing for production

~9 min/job at 1 seed on 40 GB gpumem, 1283-1285 residues. 25 seeds is roughly 25× that
per job before batching, so run variants as **parallel jobs**, not sequentially in one —
sequential replicates in one job has twice exceeded walltime in this project.
