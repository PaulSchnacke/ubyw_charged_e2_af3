# MD round 2: AF3 as prior, MD as the measurement

Six systems requested. Three had force-field parameters ready and are running;
three needed new ones and are best-effort, at your direction.

## Starting structures — selected from the sweep by closest approach

| system | AF3 start | source model |
|---|---|---|
| K11 XisoK | **3.78 Å** | `seed-20_sample-3` |
| K21 XisoK | **1.77 Å** | `seed-15_sample-2` |
| K11 charged | 3.78 Å | `seed-19_sample-3` |
| K11 tetrahedral | 2.50 Å | `seed-9_sample-2` |
| K11 product | 3.24 Å | `seed-20_sample-3` |
| K11 Lys control | **12.82 Å** | `seed-13_sample-1` |

### The control cannot start from near-attack, and that is a result

You asked for near-attack starts. Five systems have one. The unmodified-lysine
control does not, and it is not for want of looking: across **200 AF3 models** the
unmodified construct residue 12 (native K11) never comes closer than **12.8 Å** to
Cys91 SG. So this trajectory asks a different question — starting from where AF3
actually puts an unmodified lysine, does it approach at all under a force field?

An error worth recording: my first control pick was the model with the smallest
`nearest_lys_nz_d`, which was **4.15 Å** — but that value belongs to construct
residue 22 (native K21), a *different site*. The stored QC only kept the nearest
lysine, not residue 12. Re-measuring residue 12 specifically across all 200 models
gave the real 12.8 Å. Had I not checked, the control would have been the wrong
lysine at the wrong site.

## Running now: the three with validated parameters

`k11_xisok`, `k21_xisok`, `k11_lyscontrol` — ff19SB/OPC, minimise → heat → NPT
equilibrate → **3 × 20 ns** production each, Langevin thermostat and Monte Carlo
barostat. The modified systems use **LYQ**, the isopeptide residue built and
charge-validated in the earlier round (integer charge, amide half agreeing with
published ALY parameters, nucleophile correctly polarised). The control needs no
custom residue.

`prepare_system.py` now takes the site as an argument, so one script serves K11
(residue 12) and K21 (residue 22). The control gets its own `prepare_control.py`
rather than a flag that would weaken the ligand assertion in the other script.

## Best-effort: the three that need new parameters

Honest confidence ordering, because these are not equivalent:

| residue | chemistry | confidence |
|---|---|---|
| **LYP** (product) | isopeptide + a second ordinary amide from Ub Gly76 | **defensible** — amides are the best-covered motif in any protein force field |
| **LYT** (thioester) | Cys-S-C(=O)-Gly | **approximate** — a thioester is not standard ff19SB chemistry; the C–S bond, reduced C=O order and near-free rotation about C–S are all poorly described by transferable parameters |
| **LYX** (tetrahedral) | sp3 C bearing S, N, O⁻ and CH₂ | **not defensible as kinetics** — a classical force field cannot form or break the bonds that define an intermediate, so a trajectory shows only whether the geometry is *mechanically stable* under ff19SB |

If LYT and LYX fail parameterisation, they fail; we decide afterwards whether to
pursue QM-derived charges.

### A D-amino acid caught by the structural check

The first LYP build contained a **D-leucine**. The CIP letters looked unremarkable
(`['S','R']`) and my prose had asserted the SMILES ordering gave L at both
centres — an assumption, not a test. Computing the signed volume of the
(N, carbonyl-C, side-chain-C, H) tetrahedron at each alpha carbon and comparing
against L-Lys/L-Cys/L-Leu references showed one centre with the wrong sign.

The check is now a gate inside the builder, and it was verified in both
directions: it rejects a deliberately planted D compound and passes all three real
ones. This is the third time in this project that a CIP letter has failed as an
L/D proxy — L-Cys is (R) while L-Lys is (S), same geometry, different letter.

## One configuration failure caught

The first MD submission died at exit 1 with an empty log. Cause: `$SLURM_SUBMIT_DIR`
is **not set** in this wrapper (it runs the script directly rather than through
`sbatch`), so `cd $W/md` moved away from the staged inputs and the copy became a
silent no-op with `2>/dev/null` hiding it. Now the source directory is captured
before any `cd`, and every input is asserted present with the source path in the
error message.

Same family as the earlier `#SBATCH`-not-first failure: this wrapper is not
`sbatch`, and assuming otherwise costs a submission each time.

## What the measurement is

For each modified system, the question MD answers that AF3 cannot: **does the
near-attack pose survive a force field, or relax away?** In the earlier round the
answer for the uncharged K11 complex was relaxation to 7.4 ± 0.3 Å with 0/1200
frames inside 4 Å. If K11 and K21 now differ in how long they hold, that is a
dynamic signal AF3's static reach metric could not see — and the reach metric was
anti-predictive (K11 ranked 7th of 8), so there is room for MD to say something
different rather than merely confirm it.
