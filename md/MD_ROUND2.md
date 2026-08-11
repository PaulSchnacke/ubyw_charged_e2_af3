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

---

## Update: the hard residues, and why "tleap accepted it" is not a pass

All three parameterised, and all three would have run. Two of them should not.

| residue | antechamber | prepgen | tleap | net charge | ATTN terms | **verdict** |
|---|---|---|---|---|---|---|
| **LYP** (product) | ok | 47 atoms | OK | −0.000002 | **0** | **usable** |
| LYT (thioester) | ok | 18 atoms | OK | 0.000000 | 6 | **NOT usable** |
| LYX (tetrahedral) | ok | 22 atoms | OK | 0.000001 | 3 | **NOT usable** |

Every signal we had said all three were fine: AM1-BCC charges assigned, integer
residue charge, tleap built a test peptide and reported `Total unperturbed charge:
0.000000`. That is a misleading pass.

`parmchk2` flagged terms as `ATTN, need revision` and set each one to a force
constant of **zero**. For LYT, all six are the thioester itself:

```
BOND   C -S        k = 0.000     <- the thioester bond
ANGLE  C -S -CT    k = 0.000
ANGLE  O -C -S     k = 0.000
ANGLE  CT-C -S     k = 0.000
DIHE   O -C -S -CT k = 0.000     <- rotation about C-S
DIHE   CT-C -S -CT k = 0.000
```

Zero force constant means **no restoring force**. The C–S bond has no equilibrium
length, rotation about it is free, and there is nothing holding the carbonyl in
plane. LYX is the same story at the tetrahedral carbon: `OH-CT-S`, `NT-CT-S` and
`NT-CT-OH` are the angles that *define* the sp3 centre, and all three are zero, so
it can flatten or invert without penalty.

A trajectory on either residue would run to completion and produce numbers that
describe nothing. This is the same class of failure as the AF3 polymer-bond drop:
exit 0, plausible output, silently meaningless.

`check_frcmod_attn.py` is now the gate. It fails when any zero-force-constant
guessed term involves an atom of the reactive centre, and it was verified on the
real files: LYP passes with 0 ATTN, LYT fails naming all six terms.

### What this means for the two intermediates

The thioester and tetrahedral states **cannot be measured by classical MD** with
transferable parameters, which was the risk flagged before running. Options, in
order of cost:

1. **Accept the limitation.** The charged state's behaviour is not answerable this
   way, and AF3's 338.6° pyramidalisation already suggested it has no clean answer
   either.
2. **QM-derived parameters** for the thioester (a Hessian-based fit at, say,
   B3LYP/6-31G* on the model compound). Days of work, and defensible.
3. **QM/MM**, treating the reactive centre quantum-mechanically. The correct
   method for a bond-forming question, and a substantially larger project.

My recommendation is (1) for now and (2) only if the thioester state becomes the
central question. The product (LYP) is genuinely usable and worth running, since
it asks a well-posed classical question: once ubiquitin is attached, is the product
complex stable, or does it fall apart?
