# Handover: modelling Ubl charging — UBE2W, UBA1, and engineered variants

For an agent picking this up cold. Everything below was established by testing, not
by reading documentation, and every trap listed here cost real time to find. Read the
**Silent failures** section before writing any code — most of the value in this
document is there.

## What this project is

UbyW installs an amino acid on a substrate lysine as an isopeptide (**XisoK**), which
presents a free α-amine as a neo-N-terminus. UBE2W then ubiquitylates that amine. The
question throughout has been: **what predicts whether a given site reacts?**

Answer so far, after five rounds: **nothing we have tried.** Sequence context gives no
motif; AF3 attack geometry is *anti*-predictive (the one working site ranked 7th of 8).
That history matters because it sets the prior for anything new: expect a negative, and
build the controls that make a negative interpretable.

Current work has moved to force-field MD, where the structure predictor supplies only
a starting geometry and the measurement is dynamic.

## Repositories

| repo | contents | state |
|---|---|---|
| **[PaulSchnacke/ubyw_charged_e2_af3](https://github.com/PaulSchnacke/ubyw_charged_e2_af3)** | Charged-E2 modelling, thioester parameterisation, MD. **Start here.** | `f5ea884`, 123 files |
| **[PaulSchnacke/ubyw_reactivity_analysis](https://github.com/PaulSchnacke/ubyw_reactivity_analysis)** | Rounds 1–5: sequence motif, XisoK ligand builder, site sweeps, first MD | `892c945`, 142 files |
| **PaulSchnacke/sortylation_success_analysis** | The prior negative on a different enzyme. Read its README for conventions. | — |

Both repos push via `./push.sh`, which **verifies the remote hash matches local HEAD**
— use it rather than bare `git push`; a silent push failure has happened.

### Files that matter most

**Building an XisoK covalent modification for AF3**
* `ubyw_reactivity_analysis/af3/make_xisok_ccd.py` — writes a `userCCD` for any
  X-isopeptide-lysine. Builds from SMILES, **carbonyl-first** (atom naming follows
  parse order, so writing it intuitively gives wrong names), validates
  stereochemistry from embedded 3D coordinates.
* `ubyw_charged_e2_af3/ccd_valence.py` — **run this on every CCD before submitting.**
  Sums declared bond orders plus the external bonds AF3 will create. It caught a
  pentavalent carbon (rendered as a cyclopropene in a viewer), a tetrahedral
  "thioester", and two over-valent nitrogens — all in files that passed AF3's own
  schema validation.
* `ubyw_charged_e2_af3/build_jobs_v2.py`, `graft_msa.py`, `qc_sweep.py`

**MD with a custom residue**
* `md/build_lyq_model.py` + `md/prep_lyq.sh` — the **LYQ** residue (isopeptide-linked
  lysine). Validated: integer charge, amide half agrees with published ALY
  (Nε-acetyllysine, from `ff19SB_modAA`) within 0.25 e, nucleophile N at −0.88 e.
* `md/make_thioester_link_frcmod.py` — the `C–S` terms for a protein–protein thioester.
* `md/prepare_charged.py`, `md/build_thioester_md.sh` — the charged-state pipeline.
* `md/check_frcmod_attn.py` — **the gate that matters.** See below.

---

## How to model a covalent connection in AF3

### The rule that governs everything

**AF3 3.0.1 silently discards covalent bonds between two POLYMER chains.** It logs
`Reducing number of bonds ... N are polymer-polymer` in
`structure_cleaning.py`, exits **0**, and leaves the atoms tens of angstroms apart
(measured: 26.25 Å). Protein→**ligand** bonds are kept.

So any inter-protein covalent bond must be expressed with at least one partner as a
ligand. For ubiquitin's C-terminus onto a cysteine:

```
protein chain U = Ub 1-74          (ends at Arg74)
ligand   chain G = GLY             (Gly75)
ligand   chain H = GLY             (Gly76)
bondedAtomPairs:
  [["U",74,"C"], ["G",1,"N"]]      extend the chain
  [["G",1,"C"],  ["H",1,"N"]]
  [["H",1,"C"],  ["B",91,"SG"]]    the thioester
```

Verified to give continuous backbone spacing (Cα–Cα 3.77/3.94/3.61 Å) and all four
bonds formed. **Always grep the log for `Reducing number of bonds`** — that one line
is the difference between a real model and a meaningless one.

### For an isopeptide (XisoK) on a lysine

Simpler: the modification is a **ligand**, so one bond suffices —
`[["A",<res>,"NZ"], ["L",1,"C01"]]`. Give the CCD via `userCCD`.

Note AF3 does this split natively for polyubiquitin: it puts Gly76 in its own
single-residue chain and encodes Gly75.C→Gly76.N as an inter-chain bond.

### Numbering

Constructs in this study carry an **N-terminal proline** to block the native
N-terminus (which is UBE2W's canonical target). Consequence: **construct residue
number = native number + 1.** SUMO2 native K11 is construct K12. This applies to
essentially every construct in the paper's supplementary table — an earlier round
applied it to SUMO2 only, which gave Ran a free native N-terminus and a competing
nucleophile. Check the shipped JSON, not the docstring.

### AF3 limitations to state up front

* **It discards ligand stereochemistry.** Opposite-enantiomer inputs give
  *bit-identical* coordinates. Any D/L comparison through AF3 is pseudo-replication.
* **It does not refine bond lengths.** Isopeptide C–N comes out at 1.07 ± 0.09 Å
  (amide is ~1.33); in one charged model 0.84 Å. Connectivity only.
* **It has no thioester chemistry.** The acyl carbon pyramidalises to 338.6 ± 1.7°
  between sp² (360°) and sp³ (328.5°) despite correct input hybridisation.
* Score **ensembles, not top-ranked models** — 5 seeds is too few; the spread is large.
  20–25 seeds × 5 samples. Bootstrap any ranking before believing it.

### Running AF3 on Euler

```bash
SIF=/cluster/apps/nss/alphafold/containers/AlphaFold3/af3.sif
DB=/cluster/project/alphafold/alphafold3
MODELS=$HOME/models
RUN=/app/alphafold/alphafold-3.0.1/run_alphafold.py
export APPTAINERENV_XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"
apptainer exec --nv --bind $DB:/db --bind $MODELS:/models --bind $WORK:/work:rw \
  $SIF python3 $RUN --json_path=/work/job.json --output_dir=/work/out \
  --db_dir=/db --model_dir=/models --flash_attention_implementation=xla
```

`--nv` is **required for inference** and must be omitted for the CPU-only MSA stage.
Account `es_lang`. MSA search is 4–6 h and I/O-bound; inference is ~10 min on a GPU.

**Reuse MSAs aggressively.** Installing a modification, changing a ligand, or moving
the site does not change the protein sequence being searched — so every site of one
protein shares one MSA. `graft_msa.py` copies them in; then `--norun_data_pipeline`
skips straight to inference. Archived searches are in `~/ubyw_msa_archive` and
`/cluster/scratch/schnpaul/ubyw/msa*` (scratch is purged after 15 days; home is not).

---

## How to model a covalent connection in Amber MD

### Thioester parameters — settled, use these

**Bonded terms: GAFF2**, transferred onto protein atom types by
`md/make_thioester_link_frcmod.py` (reads `gaff2.dat` at runtime, so no transcription
error):

```
BOND   C -S    199.66 kcal/mol/A^2   1.8104 A
ANGLE  O -C -S  76.330   123.380 deg      <- sp2 carbonyl
       CT-C -S  62.840   113.460          (also XC, 2C, 3C)
       C -S -CT 93.670    99.120          (also XC, 2C, 3C)
DIHE   X -C -S -X   2   6.200  180.0  2.0 <- holds the thioester planar
```

**Charges: Oda 2013** (*Chem. Lett.* **42**:1206–1208, doi:10.1246/cl.130517),
Figure 5, RESP/HF-6-31G(df,p) on acetylcysteine:

| atom | charge (e) |
|---|---|
| thioester **S** | **−0.307889** (a standard Cys SG is −0.1081) |
| acyl **C** | **+0.479917** |
| acyl **O** | −0.414690 |
| CT / HC(×3) | −0.156886 / +0.073065 |
| CB / HB(×2) / CA | +0.076612 / +0.074370 / +0.024622 |

The charges are not optional detail — if the question is whether a nucleophile
approaches an electrophile, the junction charges *are* the measurement.

**Why this split.** Oda also publishes bonded terms, and its equilibrium angles look
non-physical (O–C–S at 99.7°). Tested: both sets minimise to correct sp² geometry
(Oda 119.7°, GAFF2 124.4°), so Oda's are *effective* values — but Oda carries
~8 kcal/mol more residual angle strain (`ANGLE 11.10` vs `2.93`). GAFF2's bonded
terms sit in a relaxed minimum. Oda's stiffer CCSD(T)-fitted torsion (8.941 vs 6.200)
is worth testing as a variant since planarity is part of what gets measured.

**`parm19.dat` has no `C–S` bond at all** — no standard amino acid carries a
thioester. This is why an **isopeptide** link needs no new parameters (`C–N` amide
terms are standard) while a **thioester** does. If someone's polyubiquitin MD "just
worked", that is why.

### Creating the bond in tleap

```
source leaprc.protein.ff19SB
loadamberparams thioester_link.frcmod
loadamberparams frcmod.lyq
loadoff lyq.lib
sys = loadpdb system.pdb
bond sys.<UbGly76> sys.<Cys91>.SG      # explicit inter-chain bond
```

**Reverse the AF3 ligand split first.** The Gly75/76-as-ligands topology exists only
to work around AF3; Amber has no such restriction, so promote them back into the
ubiquitin chain (`md/prepare_charged.py`) and **remove Arg74's OXT** — otherwise the
carbonyl carries both OXT and the amide nitrogen. Use **CYX** (cysteine without the
SG hydrogen), or tleap caps the sulfur.

Unlike AF3, **tleap fails loudly** when a bond parameter is missing:
`Could not find bond parameter for atom types: C - S`, 7 errors, no topology written.
That failure mode cannot reach a trajectory.

### frcmod format traps

Amber rejects an entire frcmod with `Could not load parameter set` — which reads like
a *missing parameter* problem — when the formatting is wrong. Two ways to trigger it:
a **multi-line title** (exactly one line allowed), and **no `MASS` section** (must be
present even when empty). Compare against a file tleap accepts.

Also: **ff19SB does not type every sp3 carbon `CT`.** The α carbon is `XC`, β/γ are
`2C`/`3C`. Emit the angle for each variant or tleap asks for `XC-C-S` and stops.

---

## Silent failures — read this before writing code

Every one of these exited 0, or passed the check available at the time.

1. **AF3 dropping a polymer–polymer bond.** Exit 0, atoms 26 Å apart.
   → grep `Reducing number of bonds`.
2. **`parmchk2` emitting zero force constants.** It flags guessed terms `ATTN, need
   revision` and sets `k = 0.00`. For a thioester built with `-at amber`, **all six
   flagged terms were the thioester itself** — it would have simulated as two
   unconnected fragments while writing a valid trajectory. The cause was the atom-type
   flag: use **`-at gaff2`** (sulfur → `ss`, acyl C → `c`), not `-at amber`.
   → `md/check_frcmod_attn.py` refuses any residue with a zero-force-constant term at
   the reactive centre. **Run it on every custom residue.**
3. **`tleap` accepting a topology proves nothing.** Read the bond back out of the
   `parm7` and require a non-zero force constant. Instant, no queue:
   check `BONDS_WITHOUT_HYDROGEN` + `BOND_FORCE_CONSTANT` for the atom pair.
4. **`antechamber` infers element from atom name.** A two-letter carbon name can be
   read as a halogen. And `prepgen`'s `MAIN_CHAIN`/`OMIT_NAME` match **by name**, so
   RDKit's generic `C1/N1/UNL` output silently produces a garbage residue tree after
   two earlier stages reported success. Assign PDB-convention names explicitly.
5. **Case-collisions on macOS/APFS and in AF3 output directories.** Variant names
   differing only in case overwrite each other silently. Bit this project three times.
6. **CIP descriptor letters are not L/D.** L-Cys is (R) because sulfur outranks the
   carboxyl. Verify configuration by **signed volume against a reference**, never by
   the letter. This produced a D-amino acid that prose had asserted was correct.
7. **`$SLURM_SUBMIT_DIR` is unset** in this job wrapper, so `cd $W/...` moves away
   from staged inputs and `cp` silently finds nothing. Capture `SRC="$PWD"` first and
   assert every input exists. Also: `#SBATCH` lines must be first, or they become
   comments and the job runs on the login node with no GPU.
8. **The cluster's bare `python3` has no numpy or gemmi.** A QC step crashed *after*
   AF3 had written 101 models. Write cluster-side QC in the standard library only.
9. **Walltime arithmetic.** Measure ns/day on the *actual* system first (111 ns/day
   for ~110k atoms on a TITAN RTX → 20 ns ≈ 4.3 h). Replicates are independent by
   construction (same equilibrated restart, `ntx=1`, `ig=-1`), so run them as
   **parallel jobs** — sequential replicates in one job has twice exceeded the limit.
10. **Write results to `$HOME`, not scratch**, and run QC incrementally on the
    cluster. VPN drops then cost nothing.

---

## Current state and the immediate next questions

### Running now

* **Charged-state MD** (job `bf070c11`): SUMO2(K11 XisoK) + UBE2W + Ub 1-76 with a
  real thioester. First MD of the state the project actually needs — earlier runs
  were **apo** (no ubiquitin at all) and were stopped as uninterpretable, correctly:
  removing ubiquitin empties the groove the nucleophile must reach into.
* **Ub(1-75) engineering question** (job `05bbee62`): four AF3 jobs, 25 seeds each —
  charged Ub 1-76 vs 1-75, and non-covalent 1-76 vs 1-75.

**A finding already, before those land:** you cannot get Ub(1-75) by truncating the
Ub(1-76) model. Deleting Gly76 leaves the thioester at **4.20 Å** — not bonded. The
shorter variant is a genuinely different geometry and needs its own prediction. If it
turns out UBE2W cannot reach Gly75's carbonyl, that is a structural explanation for
why the engineering failed on the UBE2W side rather than at UBA1.

### Next: UBA1 with Ub(1-75)

Two sites to model, and they are chemically different:
* the **adenylation site** (Ub C-terminal carboxyl + ATP), and
* the **active-site cysteine** (thioester).

The thioester machinery above transfers directly to the second. The adenylation site
needs a different treatment — a carboxyl/AMP ester, not a thioester — and its
parameters have not been looked for yet.

### Topics deliberately not researched yet

Paul has data and specific papers in mind; **ask before searching**:
* UBE2W's lack of sequence specificity; NEDD8 also accepted.
* A recent paper on Ubl charging.
* UBA1 engineered to charge SUMO and other Ubls.

### What to do about a negative

Five rounds have been negative. The project's convention — worth keeping — is to
retract explicitly rather than quietly amend, keep the git history, and record what
an intermediate state claimed. Two retractions are recorded in
`ubyw_reactivity_analysis/README.md`. A negative with good controls is a result; a
positive from a lucky subset is a liability.
