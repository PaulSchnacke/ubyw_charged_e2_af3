# `jobs_ub75/` carries both Julian-run bugs, and three more found while fixing it

**The Ub(1-75)-vs-Ub(1-76) UBE2W comparison was run on chemistry known to be wrong.**
`jobs_ub75/` was written *after* the valence fix but reverted to the pre-fix components.
Rebuilt jobs are in `jobs_ub75_v2/` and need rerunning before those numbers mean anything.

## The two original bugs, verified

Checked with the repo's own `ccd_valence.py` under the exact external bonds each job
declares — the check that matters, since each half is fine alone and only the sum is wrong.

**Bug 1 — `LIG-1` retains the spurious `C01–H01` bond.**

```
jobs_ub75 LIG-1 : O01=C01, C01-C02, C01-H01   -> C01 at 5/4  OVER-VALENT
ccd_v2/LisoK    : O01=C01, C01-C02            -> C01 at 4/4  with the bond
```

The isopeptide acyl carbon reaches five coordination once the bond to the substrate
lysine Nζ is added. This is the bug that rendered the acyl carbon as a cyclopropene.

**Bug 2 — the terminal glycines are bare `GLY`, which retains OXT.**

```
GLY N  at 4/3  OVER-VALENT   (also takes the peptide bond from Ub Arg74)
GLY C  at 5/4  OVER-VALENT   (C, O, OXT, CA, plus the external bond to SG)
```

Note **both** atoms are over-valent, not only the acyl carbon — worse than the original
write-up described. `ube2w_ub76_charged.json` declares *two* bare `GLY` ligands, so it
carries the fault twice.

The consequence is measured: `STUB_RESULTS.md` shows 328° (sp³) with bare `GLY` versus
357–359° (planar sp²) with the no-OXT CCDs.

## What the rebuild changes

| | old | new |
|---|---|---|
| isopeptide moiety | `LIG-1` with `C01–H01` | `ccd_v2/LisoK` (no `C01–H01`) |
| Ub 1-75 tail | bare `GLY` | **`UBG1`** (Gly75, no OXT) |
| Ub 1-76 tail | bare `GLY` + bare `GLY` | **`UBGG`** (Gly75-Gly76 dipeptide, no OXT) |
| bonds declared | 4 (1-76) | 3 |

Collapsing two `GLY` ligands into one `UBGG` also removes a declared bond — the internal
Gly75–Gly76 peptide now lives inside the CCD, so it is one fewer thing for AF3 to discard.

Final audit of the rebuilt jobs, every ligand atom under its declared external bonds:

```
ube2w_ub75_charged_v2   L/LIG-1 ext={C01:1}          all within valence
                        G/UBG1  ext={N1:1, C1:1}     all within valence
ube2w_ub76_charged_v2   L/LIG-1 ext={C01:1}          all within valence
                        G/UBGG  ext={N1:1, C2:1}     all within valence
```

The `free` jobs have no ligands and no bonds; they are copied unchanged.

## Three further bugs found while fixing it — all in my own new code

Worth recording because each is the same silent-success shape, and two were caught only
because the gate ran on a job with *two* custom ligands rather than one.

**(a) `external_bonds()` pooled external bonds by atom name across all ligands.** With a
single custom ligand per job that is harmless, which is exactly how it survived the 12
UBA1 jobs. With `LisoK` + `UBGG` in one job it hands each ligand the other's atom names,
and `ccd_valence` raised `external bond names undeclared atom 'N1'`. Now keyed by ligand
chain, then atom.

**(b) A job can declare a `ccdCode` that its own `userCCD` does not define.** AF3 resolves
a ligand by the `_chem_comp.id` *inside* the CIF, not by the filename — and in `ccd_v2/`
the two differ: `LisoK_userCCD.cif` declares `_chem_comp.id LIG-1`. My first rebuild
asked for `ccdCodes: ["LisoK"]` against a block defining `LIG-1`, so AF3 would have failed
to find the component. There is now one shared `resolve_ccd()` that both the valence gate
and the embedder use, and `embed_ccd` raises if the declared code and the file's declared
id disagree.

**(c) The filename is not the component id.** Generalising (b): never derive a `ccdCode`
from a filename stem. `rebuild_ub75_jobs.py` reads the id out of the CIF.

**Regression check:** after all three fixes the 12 UBA1 production jobs rebuild
**byte-identical**, so the running jobs are unaffected.

## What needs rerunning

`jobs_ub75_v2/` — 2 charged jobs × 25 seeds. Not yet submitted; the four UBA1 batches
have the GPU allocation. The `free` jobs need no rerun.

Until then, treat any Ub(1-75)-vs-Ub(1-76) *geometry* from `jobs_ub75/` as invalid. Reach
distances measured far from the reactive centre are less affected, but the thioester
geometry and anything within a bond or two of the isopeptide or thioester carbon is not
usable — and per `STUB_RESULTS.md`, AF3's distortion reaches one bond further than the
atoms it was asked to bond.
