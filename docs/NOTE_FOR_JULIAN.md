# Note for Julian — the four AF3 jobs

Short version: the two `_charged` files needed a non-obvious fix, because **AF3
3.0.1 silently discards covalent bonds between two polymer chains**. The attached
versions work; an earlier draft did not, and failed without an error.

## The limitation

A `bondedAtomPairs` entry between two protein chains (ubiquitin Gly76 `C` to
UBE2W Cys91 `SG`) is dropped by `structure_cleaning.py`:

```
structure_cleaning.py:253] Reducing number of bonds from 2 to 1,
  of which 1 are polymer-polymer bonds and 0 are bad bonds.
```

AF3 then exits **0** and writes a model with those atoms **26.25 Å apart**.
Ubiquitin floats free, the run looks successful, and the result would read as a
clean negative while meaning nothing. Not partner-order dependent — reversing the
pair is reduced identically (both tested).

Awkwardly, this is the ligand trick's own strength turned against it: ubiquitin as
a protein chain is what keeps its MSA and fold, and is exactly what makes a bond
to it polymer–polymer.

## How the attached charged jobs get around it

Ubiquitin is split: **residues 1–74 as a protein chain**, and **Gly75-Gly76 as two
`GLY` ligand copies**, with the chain rebuilt through bonds. Nothing is
polymer–polymer, so all four bonds survive.

| bond | atoms | tested |
|---|---|---|
| XisoK acyl onto substrate Lys | `A:12:NZ` → `L:1:C01` | 1.20 Å |
| Arg74–Gly75 peptide | `U:74:C` → `G:1:N` | 1.24 Å |
| Gly75–Gly76 peptide | `G:1:C` → `H:1:N` | 1.34 Å |
| **thioester** | `H:1:C` → `B:91:SG` | **1.79 Å** |

Verified with a stub-MSA run (`--norun_data_pipeline`, fake single-sequence
alignments) purely to check topology before committing to database searches. The
reconstructed chain is geometrically continuous: Cα–Cα 3.77 Å (U73–U74), 3.94 Å
(U74–Gly75), 3.61 Å (Gly75–Gly76).

The two `_product` jobs needed no change — both their bonds are polymer→ligand.

## Please sanity-check on your first job

```bash
grep -i 'Reducing number of bonds' <log>     # should find NOTHING
```

and confirm in the output that `H:1:C` sits ~1.8 Å from `B:91:SG`. If either
check fails on your AF3 version, the charged runs aren't modelling a loaded E2
and the numbers shouldn't be used. `qc_charged.py` reports `thioester_len` for
this.

## Two caveats we can't engineer away

1. **The split weakens the model exactly where we're measuring.** Gly75-Gly76
   lose polymer context and the last two columns of ubiquitin's MSA — and
   ubiquitin's C-terminal `LRLRGG` is its most constrained region. The topology
   is right; whether the *prediction* is trustworthy there is genuinely open, and
   the product jobs are the control that speaks to it.
2. **AF3 has no thioester chemistry.** The bond length is enforced but the
   surrounding geometry isn't parameterised as a thioester, so the carbonyl need
   not be planar or correctly oriented. This is a geometric proxy for a charged
   E2, not a faithful one.

## What to measure

With ubiquitin loaded the electrophile is the **thioester carbonyl carbon**, not
the cysteine sulfur. So the distance of interest is

```
L:1:N01  ->  H:1:C          (LisoK free α-amine to the thioester carbonyl)
```

`qc_charged.py` also reports the amine-to-Cys91-SG distance, purely so these runs
stay comparable with the earlier uncharged rounds.

## On the seed count

25 seeds is right, and here's the empirical reason. In the uncharged rounds we
used 5 seeds × 5 diffusion samples, and resampling those 300 models shows the
**ensemble minimum** — the statistic we compare against ~4 Å — is not converged at
5 seeds. For the one site that reached attack geometry, per-seed minima were
12.54, 3.47, 4.08, 7.35, 15.32 Å; the run-to-run SD of the minimum is **1.13 Å at
5 seeds** versus **0.07 Å at 20**. A 3.5 Å "result" from 5 seeds is one seed out
of five.

Since the MSA is computed once per sequence set (4–6 h in our hands, I/O-bound on
the database reads) and inference is ~10 min per seed, seeds are the cheap axis.

All four jobs share the same three protein sequences, so if your pipeline can
compute the data pipeline once and reuse it across them, that's most of the cost
saved.

## Full detail

`docs/AF3_POLYMER_BOND_LIMIT.md` in
[ubyw_charged_e2_af3](https://github.com/PaulSchnacke/ubyw_charged_e2_af3),
including the test models and logs under `bondtest/`.
