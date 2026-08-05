# AF3 3.0.1 silently discards polymer–polymer covalent bonds

Found by a stub-MSA topology test before spending 4–6 h on database searches.
This is the single most important thing in this repo, because it invalidates the
obvious way to build a charged E2 and it fails **without an error**.

## What happens

Given a `bondedAtomPairs` entry between two protein chains — ubiquitin Gly76 `C`
and UBE2W Cys91 `SG` — AF3 3.0.1 logs:

```
structure_cleaning.py:253] Reducing number of bonds for <job> from 2 to 1,
  of which 1 are polymer-polymer bonds and 0 are bad bonds.
```

and then **completes normally, exit 0**, writing a model in which those atoms are
**26.25 Å apart**. There is no error, no warning at the top level, and the output
looks like a successful run. Note the wording: the bond is not called "bad" — it
is dropped for being polymer–polymer.

Measured in the two test models:

| bond | type | result |
|---|---|---|
| substrate Lys12 `NZ` → LisoK `C01` | polymer → ligand | **1.30 Å, formed** |
| Ub Gly76 `C` → Cys91 `SG` | polymer → polymer | **26.25 Å, discarded** |
| Ub Gly76 `C` → LisoK `N01` | polymer → ligand | **1.20 Å, formed** |

So the product state (ubiquitin onto the substrate amine, a ligand) works, and
the pre-transfer thioester state (ubiquitin onto the enzyme cysteine, a polymer)
does not.

**It is not partner-order dependent.** Writing the bond as
`[["B",91,"SG"], ["U",76,"C"]]` instead of `[["U",76,"C"], ["B",91,"SG"]]` is
reduced identically. Both orders were tested.

## Why this bites this project specifically

The value of the collaborator's ligand trick is that **ubiquitin enters as a
protein chain**, so it keeps its MSA, its templates and its fold. That is exactly
what makes the thioester impossible: a bond to a protein chain is by definition
polymer–polymer. The trick's strength and this limitation are the same property.

## Workaround that does work: split the reactive C-terminus off

Model ubiquitin as residues **1–74 as a protein chain**, and its last two
residues (Gly75-Gly76) as **two `GLY` ligand copies**, then rebuild the chain with
bonds. Every bond is then polymer↔ligand or ligand↔ligand, and none is
polymer–polymer:

| bond | type |
|---|---|
| `A:12:NZ` → `L:1:C01` | XisoK acyl onto substrate lysine |
| `U:74:C` → `G:1:N` | tether Gly75 back onto Ub 1–74 |
| `G:1:C` → `H:1:N` | Gly75–Gly76 peptide bond |
| `H:1:C` → `B:91:SG` | **the thioester** |

Tested result: **no bond reduction, all four bonds formed.**

```
XisoK      A/K12/NZ  – L/C01     1.20 Å
tether     U/74/C    – G/N       1.24 Å
Gly75-76   G/C       – H/N       1.34 Å
THIOESTER  H/C       – B/C91/SG  1.79 Å   <- a real C–S bond length
```

And the reconstructed chain is geometrically continuous — Cα–Cα spacings of
3.77 Å (U73–U74), 3.94 Å (U74–Gly75) and 3.61 Å (Gly75–Gly76), indistinguishable
from genuine peptide geometry. The splice leaves no visible break.

## What the workaround costs — read before using it

The bond topology is correct. Whether the **prediction** is trustworthy is a
separate question, and the honest answer is that this has not been established.

1. **The last two residues lose their polymer context.** Gly75-Gly76 become
   free-floating ligands held only by bonds. AF3's ligand handling does not carry
   the sequence-based signal a polymer residue gets, and this is precisely the
   region whose placement the experiment measures. That is an uncomfortable place
   to weaken the model.
2. **The C-terminal tail loses MSA coverage at exactly the wrong residues.**
   Ubiquitin's C-terminal `LRLRGG` is the most functionally constrained part of
   the molecule and the alignment there is deeply informative; truncating at 74
   discards the last two columns of it.
3. **AF3 has no thioester chemistry regardless.** `bondedAtomPairs` enforces a
   bond *length*; the surrounding geometry is not parameterised as a thioester,
   so the carbonyl need not be planar or correctly oriented. Any charged-E2 model
   here is a geometric proxy.
4. **The stub-MSA test says nothing about placement.** These runs used fake
   single-sequence alignments, so the 30.8 Å amine-to-carbonyl distance in the
   test model is meaningless. Only the bond topology was under test.

## Consequence for the handover

The four original jobs in `jobs/` are **not all runnable as intended**:

* the two `_product` jobs are **fine** — both bonds are polymer→ligand
* the two `_charged` jobs will **silently lose the thioester** and produce a
  model with ubiquitin floating free, which would look like a clean negative
  result and mean nothing

Anyone running the charged jobs must check for the `Reducing number of bonds` line
in the log and confirm the thioester distance in the output. `qc_charged.py`
reports `thioester_len` and flags it precisely for this reason.

## How to detect it in any future job

```bash
grep -i 'Reducing number of bonds' <af3_log>     # non-empty = bonds were dropped
```

and always verify the geometry rather than trusting the run's exit code.
