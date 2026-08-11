# The real blocker: a covalent bond between two protein chains

Paul's read was right on both counts. Modelling a thioester **is** solved, and the
colleague who did lysine-linked ubiquitin chains **is** the right person to ask —
but for the topology, not the chemistry. Here is the split, established by reading
the force-field files rather than by guessing.

## Why the isopeptide case is easier than ours

| link | bond term needed | present in `parm19.dat`? |
|---|---|---|
| **Isopeptide** (Ub–lysine chains) | `C–N` | **yes** — 490.0 kcal/mol/Å², 1.335 Å, plus `O-C-N`, `CT-C-N`, `C-N-CT` |
| **Thioester** (Ub–Cys91, ours) | `C–S` | **no — zero entries** |

An isopeptide bond is an ordinary amide, so every term is already standard protein
chemistry. Our thioester needs `C–S`, and no standard amino acid carries one, so the
protein force field never defined it. `parm19` has `CT-S` (Met/Cys thioether) and
`CT-S-CT`, but nothing joining a **carbonyl** carbon to sulfur.

Same topology, one extra ingredient. That is why his method transfers and his
parameter set does not.

## What to ask him — the topology, which is the harder half

He solved exactly our structural problem: a covalent bond between two separate
protein chains, with the reactive atom mid-chain on one side (Cys91 in UBE2W) and a
terminus on the other (Ub Gly76).

1. **One custom residue spanning both chains, or an explicit `bond` in tleap** after
   loading both? We are testing the latter.
2. **Residue numbering after linking** — did the chains stay separately addressable
   for analysis, or did they merge into one unit?
3. **Pre-stripped residues** — `CYX`-style from the standard library, or did he edit
   the library?
4. **Anything that silently broke.** Given this project's history, the highest-value
   question by far.

## Chemistry: settled, no derivation needed

GAFF2 contains every thioester term with real force constants, so
`make_thioester_link_frcmod.py` transfers them onto the protein types they
correspond to (`c→C`, `ss→S`, `c3→CT/XC/2C/3C`, `o→O`). Values are extracted from
`gaff2.dat` at runtime, so there is no transcription error, and each is checked
against `parm19.dat` first so only genuinely missing terms are written.

```
BOND   C -S    199.66    1.8104
ANGLE  O -C -S    76.330   123.380      <- sp2 carbonyl
       CT-C -S    62.840   113.460      (also XC, 2C, 3C)
       C -S -CT   93.670    99.120      (also XC, 2C, 3C)
DIHE   X -C -S -X   2   6.200  180.000  2.0   <- holds the thioester planar
```

**Result: the frcmod loads cleanly and tleap builds the inter-chain thioester with
`Total unperturbed charge: 0.000000` and zero missing-parameter complaints.**

### Two format bugs, both mine, and both misleading

Amber rejects an entire frcmod with `Could not load parameter set` — which reads
like a missing-parameter problem — when the *formatting* is wrong:

1. A **multi-line title**. Exactly one title line is allowed.
2. **No `MASS` section**. It must be present even when empty. Verified against files
   tleap does accept (our own `lyq.frcmod` and the shipped `frcmod.ff19SB`): both
   have `MASS` on line 2.

The first test also caught that ff19SB does **not** type every sp3 carbon `CT` —
tleap asked for `XC-C-S` and `C-S-2C`, because ff19SB gives the alpha carbon type
`XC` and beta carbons `2C`. The generator now emits the angle for each.

## The reassuring part: this failure mode is loud

Without the `C–S` term, tleap **refuses**:

```
Could not find bond parameter for atom types: C - S
Could not find angle parameter for atom types: XC - C - S
Exiting LEaP: Errors = 7
```

No topology is written. Unlike AF3 silently discarding a polymer–polymer bond, and
unlike `parmchk2` quietly emitting zero force constants, a missing inter-chain bond
parameter **cannot reach a trajectory**. That is worth knowing for its own sake.

## CONFIRMED: the bond is real in the topology

tleap accepting a topology is not proof the bond does anything — that distinction
has cost this project real time. So rather than trust the exit status, the bond was
read straight out of `link.parm7`:

```
bonds involving SG (atom 46):
  C(21)  - SG(46)    k = 199.66 kcal/mol/Å²   r_eq = 1.8104 Å   <- the thioester
  SG(46) - CB(43)    k = 227.00 kcal/mol/Å²   r_eq = 1.8100 Å   <- Cys side chain
```

The sulfur carries **both** its bonds — to the cysteine Cβ and to the acyl carbon —
which is precisely a thioester. The force constant is non-zero, where the earlier
failure mode showed `k = 0.00`.

**What this does and does not establish.** It proves the transfer worked
mechanically: the GAFF2 value reached the topology and applies to the intended atom
pair across two chains. It does not independently validate 199.66 as the correct
force constant for this chemistry — that is GAFF2's published `c–ss` value, so the
result is exactly as defensible as GAFF2 is, which for a thioester is the standard
choice.

A minimisation test is also queued (chains translated **25 Å apart**, which a
genuine term must pull back to ~1.81 Å) as dynamic confirmation, but the topology
read already answers the question.

A first minimisation attempt returned `NaN`, traced to **the test geometry, not the
parameters**: tleap's `combine` stacked the two fragments on top of each other
(`1-4 VDW = 6267` at step 1, `NaN` in `VDWAALS`/`EEL`, while `BOND`, `ANGLE` and
`DIHED` were all finite and sane). Translating one fragment away first removes the
overlap.

## Bottom line

The inter-chain thioester is **buildable now**. Reading the topology should have been
the first check rather than the last — it is instant, needs no queue, and answers the
question that tleap's exit status does not.
