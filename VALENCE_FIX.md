# The valence fix, and how to test formatting before spending compute

Paul asked: *"is there any way to test that the formatting is right this time?"*

Yes, and it is now two layers, because formatting was never the thing that broke.
The first charged run had schema-valid JSON with every bond on the correct residue
type — and produced a pentavalent carbon.

## Layer 1: `ccd_valence.py` — catches the bug class locally, in milliseconds

Sums declared bond orders per atom, **adds the external bonds AF3 will create
from `bondedAtomPairs`**, and refuses anything over the element limit. That last
part is the whole point: each half looked fine in isolation.

Run against the ligand that shipped to Julian:

```
$ python ccd_valence.py LisoK_userCCD.cif C01=1
  C01    C   order    5 / 4 (+1 external)  <-- OVER
```

Bug 1, in one command. It is wired into `make_ccd_v2.py` and `build_jobs_v2.py`,
so no file is written and no job is emitted if a valence is exceeded.

**It found a third bug immediately** — one that shipped in Julian's run and that
nobody had flagged. In the *product* jobs, ubiquitin bonds to the XisoK
alpha-amine, but that amine was a free `-NH2` with three declared bonds, so the
external bond made nitrogen 4-valent. That is why the product isopeptide measured
1.29 Å rather than 1.33. Fixed with `LisoK_openN_userCCD.cif`, where one amine
hydrogen is removed — chemically correct, since in the product that nitrogen
genuinely is a secondary amide.

And it caught a fourth while building the Gly-Gly: the free N-terminus had the
same problem for the bond back to Ub Arg74.

## Layer 2: a stub-MSA run on Euler — catches what only AF3 knows

`--norun_data_pipeline` with single-sequence fake alignments runs inference in
~10 min of GPU instead of 4–6 h of database search. It cannot tell you anything
about *placement*, but it tells you exactly what you need here:

* did AF3 keep every bond, or silently drop some (`grep -i 'Reducing number of
  bonds'`)
* what bond lengths and angles came out
* how did AF3 renumber the chains

This is how the polymer–polymer restriction was found in the first place, on
Euler, before any real compute was spent.

## On "the ligand trick does not work on naive AF3 on Euler"

It does. Euler's AF3 3.0.1 built the split-ubiquitin thioester topology correctly
(job `3f65636f`): no bond reduction, all four bonds formed, thioester at 1.79 Å.
What failed on Euler was the *harvest* — the 34 MB augmented JSON exceeded the
transfer threshold, so the job reported failure after succeeding. The MSAs are
intact and archived in `~/ubyw_charged_msa`.

So this can be validated end-to-end without Julian. His queue is only needed if
we want the full 25-seed production runs faster than Euler's queue allows.

## The two chemistry fixes

**Bug 1 — the acyl carbon must not carry a hydrogen.** Built from the aldehyde
SMILES `O=C[...]`, RDKit added an aldehyde H to complete the valence. Correct to
omit OXT (it is the amide's leaving group); incorrect to leave the H that replaced
it. `strip_acyl_hydrogen()` removes it and leaves the valence open for the
external bond.

**Bug 2 — the Ub C-terminal Gly must not carry OXT.** The standard `GLY` CCD keeps
it, giving the thioester carbon C, O, OXT, CA before the sulfur bond. `UBGG` is a
custom Gly-Gly whose terminal carbon has exactly `=O`, `CA2`, and one open
valence.

## Four species, declared rather than inherited

| job | species | bonds |
|---|---|---|
| `sumo2_k11_xisok_ube2w` | XisoK + UBE2W, no Ub | 1 |
| `sumo2_k11_lyscontrol_ube2w` | unmodified Lys + UBE2W | 0 |
| `sumo2_k11_charged_ube2w` | Ub thioester on Cys91, amine free | 3 |
| `sumo2_k11_tetrahedral_ube2w` | amine added across the thioester | 4 |
| `sumo2_k11_product_ube2w` | Ub isopeptide on the XisoK amine | 2 |

The tetrahedral job is the interesting case. AF3 produced that species *by
accident* last time, on a carbon that was simply over-valent. Declaring it
properly requires the carbon to be **sp3 with an alkoxide oxygen** — a C–O single
bond, not C=O — otherwise it cannot accept both the sulfur and the incoming
nitrogen. That is `UBGT_userCCD.cif`. The valence checker refused the job until
the chemistry was right, which is the behaviour we want: it forced a real decision
instead of letting AF3 improvise one.

## What is still not guaranteed

Correct valences stop the carbon collapsing. They do **not** stop AF3 driving the
free amine toward the thioester carbon in the charged job — that is a placement
question, and if it happens again with correct valences it is a genuine prediction
rather than an artefact. The tetrahedral job exists partly so the two can be
compared: one where the bond is declared, one where it is not.
