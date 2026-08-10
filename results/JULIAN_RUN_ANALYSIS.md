# Analysis of Julian's run — two bugs in my ligand chemistry

100 models (4 jobs × 25 seeds). Paul spotted both problems by eye before I
measured anything. **He is right on both counts, and neither was on purpose.**

Summary: the topology is right, the *chemistry* of two carbons is wrong, and the
K11/K21 difference in the charged jobs is an artefact of the second bug rather
than a result.

## Bug 1 — the isopeptide carbon is over-valent (the "cyclopropene")

The LisoK carbonyl carbon `C01` was declared in `LisoK_userCCD.cif` with **four**
internal bonds:

```
O01 C01 DOUB      +2
C01 C02 SING      +1
C01 H01 SING      +1     <-- this one should not exist
```

AF3 then adds the external isopeptide bond to Lys `NZ`, giving `C01` **five**
connections. Carbon cannot do that, so the geometry collapses:

| measured across all 100 models | ideal |
|---|---|
| `NZ–C01` **0.93–1.06 Å** (mean by job) | amide C–N **1.33 Å** |
| `O01–C01` **1.37 Å** | carbonyl C=O **1.23 Å** |

`NZ` at 0.9 Å from `C01` with `O01` also crowded in makes a tight three-atom
cluster, which is exactly what a viewer draws as a cyclopropene ring. **It is not
a deliberate planarity constraint — it is a valence error rendered as a ring.**

Reference values from an MMFF-optimised N-methylacetamide: C=O 1.226 Å,
C–N 1.377 Å.

### Where it came from

I built the ligand from SMILES `O=C[C@@H](N)CC(C)C`, reasoning that OXT should be
absent because it is the leaving group of the amide. Dropping OXT was correct.
What I missed is that this makes the SMILES an **aldehyde**, so RDKit added an
aldehyde hydrogen to complete the carbon's valence — and that H was written into
the CCD. The acyl carbon should have an open valence for the external bond, i.e.
no `H01` atom row and no `C01 H01` bond row.

## Bug 2 — the thioester carbon keeps its OXT, so it is sp3 not planar

Paul: *"in that case the G76 thioester should definitely NOT have planar
geometry"* — and it does not, but for the wrong reason.

The two C-terminal glycines were supplied as standard `GLY` CCD ligands. Standard
`GLY` carries **`OXT`**, the carboxylate's second oxygen. So the carbon that
becomes the thioester already has C, O, OXT, CA — and adding the bond to Cys91
`SG` makes it five-coordinate again.

Measured across all 50 charged models:

| | measured | ideal |
|---|---|---|
| thioester `C–S` | 1.72–1.73 Å | 1.78 Å (acceptable) |
| **planarity sum at C** | **328° (K11), 332° (K21)** | **360° for sp2** |
| OXT still present | **50/50 models** | should be absent |

328° is textbook tetrahedral (sp3 ≈ 328.5°). A real thioester carbonyl is planar
sp2 at 360°. So the carbon is hybridised wrongly, and the cause is the retained
OXT plus the extra bond.

## Bug 3 (consequence) — AF3 built a tetrahedral intermediate, not a charged E2

This is Paul's other suspicion — *"either there is a clash of residues for the
charged module or you tried to attach XisoK covalently to the thioester"* — and
the second reading is effectively what happened, though I never declared such a
bond.

In the **K11 charged** models the LisoK free amine `N01` sits **2.12 Å mean
(1.81–2.94)** from the thioester carbon, with **23/25 models under 2.5 Å**. No
covale record connects them; AF3 simply drove them together because that carbon
had room in its (already broken) valence.

At 2.1 Å this is not attack geometry — it is shorter than the C–N distance of a
genuine transition state. Combined with the sp3 carbon and five substituents,
**AF3 has modelled the tetrahedral intermediate of the transfer reaction**, i.e.
the reaction caught mid-flight. That is a defensible thing to model deliberately,
but it is not what these jobs were meant to be, and it was not chosen.

Genuine clashes in the K11 top model (excluding peptide bonds and the three
declared bonds): **11**, the worst being

```
0.48 A  U/GLY75/OXT  -- UA/GLY1/N
1.61 A  A/LYS12/CE   -- L/LIG-1/C01
1.74 A  A/LYS12/NZ   -- L/LIG-1/O01
1.82 A  L/LIG-1/N01  -- UA/GLY1/OXT
2.08 A  B/CYS91/SG   -- UA/GLY1/OXT
```

The retained OXT atoms are directly involved in most of them.

## The K11 vs K21 difference Paul noticed is NOT a result

| job | amine → thioester C | ≤2.5 Å | >10 Å |
|---|---|---|---|
| K11 charged | mean 2.12 Å | **23/25** | 0/25 |
| K21 charged | mean 29.84 Å | 0/25 | 24/25 |

That looks like a perfect discriminator, and it is not one. The K11 lysine simply
sits closer to the active site in the docked complex (Cα-to-Cys91-SG 8.8 Å versus
21.1 Å for K21) while both substrates dock comparably (9 vs 11 interface residues
within 8 Å). Given a carbon with a spare valence, the nearer amine gets pulled in
and the farther one cannot reach. **The split measures reach, on a broken carbon
— exactly the failure mode rounds 2–5 already established, now with an extra
artefact on top.** It must not be reported as a discriminator.

## What IS sound

* The intended bond **topology** survived. All three covale records are present
  and connect the intended atoms; the split-ubiquitin workaround for the
  polymer–polymer restriction did its job.
* The correct lysine is bonded in every model: residue 12 in the K11 jobs and 22
  in the K21 jobs, 25/25 each.
* **AF3 renumbered the split.** It absorbed the first `GLY` ligand into the
  protein chain (chain `U` is 75 residues, not 74) and kept the second as chain
  `UA`. So the real thioester bond is `B/CYS91/SG – UA/GLY1/C`, not the
  `H:1:C – B:91:SG` I wrote. Any QC script must resolve chains from the covale
  records rather than assume the input naming.
* The **product** jobs have the same isopeptide over-valence (Bug 1) but no
  thioester, so they are less badly broken — Ub Gly76 C to the LisoK amine is
  1.29 Å (K11) and 1.59 Å (K21), still short of the ideal 1.33 but far closer than
  the `NZ–C01` bond.

## Fixes required before rerunning anything

1. **Remove `H01` from the LisoK CCD** — the atom row and the `C01 H01` bond row.
   Rebuild all XisoK variants the same way; `make_xisok_ccd.py` has the same
   defect for every variant.
2. **Do not use standard `GLY` for the reactive C-terminus.** Supply a custom
   two-residue CCD for Gly-Gly with **no OXT** on the terminal carbon, so that
   carbon has exactly C=O, CA and the incoming S.
3. **Decide explicitly which species to model** and enforce it:
   * *charged E2* — thioester intact, substrate amine free. Requires the amine to
     be prevented from bonding the thioester carbon, which AF3 will not guarantee
     while a spare valence exists. Correct valences are the only real defence.
   * *tetrahedral intermediate* — what these runs accidentally produced. Modelling
     it deliberately is legitimate and interesting, but then the carbon should be
     declared sp3 with all four substituents and the geometry checked against a
     real intermediate.
   * *product* — already works.
4. **Add a valence assertion to the CCD builder**: sum the declared bond orders
   per atom, add 1 for any atom named in `bondedAtomPairs`, and refuse to write a
   file where carbon exceeds 4. This would have caught both bugs before any
   compute was spent.

## On n=2

Paul is right that two sites cannot settle anything, and the site sweep is the
agreed next step. But it should not run until the valences are fixed — 79 sites
of a broken carbon would produce a large, self-consistent, meaningless dataset.
