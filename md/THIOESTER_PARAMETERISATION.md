# Parameterising the Ub–UBE2W thioester for Amber MD

**Question for you:** we need a covalent thioester between ubiquitin Gly76 and the
UBE2W catalytic Cys91 in an Amber ff19SB/OPC simulation. The mechanical problem is
now solved — but three judgement calls remain that we'd rather not guess at, and
they decide whether the result is defensible. Those are at the end.

Short version of the episode, since the failure mode is worth knowing about: our
first attempt produced a residue that Amber accepted and would have silently
simulated as **two unconnected fragments**. Cause was our own atom-type choice, not
a gap in the force field.

---

## The system

UBE2W transfers ubiquitin to a substrate's N-terminus. We're modelling the charged
E2 — ubiquitin's C-terminal Gly76 joined by a **thioester** to Cys91 Sγ — with a
SUMO2 substrate carrying an XisoK modification whose free α-amine is the incoming
nucleophile:

```
Ub Gly76 C(=O)–S–CH2–Cys91        the thioester we need to hold
                    ↑
        XisoK α-NH2 attacks here
```

Three protein chains plus a modified residue, ~110,000 atoms solvated, ff19SB/OPC,
`pmemd.cuda` at ~111 ns/day on a TITAN RTX.

## What went wrong

We built a capped model compound (ACE–X–NME) for the thioester-bearing cysteine and
ran the standard route: `antechamber -c bcc` → `prepgen` → `parmchk2` → `tleap`.

**Every check passed.** AM1-BCC charges assigned, integer residue charge, tleap
built a test tripeptide and reported `Total unperturbed charge: 0.000000`.

**But `parmchk2` had flagged 6 terms `ATTN, need revision` and set them all to zero:**

```
BOND   c -s        k = 0.000    r_eq = 0.000
ANGLE  o -c -s     k = 0.000
ANGLE  ct-c -s     k = 0.000
ANGLE  c -s -ct    k = 0.000
DIHE   o -c -s -ct k = 0.000
DIHE   ct-c -s -ct k = 0.000
```

All six *are* the thioester. And a zero-force-constant bond in Amber is worse than
a weak one: it contributes zero energy at every separation, **and** the 1–2 pair is
still excluded from nonbonded interactions — so the acyl carbon and the sulfur would
have felt no force of any kind from each other. Ubiquitin would have drifted off the
enzyme while the job exited 0 and wrote a perfectly valid trajectory.

## What we think the mistake was

We ran `antechamber -at amber`, i.e. **protein** atom types. Checking the parameter
files directly:

* `parm19.dat` contains **no `C–S` bond term at all** — no standard amino acid has a
  thioester, so the protein force field never needed one.
* `gaff2.dat` has the full set, with real force constants:

| term | GAFF2 | equilibrium | what we got with `-at amber` |
|---|---|---|---|
| BOND `c–ss` | 199.66 kcal/mol/Å² | 1.8104 Å | 0.000 / 0.000 |
| ANGLE `o–c–ss` | 76.33 | 123.38° | 0.000 |
| ANGLE `c3–c–ss` | 62.84 | 113.46° | 0.000 |
| ANGLE `c–ss–c3` | 93.67 | 99.12° | 0.000 |
| DIHE `c3–c–ss–c3` | 1.550 | 180°, n=2 | 0.000 |
| DIHE `X–c–ss–X` | 6.200 | 180°, n=2 | 0.000 |

The `o–c–ss` angle of 123.38° is the sp2 carbonyl, and the 6.2 kcal/mol barrier
about `c–ss` is what keeps the thioester planar. Those are precisely the two
features we need, and both exist already.

So this looks like our error rather than a genuine gap in the force field.

## Resolved: it was the atom-type flag

Re-running `antechamber` with **`-at gaff2`** instead of `-at amber` fixes it
completely. Confirmed on the cluster:

| | sulfur type | acyl carbon | `C–S` in the parameter set? |
|---|---|---|---|
| `-at amber` | `S` (protein) | `C` | **no** → 6 zero-force-constant terms |
| `-at gaff2` | `ss` (thioether) | `c` | **yes** |

The GAFF2 `frcmod` contains **no sulfur terms at all**, and that is the success
signal rather than a failure: `parmchk2` only emits terms it has to *guess*. `c–ss`,
`o–c–ss`, `c3–c–ss` and `c–ss–c3` are absent from the frcmod precisely because
`gaff2.dat` already defines them —

```
BOND   c -ss     237.9 kcal/mol/Å²   1.8000 Å
ANGLE  o -c -ss   63.0              123.32°     <- sp2 carbonyl
ANGLE  c3-c -ss   61.5              113.51°
ANGLE  c -ss-c3   60.9               99.16°
```

Everything the frcmod *does* contain is nitrogen (`ns`, `n8`) with **penalty score
0.0** — exact analogues of existing amide and amine terms. None is
zero-force-constant. So the thioester needs no new parameters at all.

## Remaining approach — does this match yours?

1. `antechamber -at gaff2` → `parmchk2 -s gaff2` (**done, 0 problematic terms**).
2. Keep the **backbone** at ff19SB protein types and let only the side chain carry
   GAFF2 types — the usual mixed-type modified-residue construction.
3. `prepgen` to excise the caps and force integer charge (this worked fine already).
4. Gate on ATTN before running anything: we now refuse any residue with a
   zero-force-constant term at the reactive centre (`check_frcmod_attn.py`).

**The three things we'd genuinely value your view on** — these are judgement
calls, not blockers:

* **Mixing ff19SB and GAFF2 types in one residue** — is that acceptable here, or do
  you use a different convention for the acyl-enzyme linkage? This is the one we're
  least sure about: the backbone needs ff19SB (CMAP), but the thioester needs GAFF2,
  and they meet two bonds apart.
* **Charges.** AM1-BCC on a capped model compound, or RESP from a HF/6-31G* ESP?
  For a charge-separated thioester near a catalytic site we suspect AM1-BCC may be
  too crude, but don't want to over-engineer.
* **Ub as a separate chain.** Gly76 is ubiquitin's C-terminus and Cys91 is in the
  middle of UBE2W, so the bond is inter-chain. Do you build one custom residue
  spanning both (and lose tleap's ability to treat them as separate chains), or
  create the bond in tleap with an explicit `bond` command after loading both, using
  a modified CYS with the acyl group?
(We've dropped the earlier question about whether a distance restraint would be
acceptable — with real `c–ss` parameters in hand there's no need to approximate.)

## Context on why the thioester is the whole point

We've now run four rounds of AlphaFold3 predictions on this. The static geometry is
**anti-predictive**: across all 8 SUMO2 lysines with 100 models each, the one site
that works experimentally (K11) ranked 7th of 8 on how close the nucleophile gets to
the catalytic cysteine, and AUC on the two sites with known outcomes was 0.39 — below
chance.

We also asked AF3 to model the charged state directly. It has no thioester
chemistry: the acyl carbon came out at a bond-angle sum of **338.6 ± 1.7°**, between
planar sp² (360°) and tetrahedral sp³ (328.5°), despite correct sp2 hybridisation in
the input. That's why we want a force field — GAFF2 has an explicit
barrier about `c–ss` holding the thioester planar and an `o–c–ss` equilibrium of
123.32°, neither of which AF3 respects.

MD without the thioester is pointless: with ubiquitin absent, the active-site groove
the nucleophile has to reach into is empty, so the steric competition that decides
the outcome is missing. We ran three such systems by mistake and stopped them.

## What we already have working

* **LYQ** — the isopeptide-linked lysine (XisoK on the substrate). Validated: integer
  charge, amide half agrees with published ALY (Nε-acetyllysine) from
  `ff19SB_modAA` within 0.25 e, nucleophile N at −0.88 e. Built from a capped model
  compound with a structural L-configuration check.
* **LYP** — the product state (Ub amide-linked to the XisoK amine), **0 ATTN**, since
  both linkages are ordinary amides.
* Two systems solvated and fully equilibrated, ~111 ns/day throughput measured.

The thioester is no longer blocking. What we'd like from you is a sanity check on
the three points above before we spend the GPU time, since getting the charge model
or the type-mixing wrong would produce a trajectory that looks fine and means
nothing — which is exactly the trap we just climbed out of.
