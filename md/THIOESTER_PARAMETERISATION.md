# Parameterising the Ub–UBE2W thioester for Amber MD

**Question for you:** we need a covalent thioester between ubiquitin Gly76 and the
UBE2W catalytic Cys91 in an Amber ff19SB/OPC simulation. Have you done this, and is
the approach below the one you'd use?

Short version: our first attempt produced a residue that Amber accepted and would
have silently simulated as **two unconnected fragments**. We've since found the
likely cause and a candidate fix, but would rather have your opinion than run
another 13 GPU-hours on a guess.

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

## Our candidate approach — does this match yours?

1. Re-run `antechamber` with **`-at gaff2`** on the model compound so the thioester
   atoms get `c` / `ss` types, then `parmchk2 -s gaff2`. Expect 0 ATTN.
2. Keep the **backbone** at ff19SB protein types and let only the side chain carry
   GAFF2 types — the usual mixed-type modified-residue construction.
3. `prepgen` to excise the caps and force integer charge (this worked fine already).
4. Gate on ATTN before running anything: we now refuse any residue with a
   zero-force-constant term at the reactive centre (`check_frcmod_attn.py`).

**Specific things we'd value your view on:**

* **Mixing ff19SB and GAFF2 types in one residue** — is that acceptable here, or do
  you use a different convention for the acyl-enzyme linkage?
* **Charges.** AM1-BCC on a capped model compound, or RESP from a HF/6-31G* ESP?
  For a charge-separated thioester near a catalytic site we suspect AM1-BCC may be
  too crude, but don't want to over-engineer.
* **Ub as a separate chain.** Gly76 is ubiquitin's C-terminus and Cys91 is in the
  middle of UBE2W, so the bond is inter-chain. Do you build one custom residue
  spanning both (and lose tleap's ability to treat them as separate chains), or
  create the bond in tleap with an explicit `bond` command after loading both, using
  a modified CYS with the acyl group?
* **Whether a restraint is defensible instead.** If proper parameters are a lot of
  work, would you accept an `nmropt` distance restraint holding Gly76 C to Cys91 Sγ
  at 1.81 Å? It keeps ubiquitin's bulk in the active site — the part we actually
  need — but doesn't constrain the angles or the charge distribution.

## Context on why the thioester is the whole point

We've now run four rounds of AlphaFold3 predictions on this. The static geometry is
**anti-predictive**: across all 8 SUMO2 lysines with 100 models each, the one site
that works experimentally (K11) ranked 7th of 8 on how close the nucleophile gets to
the catalytic cysteine, and AUC on the two sites with known outcomes was 0.39 — below
chance.

We also asked AF3 to model the charged state directly. It has no thioester
chemistry: the acyl carbon came out at a bond-angle sum of **338.6 ± 1.7°**, between
planar sp² (360°) and tetrahedral sp³ (328.5°), despite correct sp2 hybridisation in
the input. That's why we want a force field — it has an explicit 6.2 kcal/mol
barrier holding the thioester planar, which AF3 does not.

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

Only the thioester is blocking.
