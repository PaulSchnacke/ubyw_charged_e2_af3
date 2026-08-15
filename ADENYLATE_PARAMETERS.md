# Ub–AMP adenylate parameters: the papers, and the check that changes the plan

Supersedes the "what I need from you" section of `ADENYLATE_LITERATURE.md`. Paul supplied
Hegazy & Richards 2013 (`doi:10.1007/s00894-013-1990-x`), **its supporting information**,
and Vanommeslaeghe 2009 (the CGenFF paper). The SI is the valuable part — it contains the
complete optimised parameter set, not just the derivation.

**Headline: GAFF2 already contains every acyl-phosphate junction term.** The same
thing that happened with the thioester has happened again, and the paper's main value is
now as *independent validation* rather than as a source to translate.

## What the paper actually provides

Hegazy & Richards parameterised **acetyl phosphate** (their model compound 4) — which is
our chemistry exactly, minus the adenosine. Their method, from the main text:

* geometries and dihedral scans at **MP2/6-31+G(d)**, single points at MP2/cc-pVTZ;
* **charges fitted to water-interaction energies**, not to an ESP — the CGenFF
  convention. Water was placed at each donor/acceptor in TIP3P geometry and each
  interaction distance optimised at HF/6-31G(d);
* interaction energies **unscaled**, because both model compounds are anionic;
* they introduced one new atom type, `OG305` (the bridging ester oxygen), borrowing its
  Lennard-Jones parameters from `OG303`.

The SI gives the complete CHARMM stream file for `RESI ampa` (total charge −1):

| atom | type | charge | | atom | type | charge |
|---|---|---|---|---|---|---|
| P1 | PG1 | +1.30 | | O8 (bridge) | **OG305** | −0.38 |
| O2, O7 | OG2P1 | −0.71 | | C9 (acyl C) | CG2O2 | +0.34 |
| O3 | OG303 | −0.46 | | O10 (carbonyl O) | OG2D1 | −0.48 |

and the new bonded terms (their Table 5).

## The check I said should come first — and its result

I flagged in `ADENYLATE_LITERATURE.md` that the right first move was to ask whether
GAFF2 already has these terms, exactly as it turned out to for the thioester. Run against
`/cluster/software/manual/amber/24/x86_64/dat/leap/parm/gaff2.dat`, with the acyl
phosphate mapped to GAFF2 types `c` (acyl C) – `os` (bridging O) – `p5` (phosphorus):

**Bonds — all present.**

```
c -os   351.31  1.3621          os-p5   387.59  1.6182
c -o    590.60  1.2190          o -p5   573.03  1.4885
c -c3   234.85  1.5270
```

**Angles — all present except one.**

```
c -os-p5    56.89  122.16       os-p5-os  102.92  101.72
o -c -os   101.26  123.20       o -p5-os  100.30  115.44
c3-c -os    77.97  110.80
os-p5-o     ** MISSING **   <- see below
```

**Dihedrals — present, and specifically parameterised rather than generic.** GAFF2 has
not only the generic `X-c-os-X` and `X-os-p5-X` wildcards but the *specific* junction
torsions:

```
c3-c -os-p5   1   2.700  180.0  -2.0     Junmei
c -os-p5-o    1   0.800    0.0  -2.0     Junmei
o -c -os-c3   1   2.700  180.0  -2.0     Junmei
c3-os-p5-os   1   3.520    0.0   2.0
```

`c3-c-os-p5` and `c-os-p5-o` are the two rotatable torsions that define acyl-phosphate
conformation — precisely the ones Hegazy & Richards scanned at MP2 and fitted.

The one apparent gap, `os-p5-o`, is **not a real gap**: it is the same angle as
`o -p5-os`, which is present at `100.30 / 115.44`. Amber angle types are written in one
canonical order and my probe searched both spellings; only one can match. So the
practical answer is that **nothing is missing**.

## What this changes

My earlier framing — that this was "plausibly a harder problem than the thioester" — was
wrong in the same direction as the original claim it was correcting. The honest summary:

* **No new derivation is needed for the bonded terms.** GAFF2 covers the junction, with
  specific (not wildcard) torsions at the two rotatable bonds.
* **The paper becomes a validation target rather than a source.** GAFF2 and CGenFF
  disagree in ways worth checking, since the two force fields were fitted to different
  observables (GAFF2 to QM energies; CGenFF partly to water-interaction energies):

  | term | GAFF2 | CGenFF (Hegazy) |
  |---|---|---|
  | `os–p5` / `PG1–OG305` bond | k 387.6, r 1.618 Å | k 170, r **1.78 Å** |
  | `c–os` / `CG2O2–OG305` bond | k 351.3, r 1.362 Å | k 230, r 1.34 Å |
  | `c–os–p5` angle | k 56.9, **122.2°** | k 70, **121.5°** |
  | `o–c–os` angle | k 101.3, 123.2° | k 70, 118.0° |

  The **angles agree well** (122.2 vs 121.5°), which is the reassuring part — that is the
  junction geometry. The `os–p5` **bond length differs by 0.16 Å**, which is large enough
  to matter: GAFF2's 1.618 Å is a typical phosphate ester P–O, while 1.78 Å is long.

  **Resolved against 4NNJ (2.4 Å, Uba1 with ubiquitin-AMP).** Measured directly from the
  deposited coordinates, both AMP copies:

  ```
  copy      P–O5' (ester)   P–O1P   P–O2P   P–O3P
  B:101         1.589       1.515   1.481   1.503
  D:101         1.586       1.512   1.468   1.502
  ```

  The bridging phosphate-ester P–O is **1.586–1.589 Å**, against GAFF2's 1.618 Å and
  CGenFF's 1.78 Å. **GAFF2 is right and the CGenFF value is long by ~0.19 Å** for this
  chemistry. Use GAFF2's bonded terms; do not port the CGenFF bond length across.

  One caveat on that measurement: in 4NNJ the adenylate is deposited as AMP and ubiquitin
  **with no covalent linkage modelled** — the Gly76 carbonyl carbon sits 2.61–2.66 Å from
  P, and 1.49–1.50 Å from O1P, i.e. the acyl–phosphate bond is present in the density but
  not built into the deposited model. So the number above is the *ribose* ester P–O,
  which is the correct like-for-like proxy for a phosphate-ester bond but not literally
  the acyl–O–P bond. The conclusion holds either way, since 1.78 Å is far outside the
  range of any P–O ester bond in the structure.
* **Charges still need deriving**, and this is now the *only* real work. Hegazy's charges
  are CGenFF-convention (fitted to water-interaction energies, `OG305` borrowed LJ from
  `OG303`) and are **not portable to Amber**, which expects RESP at HF/6-31G* —
  deliberately over-polarised to mimic condensed phase. The same reasoning that made Oda
  2013 the right choice for the thioester applies: the junction charges *are* the
  measurement, so AM1-BCC is too crude.

  Their charges remain useful as a **sanity check on sign and magnitude** after our RESP
  fit — if our bridging oxygen comes out at −0.05 or −0.9 rather than in the −0.4 to −0.5
  region, something is wrong.

## Revised plan

1. **Skip the bonded-term derivation.** Confirm with `parmchk2` on the actual Ub–AMP
   residue that no term falls back to a zero force constant, then gate with
   `check_frcmod_attn.py` exactly as for the thioester — that gate catches the failure
   where a missing term at the reactive centre silently simulates the adenylate as two
   unconnected fragments while writing a valid-looking trajectory.
2. ~~Resolve the `os–p5` bond length against 4NNJ.~~ **Done — GAFF2 is correct**
   (1.618 Å vs 1.586–1.589 Å measured; CGenFF's 1.78 Å is long by ~0.19 Å).
3. **Derive RESP charges** at HF/6-31G* on a capped model compound (acetyl-AMP, or acetyl
   phosphate to match theirs directly), then compare against Table 1 for sign and
   magnitude.
4. **Decide the Mg²⁺ treatment.** This is the part with no help from the paper — their
   model compounds are bare anions, whereas our adenylate is Mg-coordinated. Which Mg
   model and how many first-shell waters becomes part of the answer, and it is likely a
   larger source of uncertainty than any junction term.
5. **Run the minimisation comparison** exactly as was done for the thioester: GAFF2
   bonded terms vs the CGenFF-derived values converted to Amber form, and compare
   residual angle strain. That test settled the thioester question and is the right way
   to settle this one.

## Note on the CGenFF paper

Vanommeslaeghe 2009 is the methodology reference for *how* CGenFF parameters are derived
(the penalty scoring, the analogy-based assignment, the water-interaction charge
protocol). It is worth citing for why Hegazy's charges are not Amber-portable, but it
contains nothing to transfer directly. It does not need reading in depth for this task.
