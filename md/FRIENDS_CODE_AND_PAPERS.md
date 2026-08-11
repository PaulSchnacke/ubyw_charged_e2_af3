# Your friend's code, and the two papers

Both papers gave up what we needed, and the code answers the topology question —
though not the way I expected. Summary first, then details.

| source | verdict |
|---|---|
| **Oda 2013** | **Usable and directly on target.** Full parameter table + RESP charges extracted below. |
| **Tieleman 2017** | Not a new parameter set. Its all-atom values *ship in CHARMM36* already. |
| **Friend's code** | Solves our exact topology problem — in **OpenMM/CHARMM36m**, not Amber, with an approach that has one deliberate approximation we cannot accept. |

---

## 1. Oda 2013 — the parameters, extracted and verified

**Table 1, AMBER force field parameters for the thioester moiety** (ff99SB and
ff12SB; the paper states the values are identical in both):

| term | K | equilibrium |
|---|---|---|
| BOND `C–S` | 208.865 kcal mol⁻¹ Å⁻² | 1.769 Å |
| ANGLE `C–S–C` | 30.999 kcal mol⁻¹ rad⁻² | 81.535° |
| ANGLE `C–C–S` | 44.634 | 104.807° |
| ANGLE `O–C–S` | 32.521 | 99.676° |
| DIHEDRAL `X–C–S–X` | Vₙ/2 = 8.941 kcal mol⁻¹ | 180°, n = 2, 2 paths |

**Figure 5, RESP charges for acetylcysteine** (HF/6-31G(df,p); backbone amide
charges set to standard ff99SB values, the rest from the RESP fit):

| atom | charge (e) | | atom | charge (e) |
|---|---|---|---|---|
| **S (thioester)** | **−0.307889** | | CA | +0.024622 |
| **C (acyl carbonyl)** | **+0.479917** | | C (backbone) | +0.5973 |
| **O (acyl carbonyl)** | **−0.414690** | | O (backbone) | −0.5679 |
| CT (acetyl methyl) | −0.156886 | | N | −0.4157 |
| HC (×3) | +0.073065 | | H (amide) | +0.2719 |
| CB | +0.076612 | | HB (×2) | +0.074370 |

Torsion fitted against **CCSD(T)/aug-cc-pVTZ**; validated on Ace–Ala–Cys–Ala–Nme
conformer energies, QM vs MM, with a better r² than GAFF.

**Why the charges matter as much as the bonded terms.** The thioester sulfur at
**−0.308 e** is roughly three times more negative than a standard Cys SG (−0.108 e
in ff19SB), and the acyl carbon at **+0.480 e** is the electrophile our XisoK amine
has to attack. Any approach that leaves standard Cys charges in place gets both
wrong — which is exactly the approximation in the code below.

### Verification

I extracted these from the text layer, then checked them against the paper's own
four comparative statements, all of which hold:

* "req … similar to the GAFF value of req = 1.762 Å" → 1.769 vs 1.762, diff 0.007 Å
* "force constants Kr and K_θ … were smaller than those of GAFF" → 31.0 < 60.9,
  44.6 < 61.5, 32.5 < 63.0
* "the equilibrium bond angles were smaller" → 81.5 < 99.2, 104.8 < 113.5, 99.7 < 123.3
* "Vₙ/2 for X–C–S–X was larger than that of X–c–ss–X in GAFF (6.200)" → 8.941 > 6.200

That last one independently confirms the GAFF value I read off the cluster.

### The caveat, tested and resolved

Three of those equilibrium **angles are not physical as geometries**: a real sp²
thioester carbonyl has O=C–S near 123°, not 99.7°. That pattern suggests *effective*
parameters fitted against energy profiles that still contain the 1–3 nonbonded
repulsion — in which case minimisation still lands at the physical geometry.

**Tested rather than assumed.** Built the same capped model system with each
parameter set (fragments translated 25 Å apart, so a real bond has to pull them
together) and minimised:

| term | physical | Oda 2013 | GAFF2 transfer |
|---|---|---|---|
| C–S | ~1.78 Å | **1.794 Å** | **1.822 Å** |
| O=C–S | ~123° | **119.7°** | **124.4°** |
| C–C–S | ~116° | 114.7° | 114.7° |
| C–S–C | ~100° | 97.4° | 100.9° |

**Both give a correct planar sp² thioester.** Oda's odd equilibrium values are
confirmed effective parameters and are safe to use — the 99.7° equilibrium
minimises to 119.7° once the nonbonded terms act.

But the energy decomposition at the minimum favours GAFF2:

```
Oda     BOND 1.972   ANGLE 11.098   DIHED 6.295
GAFF2   BOND 1.835   ANGLE  2.933   DIHED 6.291
```

GAFF2 sits in a relaxed minimum; Oda carries **~8 kcal/mol more angle strain**,
because its equilibria are far from the geometry the nonbonded terms enforce. That
residual strain is real — it pulls on the surrounding geometry — and it is the
argument for using GAFF2's bonded terms as primary.

**Recommendation: GAFF2 bonded terms + Oda's RESP charges.** The charges are where
Oda is unambiguously better (there is no GAFF2 equivalent — AM1-BCC on a fragment is
not a substitute for a RESP fit on the actual acetylcysteine), and GAFF2's angles
give the relaxed geometry. Oda's stiffer CCSD(T)-fitted torsion (8.941 vs 6.200) is
worth testing as a variant, since planarity is part of what we measure.

**Other caveats, from the paper:** derived against GAFF (not GAFF2) and for
ff99SB/ff12SB (not ff19SB); and the thioester sulfur's **vdW parameters were not
refitted** — they were taken from the existing Met/Cys `S`/`SH` types.

## 2. Tieleman 2017 — resolves my open question, negatively

The headline force field is Martini (coarse-grained), unusable for us. But the open
question was whether it also derived all-atom parameters. It did not: the parameters
were *"obtained directly from the February 2016 GROMACS implementation of the
CHARMM36 force field … available from the Mackerell group website"*, with the
lipidated amino acids *"published explicitly in the CHARMM format"*.

**So no new numbers — but a useful fact:** palmitoyl-cysteine, which *is* a
thioester, is a shipped CHARMM36 residue. That matters for the route below.

---

## 3. Your friend's code — solves the topology, with one approximation we can't take

`src_md` screens AF3 predictions of HDAC6 bound to K48/K63 polyubiquitin under MD.
Different biology, but **structurally our exact problem**: a covalent bond joining
two separate protein chains, parsed out of an AF3 `_struct_conn` record.

### What transfers directly, and is genuinely valuable

**The strip-and-re-add protocol.** Both `addHydrogens` and `createSystem` do residue
template matching that fails when Lys NZ carries an external bond:

> "No template found for residue LYS. The externally bonded atoms has 1 N atom too many."

His solution: remove all inter-chain bonds → add hydrogens → re-add the bonds. Same
again around `addSolvent`. That is the OpenMM analogue of the problem I hit in tleap,
and it's a cleaner answer than mine.

**Reading the linkage from `_struct_conn`.** He parses AF3's covale records to find
the bond, exactly as our QC does — and independently discovered the same AF3
artefacts we did:

* AF3 splits Gly76 into its **own single-residue chain** and encodes Gly75.C→Gly76.N
  as an inter-chain bond. He merges it back into the donor chain. *This is the same
  split-Gly topology I had to invent to dodge AF3's polymer–polymer bond removal —
  he found AF3 doing it natively for isopeptides.*
* Gly75's **OXT has bad coordinates** in AF3 CIFs (0.68–0.77 Å from C). He removes it.
  We hit the same class of problem from the other direction.
* After stripping, Gly76 looks like a free C-terminus so `addHydrogens` adds an
  **OXT that must be deleted**, or the acyl carbon ends up bonded to both OXT and the
  nucleophile.
* `addHydrogens` places **HZ1/HZ2/HZ3 at free-amine geometry**, leaving one H ~0.3–0.5 Å
  from the linker carbon — severe clash, NaN on minimisation. He repositions them
  tetrahedrally rather than deleting (the template requires exactly 3).

Every one of those is a silent-failure trap of the kind that has cost us days. Worth
having regardless of which force field we use.

**The bonded patch.** After `createSystem` he adds the missing terms by hand
(`_ensure_isopeptide_bonded_terms`): one harmonic bond (r₀ = 0.1335 nm,
k = 300 000 kJ mol⁻¹ nm⁻²), three angles (CA–C–N 116°, O–C–N 123°, C–N–CE 121°),
one torsion (180°, n = 1, 12 kJ mol⁻¹), **plus nonbonded exclusions for the new 1–2
and 1–3 pairs**. That exclusion step is the part most people forget.

### The approximation we cannot accept

His own docstring is explicit:

> "Standard CHARMM36m assigns LYS NZ the charge of a protonated amine (+0.26 e on N,
> standard for charged Lys) rather than the neutral amide nitrogen charge of the true
> isopeptide bond … For interface STABILITY SCREENING (does the complex stay bound?)
> this approximation is acceptable. For binding FREE ENERGY calculations … use
> CHARMM-GUI polyubiquitin PATCH parameters on top candidates."

He is asking "does the complex stay bound", so a charge error at the junction is
tolerable. **We are asking whether a nucleophile approaches an electrophile** — the
charges at the junction are the measurement. With standard Cys charges the sulfur
would sit at −0.108 e instead of −0.308 e and the acyl carbon would not carry its
+0.48 e at all.

So we take his topology handling and reject his charge treatment. Oda's RESP charges
are exactly the missing piece — and he tells us where to look if we went the CHARMM
route: CHARMM-GUI patches, or the shipped palmitoyl-Cys residue that Tieleman points at.

### Amber or OpenMM?

Sticking with Amber, because our validated LYQ isopeptide residue, the equilibrated
systems, and the analysis scripts are all Amber. His pipeline is the fallback if the
tleap route fails, and his traps apply either way.

One thing worth stealing outright: his interface observables (interface RMSD +
**native contact fraction Q**) are better than backbone RMSD for "did it stay
engaged", which we learned the hard way in an earlier round.

---

## Status

* **Parameters: settled.** GAFF2 bonded terms (verified to give relaxed, correct sp²
  geometry) plus Oda's RESP charges (the thioester S at −0.308 e and acyl C at
  +0.480 e that no standard residue provides).
* **Topology: solved**, both in tleap (verified: `C–SG` bond live in the topology at
  k = 199.66, r_eq = 1.8104 Å) and now with an independent OpenMM route.
* **Charges: Oda's RESP set** replaces the approximation both my GAFF2 transfer and
  your friend's pipeline would otherwise carry.
