# Published thioester parameters: what I need from you

Searched Europe PMC, CrossRef and Semantic Scholar across eight query families
(thioester force fields, S-palmitoylation, acyl-CoA, acyl-enzyme intermediates,
E2~Ub simulations, HECT/RBR transfer, acyl carrier proteins, modified-residue
parameterisation). 453 unique records screened on title and abstract.

**Two papers are paywalled and I need them. Everything else I could already read.**

---

## Need from you

### 1. Oda, Fukuyoshi, Nakagaki & Takahashi (2013) — the target paper

*Determination of AMBER Force Field Parameters for Thioester by Quantum Chemical
Calculations*
**Chemistry Letters** 42(11):1206–1208 · doi:[10.1246/cl.130517](https://doi.org/10.1246/cl.130517)

This is exactly our problem, solved and published. From the abstract: the AMBER
parameters "around the sulfur atom of the thioester moiety" were determined by
high-accuracy QM, and atomic charges were computed **for acetylcysteine** — an acyl
group on a cysteine sulfur, which is our linkage. They then validated on short
peptides containing acetylcysteine, QM against MM.

**What I need from it:** the table of bond, angle and dihedral terms around the
sulfur, and the acetylcysteine partial charges. That would replace my
GAFF2-to-protein-type transfer with published values derived for protein atom types.

Only 5 citations, which is why it did not surface earlier — small journal, 3-page
paper, and the title says "thioester" rather than anything about proteins.

### 2. Atsmon-Raz & Tieleman (2017) — palmitoylated cysteine

*Parameterization of Palmitoylated Cysteine, Farnesylated Cysteine,
Geranylgeranylated Cysteine, and Multiple Lipidated Cysteines*
**J. Phys. Chem. B** 121(49):11132–11143 ·
doi:[10.1021/acs.jpcb.7b10175](https://doi.org/10.1021/acs.jpcb.7b10175)

S-palmitoylation **is** a cysteine thioester — the same C(=O)–S linkage as ours with
a lipid tail instead of ubiquitin. 43 citations, so its parameters are the de facto
standard for acylated cysteine.

**Caveat I could not resolve without the text:** the headline force field is Martini
(coarse-grained), which cannot drive our atomistic MD. Such papers usually report
the all-atom reference they mapped from. **What I need:** whether it contains
all-atom (CHARMM36 or AMBER) palmitoyl-Cys parameters, and the S–C(=O) terms if so.
If it is Martini-only, it is not useful to us and we drop it.

---

## Already have — and two of them matter

### 3. Zhao, Schaub, Tsai & Luo (2021) — independent validation of the GAFF2 route

*Development of a Pantetheine Force Field Library for Molecular Modeling*
**J. Chem. Inf. Model.** 61(2):856–868 · doi:10.1021/acs.jcim.0c01384 · PMC8266206

Builds acyl-CoA and phosphopantetheine **thioester** parameters on **gaff2 + ff14SB**
— the same combination I arrived at — and validates an *"acetyl cysteamine fragment
… as a representative thioester extending unit"* against MP2/aug-cc-pVDZ,
reproducing C–S stretches at 645 and 749 cm⁻¹, the O–C–S bend at 439 cm⁻¹ and the
thioester carbonyl at 1720 cm⁻¹.

**Why this matters:** it is published evidence that GAFF2 describes a thioester
correctly at QM level. If Oda 2013 turns out to be hard to apply, this justifies the
GAFF2 transfer on its own.

### 4. Elftmaoui & Bignon (2023) — methodology template only

*Robust AMBER Force Field Parameters for Glutathionylated Cysteines*
**Int. J. Mol. Sci.** 24(19):15022 · doi:10.3390/ijms241915022 · open access

**Not our chemistry** — glutathionylation is a disulfide (S–S), not a thioester
(C–S), so its numbers are not transferable. Useful as a worked example of
parameterising a modified cysteine for AMBER and validating against experiment
across 33 µs of MD, i.e. the standard to hold ourselves to.

### 5. Parkin transthiolation intermediate (2024) — what the field actually does

*Capturing the catalytic intermediates of parkin ubiquitination*
**PNAS** 121:e2403114121 · doi:10.1073/pnas.2403114121 · open access

The closest published precedent to our system: the E2~Ub → E3~Ub transthiolation
intermediate, modelled and validated. **They do not parameterise a thioester at
all.** The covalent linkage is imposed as *"a single unambiguous restraint … between
the catalytic cysteine"* in HADDOCK, and the MD is plain **AMBER99SB** with no
custom residue.

Worth knowing before we spend effort: the leading structural paper on exactly our
intermediate sidestepped the parameter problem. That does not mean we should — a
restraint fixes the distance but not the sp² geometry or the charge distribution,
and the planarity is part of what we want to measure — but it means a restrained run
is defensible, published practice rather than a shortcut.

---

## Where this leaves the parameters

I do not think we are blocked either way:

* **GAFF2 already contains every thioester term** with real force constants
  (`c–ss` 199.66 kcal/mol/Å², 1.8104 Å; `o–c–ss` 123.38°; `X–c–ss–X` 6.2 kcal/mol),
  and I have transferred them onto protein types and confirmed the bond is live in
  the topology (`C–SG` at k=199.66, r_eq=1.8104 Å, with the sulfur carrying both its
  bonds).
* Zhao 2021 validates that GAFF2 thioester description against MP2 QM.
* Oda 2013 would let us cite parameters derived *for protein atom types on
  acetylcysteine* rather than justify a transfer, which is strictly better.

So Oda 2013 is an upgrade in defensibility, not a prerequisite. If it turns out to
disagree with the GAFF2 values, that disagreement is itself worth knowing before we
commit GPU time.
