# Does UBE2W refuse Ub(1-75)? Three of four variants in

378 AF3 models, 126 per variant, 25 seeds each. The fourth (non-covalent Ub 1-75) is
still in its MSA search — and it is the one that matters most.

## The charged comparison: no evidence of a UBE2W-side block

| variant | thio min | thio med | formed | clash med | pLDDT tail | iptm |
|---|---|---|---|---|---|---|
| Ub(1-76) charged | 1.59 Å | 1.73 Å | **126/126** | 5.0 | 84.2 | 0.640 |
| Ub(1-75) charged | 1.60 Å | 1.68 Å | **126/126** | 4.0 | 83.9 | 0.610 |
| Ub(1-76) non-covalent | 3.23 Å | 3.48 Å | 0/126 | 2.0 | 74.1 | 0.800 |

AF3 accommodates the shorter donor **just as readily**: every model formed the
thioester, the median distance is marginally shorter, and the short variant has
*fewer* clashes (4 vs 5) — the opposite of a steric block. iptm 0.610 vs 0.640 is not
a meaningful difference.

**But the bond was IMPOSED in both charged jobs.** AF3 was told to bond Gly75 (or
Gly76) to Cys91 and complied 126/126 times. A declared bond forming is not evidence
that the geometry is accessible, so this comparison cannot on its own distinguish a
tolerated donor from an intolerable one.

The non-covalent row is the informative one, and it is sobering: with no bond imposed,
even the **wild type** reaches only 3.23 Å at closest and 0/126 models are at bonding
distance. AF3 does not spontaneously produce the charged geometry for either variant.

## Contact analysis: a finding I nearly reported, and the control that killed it

The differential contact list looked dramatic — Pro82 and Tyr85 at **100% → 0%**.
Before reporting it, I asked which *ubiquitin* residue makes each contact:

| UBE2W residue | Ub76 → Ub75 | contacted by | verdict |
|---|---|---|---|
| Pro82 | 100% → 0% | Gly76 ligand only | **trivial** — the residue is absent |
| Tyr85 | 100% → 0% | Gly76 ligand only | **trivial** |
| Trp144 | 100% → 10% | Gly76 ligand only | **trivial** |
| Ile94 | 88% → 4% | **Leu73**, present in both | real |
| Ser118 | 0% → 52% | **Leu73**, present in both | real |

Three of the five largest changes are bookkeeping: those contacts are made
*exclusively* by the Gly76 ligand, which does not exist in the shorter variant. Its
contacts vanish because the atom vanishes.

**What survives the control is small but real.** Leu73 exists in both variants and it
*moves*: in the wild type it contacts Ile94 (68 counts / 30 models), in the short
variant that falls to 17 and it instead contacts Ser118 (17 counts, 0 in wild type).
The tail slides along the groove rather than being excluded from it.

Mean contacting UBE2W residues per model drops 20.6 → 15.8 (−23%). Part of that is
again the missing residue's own contacts, so the honest statement is *reduced
engagement, largely but not entirely explained by one fewer residue*.

## What this means for the engineering question

On this evidence, **AF3 offers no structural explanation for why Ub(1-75) failed
experimentally**. UBE2W's active site does not refuse the shorter donor: no clash
increase, no confidence collapse, no specific residue that blocks it.

That does **not** implicate UBA1 by elimination. It means this assay cannot
distinguish the two, for a specific reason: AF3 forms whatever bond it is told to
form. Two things would change the picture —

1. **the non-covalent Ub(1-75) run** (in progress): if the shorter tail cannot reach
   Cys91 where the wild type can, that is a UBE2W-side answer;
2. **MD on the charged Ub(1-75) state**: the thioester machinery now works
   (`md/CHARGED_MD_FIRST_RESULT.md`), so whether the shorter donor *holds* attack
   geometry is directly testable and is a dynamic question AF3 cannot address.

One prior worth stating: reduced engagement of the kind seen here (Leu73 sliding,
23% fewer contacts) is the sort of effect that a force field would either amplify into
a dissociation or wash out entirely. Which of those happens is the measurement.

## Numbers on which no interpretation rests

`pLDDT` and `iptm` are reported for model quality only. AF3 confidence has been
anti-predictive for reactivity across four rounds in this project, so the 0.800 iptm
of the non-covalent job — the highest of the three — should not be read as it being
the most "correct" model.
