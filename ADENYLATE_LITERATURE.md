# Ub–AMP adenylate parameters for MD: what I found, and what I need from you

Searched Europe PMC across 16 query families (acyl adenylate / acyl phosphate force
fields, aminoacyl-tRNA synthetase MD, acyl-CoA and NRPS adenylation domains, CGenFF
acylphosphate, luciferase, E1 transthiolation, ATP/phosphate AMBER parameters).

**Not blocking the AF3 jobs** — the adenylation site is modelled non-covalently there,
by design. This is for the MD stage.

## First, a correction to something I told you

I said the Ub–AMP linkage is *"better covered in the literature than the thioester was —
probably a smaller problem."* **The search does not support that**, and I should not have
said it before looking. What I actually found:

* There is **no published AMBER parameter set for an acyl adenylate / acyl phosphate on
  a protein C-terminus** that I could locate — nothing equivalent to Oda 2013 for the
  thioester.
* The one directly on-target hit is **CGenFF, not AMBER** (Hegazy & Richards 2013,
  below), so it needs translating across force fields rather than transferring.
* The aminoacyl-tRNA synthetase MD literature exists but the papers surfacing are
  mostly *applications* (inhibitor docking, specificity prediction) that use adenylate
  analogues without publishing derived parameters.

So this is plausibly a **harder** problem than the thioester, not an easier one — the
thioester at least had GAFF2 containing every term natively. Corrected here rather than
left standing.

## The one paper I most want

**Hegazy L, Richards NG (2013). "Optimized CGenFF force-field parameters for
acylphosphate and N-phosphonosulfonimidoyl functional groups."**
*J Mol Model* 19(11):4907-4914 · doi:[10.1007/s00894-013-1990-x](https://doi.org/10.1007/s00894-013-1990-x)
— paywalled, 2 citations.

This is the closest thing to our chemistry in the literature: an **acylphosphate**
functional group, deliberately parameterised, with the QM target data. Our Ub–AMP is an
acyl phosphate (a mixed carboxylic–phosphoric anhydride).

**What I need from it:** the bond, angle and dihedral terms around the
C(=O)–O–P linkage, the partial charges, and **which QM level the ESP/charges were fitted
at**. It is CGenFF, so the charges are not directly portable to AMBER (different
scaling and charge philosophy) — but the *bonded* terms and especially the equilibrium
geometries and the QM reference data are what I would use to build and then validate an
AMBER-compatible set. Knowing their QM target also tells me what to fit against.

## Two more worth having if easy

**Truongvan N, et al. (2022) — UBA6 / dual specificity**, *Nat Commun*, PDB `7ZH9`
("Uba1 in complex with ATP", 1.72 Å). Open access, so I can likely reach the text —
what I want is the **ATP-site geometry** as a reference for the reach measurement: the
real distance from the Ub C-terminal carboxyl carbon to ATP Pα in a genuine
adenylation-competent complex. That is the number our AF3 reach values should be judged
against. `4NNJ` (Schäfer et al. 2014, *Acta Cryst D*, "Uba1 loaded with two ubiquitin
molecules", 2.4 Å, contains **AMP**) is arguably better still — it is the doubly-loaded
state with the adenylate site occupied.

**Any recent paper you have on Ubl charging / UBA1 engineered to charge SUMO.** Your
handover lists this as a topic to ask about rather than search, so I have not gone
looking. Two things surfaced anyway that may or may not be what you had in mind:
- *"A modular Uba1-nanobody fusion enables selective ubiquitin transfer to tagged E2
  enzymes"* (2025, *JBC*, doi:10.1016/j.jbc.2025.110910, open access)
- *"Cryo-EM structures reveal the molecular mechanism of SUMO E1-E2 thioester transfer"*
  (2025, *Nat Struct Mol Biol*, doi:10.1038/s41594-025-01681-8)

Tell me which papers you actually meant and I will read those instead of guessing.

## Structures I can already use (no access needed)

From an RCSB search, ranked by usefulness for the reach reference:

| PDB | res | contents | why |
|---|---|---|---|
| **4NNJ** | 2.4 Å | Uba1 + **ubiquitin-AMP** + thioesterified Ub (yeast) | the doubly-loaded state: **both** our sites occupied at once |
| **6DC6** | 3.14 Å | **human** UBA1 + ubiquitin, MG + POP | human sequence, matches our construct |
| **7ZH9** | 1.72 Å | Uba1 + ATP, K + MG | highest resolution ATP-site geometry |
| **9MC4-9MC9** | 3.3 Å | **human** UBA1-UBE2O-Ub, transthiolation states 1-4 | the UBE2W-present question, cryo-EM, several states |
| **3CMM** | 2.7 Å | Uba1-ubiquitin (Lee & Schindelin, *Cell* 2008) | the classic reference |

Note **4NNJ, 7ZH9 and 3CMM are yeast**; 6DC6 and the 9MC* series are human. Our jobs use
human UBA1 (P22314, 1058 aa), so 6DC6 and 9MC4-9MC9 are the sequence-matched comparators.

The 9MC* series is directly relevant to the UBE2W-present arm of this design — it is
human UBA1 with an E2 bound in four resolved transthiolation states. Worth reading
before interpreting our `with_ube2w` jobs.

## My provisional plan for the MD parameters

Not started, and I would rather you sanity-check the route than have me build it:

1. **Split the problem.** The AMP/phosphate half is well covered (nucleotide force
   fields; Meagher et al. 2003 for phosphates). The novel piece is only the
   **C(=O)–O–P** junction — the same shape of problem as the thioester, where only `C–S`
   was missing.
2. **Check GAFF2 first, as we did for the thioester.** GAFF2 turned out to contain every
   thioester term with real force constants, which made that problem evaporate. The
   equivalent check here is whether `gaff2.dat` has the acyl-phosphate ester terms
   (`c–os–p5` and neighbours). This is a 10-minute check on the cluster and it should be
   done before any derivation — I have not run it yet.
3. **Charges will need deriving.** There is no Oda-equivalent published set. RESP on a
   capped model compound (acetyl-AMP or acetyl phosphate) at HF/6-31G*, matching AMBER
   convention, is the defensible route. AM1-BCC would be too crude here for the same
   reason it was too crude for the thioester: the junction charges *are* the
   measurement.
4. **Gate on ATTN exactly as before.** `check_frcmod_attn.py` refuses a residue with a
   zero-force-constant term at the reactive centre — the failure that would otherwise
   simulate the adenylate as two unconnected fragments while writing a valid trajectory.

One thing that argues for doing this properly: the adenylate is **charge-separated and
Mg-coordinated**, so it is a worse case for transferable parameters than the thioester
was, and the Mg²⁺ treatment (which Mg model, how many waters in the first shell) becomes
part of the answer.
