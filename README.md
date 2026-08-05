# Charged-E2 co-folds: UBE2W~Ub with a XisoK substrate

**Handover package.** Four AF3 job JSONs, ready to run, plus the QC script for
the outputs. Nothing here needs the Euler setup — any AF3 3.0.x installation
will do.

## The question

Earlier rounds of this project ([ubyw_reactivity_analysis](https://github.com/PaulSchnacke/ubyw_reactivity_analysis))
asked whether AF3 can predict which XisoK sites UbyW modifies. All of them
measured the LisoK neo-N-terminus against a **bare** catalytic cysteine, and all
of them were negative:

| round | test | result |
|---|---|---|
| 1 | sequence motif around the XisoK | no motif; 84 tests, min *q* = 0.91 |
| 2 | amine reaches Cys91, 14 sites | **AUC 0.367** — below chance |
| 3 | 12 XisoK chemistries vs measured yield | correlation inverted, confounded by side-chain bulk |
| 4 | proteome-wide, 53,448 sites | dominated by a tryptic artefact; small flexibility trend survives |
| 5 | lysine-only control | the **modification drives the docking** (large, significant) |
| MD | 3 × 20 ns on the AF3 pose | relaxes to 7.4 ± 0.3 Å; **0 / 1200 frames ≤ 4 Å** |

Round 5 and the MD together diagnose the problem: AF3 docks the modification
into an *empty* active site regardless of whether the site is reactive, and the
sub-4 Å contacts that looked like attack geometry are not a force-field minimum.
An empty cysteine is too permissive a target.

**Loading ubiquitin changes the question** from "can the amine reach a cysteine"
to "can it reach the carbonyl it must actually attack".

## The chemistry, and why it moves the measured atom

An E2~Ub thioester is ubiquitin's C-terminal Gly76 carboxyl joined to the
catalytic cysteine through sulfur:

```
Ub Gly76 C(=O)–S–Cys91          the loaded enzyme
```

Aminolysis transfers ubiquitin to the substrate amine:

```
substrate–NH2  +  Ub–CO–S–Cys   →   substrate–NH–CO–Ub  +  Cys–SH
```

So with Ub loaded:

* the **electrophile** is **Ub Gly76 C**, *not* Cys91 SG
* the **nucleophile** is unchanged — the LisoK leucyl α-amine, ligand atom `N01`
* **the distance to measure becomes `N01` → `U:76:C`**

That last point is the whole reason this round can say something the others
couldn't. Measuring to Cys91 SG in a charged complex would be measuring to an
atom that is no longer the target.

> **On `NZ`:** your friend's examples bond Ub Gly76 C to a lysine **NZ** — the
> ε-amine of a lysine side chain. That is the *product* isopeptide state of
> conventional ubiquitylation. We attach to **Cys91 SG** instead, which is the
> *pre-transfer* thioester state. Same mechanism, different point along it. The
> pre-transfer state is the one that can discriminate, because in it the
> substrate amine is still free and has to reach.

## The four jobs

Each is 323 protein residues (SUMO2 construct 96 + UBE2W 151 + ubiquitin 76)
plus the LisoK ligand, with 5 seeds.

| file | site | state | bonds |
|---|---|---|---|
| `sumo2_k11lisok_ube2w_ub_charged.json` | K11 — **works** | thioester on Cys91 | A/K12 NZ–L C01; U/G76 C–B/C91 SG |
| `sumo2_k21lisok_ube2w_ub_charged.json` | K21 — **fails** | thioester on Cys91 | A/K22 NZ–L C01; U/G76 C–B/C91 SG |
| `sumo2_k11lisok_ube2w_ub_product.json` | K11 | product isopeptide | A/K12 NZ–L C01; U/G76 C–L N01 |
| `sumo2_k21lisok_ube2w_ub_product.json` | K21 | product isopeptide | A/K22 NZ–L C01; U/G76 C–L N01 |

**Why both states.** The charged jobs are the experiment. The product jobs are a
positive control: if AF3 cannot build a plausible product complex even when told
the isopeptide exists, then a null in the charged state is uninterpretable rather
than informative. This is the lesson from the sortylation round, where a negative
was unreadable because its positive control also failed.

**Why K11 and K21.** They are on the same protein, in the same construct, and
share one MSA — K11 is modified by UbyW, K21 is not. It is the cleanest
within-protein pair in the dataset, and in the uncharged rounds the two were
indistinguishable (3.9 Å vs 3.7 Å closest approach). If loading ubiquitin
separates them, that is the first positive structural signal in the project.

## Construct numbering — please don't "fix" this

The SUMO2 chain carries an **N-terminal proline** and is therefore 96 aa, not 95.
This is from the paper, not a mistake: *"As UBE2W is known to ubiquitylate
SUMO2's native N-terminus we introduced an N-terminal proline residue to prohibit
its modification."* (Results, p. 3.)

That Pro shifts every residue by +1, so **native K11 is construct residue 12**
and native K21 is residue 22. The `bondedAtomPairs` entries use construct
numbering, which is what AF3 sees.

## Running them

Standard AF3, no special flags. The `userCCD` block is embedded in each JSON, so
there is no separate ligand file to pass:

```bash
python run_alphafold.py --json_path=sumo2_k11lisok_ube2w_ub_charged.json \
                        --output_dir=out/ --model_dir=<weights>
```

MSAs: all four jobs share the same three protein sequences, so if you have a way
to compute the data pipeline once and reuse it, that saves three quarters of the
search time. In our hands jackhmmer against the full AF3 database set took 4–6 h
per job because it is I/O-bound on the database reads, while inference itself is
~10 min. If reuse is awkward, just run all four — it is 4 searches, not 20.

## QC

```bash
python qc_charged.py <output_dir> qc_charged.csv
```

Reports per model:

| column | meaning |
|---|---|
| `n01_ub_c` | amine → thioester carbonyl carbon — **the measurement** |
| `n01_cys_sg` | amine → Cys91 SG, for comparability with rounds 2–5 |
| `thioester_len` | Gly76 C → Cys91 SG: did AF3 honour the bond at all |
| `isopeptide_ok` | the XisoK acyl bond onto the substrate lysine formed |
| `nearest_other_C` | specificity — is the amine near *this* carbonyl or just near something |

The script distinguishes "no ubiquitin chain present" (wrong input file) from
"ubiquitin present but the bond did not form" (AF3 declined the constraint),
because those need different fixes.

Tested against a structure with a known ground-truth thioester distance: it
recovers 1.80 Å exactly and reproduces the previously measured 4.86 Å Cys91 SG
distance from the uncharged co-fold.

## The honest caveat

**AF3 has no thioester chemistry.** A `bondedAtomPairs` entry between Gly76 C and
Cys91 SG tells AF3 those atoms are bonded, and it will place them at bonding
distance — but the geometry around that bond is not parameterised as a thioester,
so the carbonyl will not necessarily be planar or correctly oriented. This is a
**geometric proxy** for a loaded E2, not a chemically faithful model of one.

It is still a much better proxy than a bare cysteine, which is what rounds 2–5
used and which we now know is too permissive. But a positive result here would
want confirming with something that does have the chemistry — the MD pipeline in
`ubyw_reactivity_analysis/md/` already has a validated isopeptide residue (`LYQ`)
and could be extended to a thioester.

## What would count as an answer

* **K11 reaches the thioester carbon and K21 does not** → the first structural
  discriminator in this project, and worth extending to more site pairs.
* **Both reach** → the loaded active site is *also* permissive, and reachability
  is confirmed as the wrong observable. Round 2's conclusion stands.
* **Neither reaches, but the product controls build cleanly** → AF3 cannot model
  the pre-transfer geometry, and the approach is exhausted.
* **Neither reaches and the product controls also fail** → uninterpretable; the
  ubiquitin placement itself is the problem, not the substrate geometry.

## Provenance

Sequences are fetched from UniProt at build time and asserted, not typed from
memory — an earlier draft of `build_charged_jobs.py` had a hand-typed UBE2W that
was 155 aa instead of 151, and the assertion caught it before any JSON was
written. Rebuild with:

```bash
python build_charged_jobs.py jobs --ccd LisoK_userCCD.cif
```

`LisoK_userCCD.cif` is the validated ligand from the main repo: L-leucine
isopeptide-linked at `C01`, free α-amine at `N01`, verified L-configuration by
structural comparison against the authoritative L-leucine SMILES rather than by
CIP letter.

**Note on stereochemistry:** we established in round 3 that AF3 **discards
ligand stereochemistry** — L- and D-XisoK inputs produced bit-identical output
coordinates. That does not affect these jobs (LisoK only), but it means any
future stereochemical comparison cannot be done this way.
