# ubyw_charged_e2_af3

AF3 handover jobs: **UBE2W charged with ubiquitin** (via the ligand trick)
co-folded with a XisoK-bearing SUMO2 substrate.

Companion to [ubyw_reactivity_analysis](https://github.com/PaulSchnacke/ubyw_reactivity_analysis),
which established that a **bare** catalytic cysteine is too permissive a target:
AF3 docks the modification into an empty active site whether or not the site is
reactive (AUC 0.367 across 14 sites), and MD showed that pose is not a
force-field minimum (7.4 ± 0.3 Å, 0/1200 frames ≤ 4 Å).

**Full documentation: [HANDOVER.md](HANDOVER.md).**
Ligand-trick adaptation notes: [docs/LIGAND_TRICK_ADAPTED.md](docs/LIGAND_TRICK_ADAPTED.md).

## Run these four

```bash
python run_alphafold.py --json_path=jobs/sumo2_k11lisok_ube2w_ub_charged.json \
                        --output_dir=out/ --model_dir=<weights>
```

| file | site | state |
|---|---|---|
| `jobs/sumo2_k11lisok_ube2w_ub_charged.json` | K11 (works) | thioester on Cys91 — **the experiment** |
| `jobs/sumo2_k21lisok_ube2w_ub_charged.json` | K21 (fails) | thioester on Cys91 — **the experiment** |
| `jobs/sumo2_k11lisok_ube2w_ub_product.json` | K11 | product isopeptide — positive control |
| `jobs/sumo2_k21lisok_ube2w_ub_product.json` | K21 | product isopeptide — positive control |

Each is self-contained: 323 protein residues, 5 seeds, `userCCD` embedded.

## Then QC

```bash
python qc_charged.py out/ qc_charged.csv
```

**The measurement is `L:1:N01` → `U:76:C`** — the LisoK α-amine to ubiquitin's
Gly76 thioester carbonyl. With ubiquitin loaded, that carbon is the electrophile;
Cys91 SG no longer is.
