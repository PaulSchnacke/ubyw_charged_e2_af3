# Where the UBA1 models live, and how to archive them

Measured 2026-08-16, not estimated.

## They are NOT on scratch

`~/ubyw_uba1_models/` is on `/cluster/home/schnpaul` — a 16 TB NFS volume, backed up, no
purge timer. This was deliberate: the submit scripts harvest each variant into `$HOME` as
it finishes, precisely so that a purge, a VPN drop, or a killed job cannot lose work.

Scratch holds only the job *working* directories (`/cluster/scratch/$USER/.claude-science/
jobs`, 8.4 GB across 453 dirs). Verified explicitly: of the UBA1 `model.cif` files there,
**every production model also exists in `$HOME`**. The only 15 that do not are the
single-sequence **stub validation** models from job `f2ea3a80`, which were already
transferred to the analysis workspace and whose findings are in `STUB_RESULTS.md`. Losing
them to the purge costs nothing.

**One thing worth acting on:** scratch also holds ~800 models from a *different* project
(`1t2w_*`, `7s4o_*`, `8t8g_*`, `7s51_*` — cognate / scramble / wrongclass / apo series).
Those are **not** mirrored in `$HOME` as far as this check went, and scratch purges after
15 days. If they matter, move them.

## Sizes, measured

| | |
|---|---|
| models, all 12 variants (1500 CIF) | **1.22 GB** |
| MSAs (3 donor sets) | **0.32 GB** |
| per variant (125 models) | ~94–116 MB |
| per model | ~760–930 KB, 8,900–10,100 atoms, 12,500 lines |
| home quota | 21.6 GB of 45 GB soft — currently no pressure |

Each CIF carries the full atom table plus a per-atom pLDDT column, which is where most of
the bulk sits.

## Compression, measured on one full variant (125 models, 94.0 MB raw)

| method | size | ratio |
|---|---|---|
| `tar.gz -6` | 24.3 MB | 3.9× |
| `tar.zst -19` | 12.5 MB | 7.5× |
| **`tar.xz -6`** | **11.0 MB** | **8.5×** |

CIF is plain text with highly repetitive coordinate formatting, so it compresses well.
Extrapolating to the whole set: **1.22 GB → ~145 MB as `tar.xz`**, or ~165 MB with `zstd`.

## Recommendation

**Archive as one `tar.xz` per variant (~11 MB each, ~145 MB total).** Per-variant rather
than one monolithic archive so you can retrieve a single condition without unpacking 1.5 GB.

```bash
cd ~/ubyw_uba1_models
for d in uba1_*; do
  tar -cf - "$d" | xz -6 -T4 > "$HOME/ubyw_uba1_archive/${d}.tar.xz"
done
sha256sum ~/ubyw_uba1_archive/*.tar.xz > ~/ubyw_uba1_archive/CHECKSUMS.sha256
```

Then verify before deleting anything:

```bash
cd ~/ubyw_uba1_archive && sha256sum -c CHECKSUMS.sha256
tar -tJf uba1_cys_ub76_no_ube2w.tar.xz | wc -l    # expect 127 (125 models + top-ranked + log)
```

At 145 MB the whole set fits comfortably anywhere — institutional storage, a Zenodo deposit
(50 GB limit), even a repo release asset (2 GB limit). **Do not commit it to the git repo
itself**: git stores every version in full, so a re-run would double the clone size forever.

**Keep uncompressed indefinitely** the three derived CSVs — `results/qc_uba1.csv`,
`results/uba1_aden_reach.csv`, `results/uba1_differential_contacts.csv`. They are ~300 KB
together and carry every measurement the report depends on, so the archive only needs
unpacking if someone wants to re-measure something new.

**Also archive the MSAs (0.32 GB).** They are the expensive part to regenerate — a fresh
UBA1 search is hours of CPU — and they compress poorly (already near-random alignment
text), so budget the full 0.32 GB. Given AF3 is deterministic per seed, MSAs + job JSONs +
seeds is in principle enough to regenerate every model, which makes them arguably more
valuable per byte than the models themselves.

## Representative models

Five structures, one per chemistry in the run, top-ranked model (`seed-1_sample-0`) each:

| file | what it shows |
|---|---|
| `uba1_aden_ub76_no_ube2w.cif` | adenylation site with ATP·Mg, full-length tail |
| `uba1_aden_ub75_no_ube2w.cif` | same, short tail — the direct comparison |
| `uba1_thio_ub76_with_ube2w.cif` | covalent thioester at Cys632, full tail, UBE2W present |
| `uba1_thio_ub75_with_ube2w.cif` | same, short tail |
| `uba1_cys_ub76_with_ube2w.cif` | free Cys632, non-covalent — no declared bond |

Reminder when inspecting these: geometry **at or adjacent to a declared bond is not
physical** (AF3 treats `bondedAtomPairs` as connectivity — the thioester C–S comes out at
~1.5 Å where 1.78–1.81 is correct, and the bonded cysteine's own Cβ–Sγ is distorted too).
The placement is meaningful; the junction bond lengths are not.
