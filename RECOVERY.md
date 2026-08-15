# Recovery: resuming the UBA1 jobs after a VPN drop or a new session

Written 2026-08-15 while two jobs were live on Euler. **Nothing below depends on the
session that submitted them.** Slurm owns the jobs; a dropped VPN, a closed laptop or a
dead kernel cannot affect them.

## Live jobs

**DONE** — `f2ea3a80` (stub validation) and `4e124915` (MSA search) both succeeded.
Findings in [STUB_RESULTS.md](STUB_RESULTS.md); MSAs complete in `~/ubyw_uba1_msa/`.

**RUNNING** — production inference, 12 jobs × 25 seeds, split into 4 parallel batches
(sequential replicates in one job has twice exceeded walltime in this project):

| job id | batch | jobs |
|---|---|---|
| `afeb8aec-769d-4154-b6a3-383a1eca64fc` | 1/4 | `aden_ub75_no_ube2w`, `cys_ub75_no_ube2w`, `thio_ub75_no_ube2w` |
| `4ad60257-502f-49a2-a601-6c0e8718f3b1` | 2/4 | `aden_ub75_with_ube2w`, `cys_ub75_with_ube2w`, `thio_ub75_with_ube2w` |
| `25e500e4-b6e2-4d6e-8bea-d50e99b25269` | 3/4 | `aden_ub76_no_ube2w`, `cys_ub76_no_ube2w`, `thio_ub76_no_ube2w` |
| `309c38e2-4436-4545-a974-142ab1764d27` | 4/4 | `aden_ub76_with_ube2w`, `cys_ub76_with_ube2w`, `thio_ub76_with_ube2w` |

All write their real results into **`$HOME`** on Euler, never scratch:

* MSAs → `~/ubyw_uba1_msa/msa_donor_ub{76,75,74}_data.json` (**complete**)
* production models → `~/ubyw_uba1_models/<job>/<job>__seed-N_sample-M__model.cif`,
  harvested per variant **as each finishes**, so a kill costs at most the one variant in
  flight

> **Resume threshold — read before re-running.** AF3 writes **5 samples per seed**, so
> 25 seeds = **125 models** per variant, not 100. The submitted script's skip test is
> `-ge 100`, which would treat a variant harvested at 100 models (20 seeds) as complete
> and silently leave it 20% short. This does **not** affect the current run — nothing was
> pre-existing, so nothing was skipped — but **any re-run must use `-ge 125`**, or
> better, delete the incomplete variant's directory and let it rebuild. Check with:
> ```
> for d in ~/ubyw_uba1_models/uba1_*; do echo "$(basename $d) $(ls $d/*.cif|wc -l)/125"; done
> ```
> A variant showing anything other than 125 is incomplete regardless of what the job's
> exit status said.
* only logs and a `MANIFEST.txt` transfer back — 1200 CIFs would be far too much, and an
  oversized transfer has previously made a job that *succeeded* report failure

## State as of the last check (2026-08-15, ~15:50)

| | |
|---|---|
| MSAs | **complete**, 110-111 MB each, in `~/ubyw_uba1_msa/` |
| batches 1-3 | **running**, ~1 h elapsed, **~21.8 h walltime remaining** each |
| batch 4 (`309c38e2`) | **queued**, blocked on `QOSMaxGRESPerUser` — normal, it starts as a slot frees |
| models harvested | 0 so far (first variant of each batch still in inference) |
| home quota | 20.2 GB of 45 GB soft — 1200 CIFs will not threaten it |

Walltime headroom is the number that matters overnight: ~21.8 h remaining against a job
that needs roughly 3 × 25 × 9 min ≈ 11 h. Comfortable, but if a batch does hit the wall,
the per-variant harvest means the finished variants are already safe in `$HOME` and only
the one in flight is lost.

## Telling a VPN drop apart from a job failure

Seen for real on 2026-08-15 at ~22 h of job wall time. The signature is unambiguous:

```
ssh: connect to host euler.ethz.ch port 22: Operation timed out
```

That is **the tunnel, not the jobs.** Slurm keeps running them; you have simply lost the
ability to look. Reconnect the VPN and re-run the same command.

**Two diagnostic traps, both of which cost me a turn:**

1. **A failed remote command can return empty stdout with an apparently clean status.**
   The first call after the drop printed *nothing at all* — no error, no rows — which
   reads like "the queue is empty and the model directories are gone", i.e. catastrophe,
   when the truth was "the command never ran". Always check `stderr` and the exit code,
   not just stdout. An empty result from a command that should always print something
   (`hostname`, `ls -d ~`) means the connection failed, not that the data vanished.
2. **The ledger's `state` is the last *polled* state, not live truth.** It said `running`
   for all four batches — correct as of the last successful poll, but the daemon's poller
   cannot reach the host either while the tunnel is down. Treat ledger state during an
   outage as "last known", and re-confirm against `squeue` once reconnected.

## Why a VPN drop is safe

1. **Slurm keeps running the job.** The submission is complete; nothing streams from
   this machine. Losing the tunnel loses only *your view*.
2. **Results land in `$HOME`**, which is not purged. (`/cluster/scratch` is purged after
   15 days — that is why nothing important is left there.)
3. **The MSA job is incremental and resumable.** Each donor's `_data.json` is copied to
   `~/ubyw_uba1_msa/` the moment that donor finishes, and on re-run any donor whose file
   is already present and non-empty is skipped. So even a hard kill costs at most the
   one search that was mid-flight.
4. **The big files are deliberately not transferred.** The augmented JSONs are ~34 MB
   each; a 34 MB transfer previously made a job that had *succeeded* report failure. Only
   small text summaries come back, so the harvest cannot fail on size.

## To resume — the short version

```python
# repl tool
c = host.compute.create("euler")
for j in host.compute.ledger().jobs:
    print(j)                      # states of everything this conversation owns
```

Then, per job:

```python
res = c.attach_job("f2ea3a80-9f4e-4a55-8e7d-a6182639ab3b").result()
print(res.state, res.exit_code, res.stdout_tail)
```

`.result()` **raises `JobPending` until the job is terminal** — that means "still
running", not an error. Do not poll it in a loop; park on `wait_for_notification`
instead, which is how the completion actually arrives.

If this conversation's ledger is empty (a genuinely new session), the jobs are still
findable on the host directly:

```python
c.call_command("squeue -u $USER; ls -lh ~/ubyw_uba1_msa/", intent="check UBA1 job state")
```

**First `call_command` after reconnecting the VPN may time out at 60 s while the tunnel
settles.** That is expected — retry after ~2 min rather than treating it as a
misconfiguration.

## The three questions the stub must answer

Read these from `harvest/*.log` and the returned `.cif`:

1. **`grep -i "Reducing number of bonds"` → must be absent.** If present, AF3 discarded a
   bond and the model is meaningless. Exit code is 0 either way.
2. **Thioester geometry at the acyl carbon:** C–S length (~1.78-1.81 Å) and the
   bond-angle sum. sp2 is 360°, sp3 is 328.5°. The earlier charged run gave 338.6°
   *with* a valence bug; with `UBGG`/`UBG1` (no OXT) it should sit much closer to planar.
   If it still collapses toward 328°, that is AF3's prior and not our chemistry — which
   is itself the argument for doing this in MD.
3. **Chain renumbering.** AF3 absorbs a leading single-residue ligand into the protein
   chain (last time chain `U` came back as 75 residues, not 74, with the second glycine
   as chain `UA`). **Resolve chains from the `covale` records, never from input naming.**

## Then: production inference (not yet submitted)

Gated on the stub being clean *and* the MSAs existing. Sequence:

1. Confirm `~/ubyw_uba1_msa/` holds all three `_data.json` files, non-empty.
2. Graft the MSAs onto the 12 production jobs **on the cluster** (the donors are there,
   and they are too big to move):
   `python graft_msa.py jobs_uba1/<job>.json ~/ubyw_uba1_msa/msa_donor_ub<N>_data.json out.json`
   — grafting matches **by sequence, not chain id**, and refuses if any recipient chain
   finds no exact match. Pick the donor by the job's ubiquitin length: `thio` jobs use
   Ub 1-74, `aden`/`cys` jobs use 1-76 or 1-75 as named.
3. Submit inference with `--norun_data_pipeline` and `--nv` (required for inference;
   must be *omitted* for the CPU MSA stage).
4. Size GPU memory from what the stub actually used — these are up to 1285 residues,
   4× the earlier charged jobs, so the old `gpumem:24g` guidance does not carry over.
   The stub requested 40 GB.
5. Run replicates as **parallel jobs**, not sequentially in one job: walltime arithmetic
   has twice exceeded the limit in this project.
6. Harvest with `<job>__<seed-sample>__model.cif`. Every AF3 model is named `model.cif`
   inside its own `seed-N_sample-M` directory, so a flat `cp` collapses 126 files onto
   one name — and a failing `cp` exit status has previously stopped the remaining
   variants. Return 0 unconditionally from the harvest and detect
   already-complete variants so they are harvested rather than recomputed.

## When the models land — the pipeline is already built and tested

Both scripts were written and run against the stub output, so they are known to work
before real data arrives. On the cluster (stdlib only, no numpy/gemmi needed):

```bash
python qc_uba1.py ~/ubyw_uba1_models qc_uba1.csv
```

Then locally, where pandas/matplotlib are available:

```bash
python analyse_uba1.py qc_uba1.csv results_uba1/
```

`qc_uba1.py` records `acyl_source` on every row and **raises** on an unrecognised layout
rather than falling back — the failure that made `analyse_ub75.py` report a thioester at
6.80 Å where the truth was 1.68 Å. `analyse_uba1.py` prints the AUC direction in the same
breath as the number, prints the condition count and measured ICC beside the model count,
and reports medians with bootstrap CIs and fraction-within-cutoff rather than minima.

Sanity checks to run on the QC output **before** interpreting anything:

* `acyl_source` distribution matches the job types — 4 variants should show
  `UBGG:C2`/`UBG1:C1` (the covalent thioester jobs), the other 8 `protein_Cterm:*`;
* `cys632_is_cys` true in every row;
* declared thioester present in ~100% of covalent models, planarity near 357-360°;
* 125 models per variant.

## Local state, if the workspace is swept

The workspace is scratch and gets swept; the artifacts and the git branch are durable.

* Branch **`uba1-sites`** of `ubyw_charged_e2_af3` holds everything:
  `make_uba1_ccd.py`, `build_uba1_jobs.py`, `make_stubs.py`, `ccd_uba1/`,
  `jobs_uba1/` (12), `stubs_uba1/` (3), `msa_donors/` (3), `seqs_uba1.json`,
  `UBA1_JOBS.md`, this file.
* Regenerate everything deterministically:
  ```bash
  PYTHONPATH=. python make_uba1_ccd.py ccd_uba1
  cp ccd_v2/UBGG_userCCD.cif ccd_uba1/
  PYTHONPATH=. python build_uba1_jobs.py jobs_uba1 --seeds 25
  PYTHONPATH=. python make_stubs.py jobs_uba1/uba1_thio_ub76_with_ube2w.json \
      jobs_uba1/uba1_thio_ub75_with_ube2w.json \
      jobs_uba1/uba1_aden_ub76_with_ube2w.json stubs_uba1
  ```
  Needs `rdkit` (conda-forge). The scripts import `ccd_valence` from the repo root, and
  the kernel here strips the script directory from `sys.path`, so **`PYTHONPATH=.` is
  required**.

## Open item for the user

The adenylate parameterisation for MD is **not** blocking these AF3 jobs and was
deliberately deferred. The literature search results and the papers to obtain are in
[ADENYLATE_LITERATURE.md](ADENYLATE_LITERATURE.md).
