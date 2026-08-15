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
  harvested per variant **as each finishes**; a variant with ≥100 models already present
  is skipped on re-run, so a kill costs at most the one variant in flight
* only logs and a `MANIFEST.txt` transfer back — 1200 CIFs would be far too much, and an
  oversized transfer has previously made a job that *succeeded* report failure

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
