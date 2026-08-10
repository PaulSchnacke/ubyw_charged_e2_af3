# Overnight sweep — what protects it from a dropped VPN

SLURM jobs are unaffected by the tunnel going down; only *polling* is. The real
risks are different, and each has a specific guard.

## Jobs in flight

| group | jobs | id |
|---|---|---|
| `k11species` | 5 species at K11: xisok, lyscontrol, charged, tetrahedral, product | `44db2930` |
| `sites` | K5, K7, K21, K33, K35, K42, K45 (xisok) | `423495fb` |

20 seeds x 5 samples = 100 models per job, 1200 total. **No database search:** all
three MSA sets were already computed in earlier rounds and are grafted on, so the
whole night is inference (~10 min per job on this hardware).

## The five safeguards

1. **Results land in `$HOME`, not scratch.** Everything writes to
   `$HOME/ubyw_sweep_results/<group>/` — progress log, per-job AF3 logs, and the
   QC CSV. Immune to the 15-day scratch purge, and readable next time regardless
   of what happened to the connection.

2. **QC runs after *every* job, not at the end.** `qc_partial.csv` is rewritten
   each time a job finishes, so a night that dies at job 4 of 7 still leaves
   analysed results for 1-4 rather than a directory of unparsed CIFs.

3. **Analysis happens on the cluster.** 1200 CIFs is far past the harvest
   threshold — a 34 MB file already caused a *successful* job to report failure
   earlier in this project. The transferred artefact is a few hundred kB of CSV.

4. **No `set -e` in the loop.** Each job is wrapped in `run_one` which always
   returns 0, so one failure costs one job rather than the night. Failures are
   recorded in the progress log with their cause.

5. **The silent failure mode is checked explicitly.** After every run the log is
   grepped for `Reducing number of bonds` — AF3 3.0.1 drops polymer-polymer bonds
   with exit 0 and a plausible-looking model. A dropped bond is written to the
   progress log as a WARNING rather than being discovered weeks later.

Plus: MSA grafts are matched **by sequence, not chain id**, and refuse to write a
job with partial alignments. Verified before launch — all four donor/recipient
pairs matched exactly (`MISSING []`).

## Two failures caught before you left

**Exit 9 on the first attempt: `#SBATCH` directives must be the FIRST lines.**
I put a `cp` above the block, which demoted every directive to an ordinary
comment. The job then ran on the **login node** with no GPU, and `nvidia-smi`
failing killed it. The working validation job had its directives at the top.
Any command above the `#SBATCH` block silently costs you the allocation.

**`nvidia-smi` is now guarded.** Its failure was what turned a misconfiguration
into a dead job. The GPU query redirects errors and defaults to *passing*
`--flash_attention_implementation=xla` when the compute capability cannot be
read — harmless on A100, required on anything older.

Related, from the validation run: AF3 3.0.1 **refuses to start** on compute
capability < 8.0 (V100, TITAN RTX) without that flag. The scheduler allocates by
availability — this sweep got a TITAN RTX (7.5), so the flag is load-bearing, not
defensive. Earlier jobs in this project got A100s by luck.

## Reading the results afterwards

```bash
ssh euler
tail -20 ~/ubyw_sweep_results/k11species/progress.log
grep -i 'WARNING\|GRAFT FAILED' ~/ubyw_sweep_results/*/progress.log
column -s, -t ~/ubyw_sweep_results/*/qc_*.csv | less -S
```

The QC resolves chains **from the covale records**, not from the input naming:
AF3 renumbers split chains (it absorbed a GLY ligand into the ubiquitin chain
last time, making it 75 residues rather than 74) and renames ligand chains.
Reading the input spec would measure the wrong atoms.

## What to look at first

For the K11 species, the question is whether the corrected chemistry changes the
picture:

* `charged` — is the free amine still 26 A from the thioester carbon across 100
  models, or does it approach now that the carbon has a correct valence? Approach
  with correct valences would be a *prediction*; before, it was an artefact.
* `tetrahedral` vs `charged` — the same complex with and without the bond
  declared, so the comparison isolates what declaring it does.
* `lyscontrol` — the unmodified lysine, which in the earlier round sat far from
  the active site while the modified one docked. Expect that to hold.

For the sites, the honest framing: K11 works and K21 fails per the curated table;
the other six are **untested**, so they are a blind set. Rounds 1-5 were all
negative and the prior is that reach does not predict reactivity. If the six
untested sites split cleanly from K11 that is interesting; if they do not, it is
the fourth consistent negative and worth accepting as an answer.
