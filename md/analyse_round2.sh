#!/bin/bash
# Analyse the round-2 MD: does the AF3 near-attack pose survive a force field,
# and do the reactive and unreactive sites differ?
#
# Run on Euler where the trajectories live:
#   bash analyse_round2.sh
#
# Observables, chosen as the ones AF3 cannot answer:
#
#   1. Nucleophile-to-Cys91-SG distance per frame. The direct question: does the
#      attack geometry hold, and for what fraction of the time?
#   2. Interface contact count. Distinguishes "the substrate came off" from "the
#      interface rearranged while staying bound" -- backbone RMSD alone conflates
#      these, which changed a panel's meaning in the previous round.
#   3. Per-residue CA RMSF, to see whether UBE2W's flexible C-terminus is what
#      moves, as Paul expected.
#   4. Distance to Cys135 and Cys151 as well. The AF3 round showed Cys151 nearly
#      as close as the catalytic Cys91 in the closest model, so "near a cysteine"
#      is not the same as "near THE cysteine".
#
# ATOM SELECTION. The topology is SUMO2 (96 residues) then UBE2W, so UBE2W
# residue n is topology residue 96+n:
#   catalytic Cys91  -> :187@SG
#   Cys135          -> :231@SG
#   Cys151          -> :247@SG
# For the modified systems the nucleophile is the LYQ free alpha-amine, atom NI.
# For the control there is no LYQ: the nucleophile is the standard lysine NZ.
set -o pipefail
W=/cluster/scratch/schnpaul/ubyw/md
OUT=$HOME/ubyw_md_results/analysis
mkdir -p "$OUT"
cd "$W" || exit 1

module load stack/2024-06 amber/24 2>/dev/null || module load amber/24 2>/dev/null
command -v cpptraj >/dev/null || { echo "no cpptraj"; exit 1; }

# Sanity-check the residue offset against the topology rather than trusting it.
python3 - <<'PY' | tee "$OUT/topology_check.txt"
import glob
import re
for p in sorted(glob.glob("*.parm7")):
    txt = open(p, errors="ignore").read()
    m = re.search(r"%FLAG RESIDUE_LABEL.*?%FORMAT\(\S+\)\s*(.*?)%FLAG", txt, re.S)
    if not m:
        print(f"{p}: no RESIDUE_LABEL flag")
        continue
    labels = m.group(1).split()
    lyq = [i + 1 for i, x in enumerate(labels) if x == "LYQ"]
    cys = [i + 1 for i, x in enumerate(labels) if x in ("CYS", "CYX")]
    print(f"{p}: {len(labels)} residues; LYQ at {lyq or '-'}; "
          f"CYS/CYX at {[c for c in cys if c > 96][:6]}")
PY

analyse() {
  local tag="$1" nuc="$2"        # nuc: the nucleophile atom mask
  local reps=()
  for r in 1 2 3; do
    [ -s "${tag}_prod${r}.nc" ] && reps+=("$r")
  done
  if [ ${#reps[@]} -eq 0 ]; then
    echo "  $tag: no trajectories yet" | tee -a "$OUT/analysis.log"
    return 0
  fi
  echo "=== $tag (replicates: ${reps[*]}) ===" | tee -a "$OUT/analysis.log"

  for r in "${reps[@]}"; do
    cat > "cpptraj_${tag}_${r}.in" <<CEOF
parm ${tag}.parm7
trajin ${tag}_prod${r}.nc
autoimage
rms first :1-247@CA out $OUT/rmsd_${tag}_${r}.dat
distance nuc_cat  $nuc :187@SG out $OUT/dist_cat_${tag}_${r}.dat
distance nuc_c135 $nuc :231@SG out $OUT/dist_c135_${tag}_${r}.dat
distance nuc_c151 $nuc :247@SG out $OUT/dist_c151_${tag}_${r}.dat
nativecontacts :1-96 :97-247 distance 4.0 out $OUT/contacts_${tag}_${r}.dat
atomicfluct out $OUT/rmsf_${tag}_${r}.dat :1-247@CA byres
go
quit
CEOF
    cpptraj -i "cpptraj_${tag}_${r}.in" > "$OUT/cpptraj_${tag}_${r}.log" 2>&1
    if [ -s "$OUT/dist_cat_${tag}_${r}.dat" ]; then
      echo "  rep$r ok: $(wc -l < "$OUT/dist_cat_${tag}_${r}.dat") frames" | tee -a "$OUT/analysis.log"
    else
      echo "  rep$r FAILED" | tee -a "$OUT/analysis.log"
      grep -iE 'error|not found' "$OUT/cpptraj_${tag}_${r}.log" | head -3 | tee -a "$OUT/analysis.log"
    fi
  done
}

# The modified systems carry LYQ, whose free alpha-amine (NI) is the nucleophile.
analyse k11_xisok      ":12@NI"
analyse k21_xisok      ":22@NI"
# The control has a standard lysine: NZ is the nucleophile.
analyse k11_lyscontrol ":12@NZ"

# Summarise. FRAME INTERVAL: ntwx=25000 x dt=0.002 ps = 50 ps/frame, verified
# against the mdin files rather than assumed -- mistaking this interval produced a
# 5x error in a previous round.
python3 - <<'PY' | tee "$OUT/summary.txt"
import glob
import os
import re
import statistics as st

ATTACK = 4.0


def frame_ns(tag):
    for p in glob.glob(f"prod_{tag}_*.in") + glob.glob(f"prod_{tag}*.in"):
        t = open(p).read()
        mx = re.search(r"ntwx\s*=\s*(\d+)", t)
        dt = re.search(r"dt\s*=\s*([\d.]+)", t)
        if mx and dt:
            return int(mx.group(1)) * float(dt.group(1)) / 1000.0
    return None


def col2(path):
    v = []
    for ln in open(path):
        if ln.startswith("#"):
            continue
        f = ln.split()
        if len(f) >= 2:
            try:
                v.append(float(f[1]))
            except ValueError:
                pass
    return v


A = os.path.expanduser("~/ubyw_md_results/analysis")
for tag in ("k11_xisok", "k21_xisok", "k11_lyscontrol"):
    fns = frame_ns(tag)
    if fns is None:
        print(f"{tag}: NO mdin found -- refusing to guess the frame interval")
        continue
    per_rep, fracs = [], []
    print(f"\n=== {tag}  ({fns} ns/frame, verified from mdin) ===")
    for r in (1, 2, 3):
        p = f"{A}/dist_cat_{tag}_{r}.dat"
        if not os.path.exists(p):
            continue
        d = col2(p)
        if not d:
            continue
        q = d[int(len(d) * 0.75):]          # final quarter = settled behaviour
        per_rep.append(st.mean(q))
        fr = sum(1 for x in d if x <= ATTACK) / len(d)
        fracs.append(fr)
        print(f"  rep{r}: {len(d) * fns:5.1f} ns  start {d[0]:5.2f}  "
              f"final-quarter {st.mean(q):5.2f} A  <={ATTACK} A in {fr:6.1%} of frames")
    if len(per_rep) >= 2:
        print(f"  across replicates: {st.mean(per_rep):.2f} +- "
              f"{st.stdev(per_rep):.2f} A   attack fraction "
              f"{st.mean(fracs):.1%} +- {st.stdev(fracs):.1%}")
    # is the catalytic cysteine actually the nearest one?
    for lab, key in (("Cys135", "c135"), ("Cys151", "c151")):
        v = []
        for r in (1, 2, 3):
            p = f"{A}/dist_{key}_{tag}_{r}.dat"
            if os.path.exists(p):
                v += col2(p)
        if v:
            print(f"  {lab}: median {st.median(v):.2f} A")
PY

tar czf md_round2_analysis.tar.gz -C "$OUT" . 2>/dev/null
ls -lh md_round2_analysis.tar.gz
