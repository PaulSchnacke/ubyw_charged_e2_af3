#!/bin/bash
# Analyse all three production replicates of the SUMO2(K11LisoK)+UBE2W test MD.
#
# Observables, chosen as the ones AF3 cannot answer:
#   1. LisoK amine (LYQ NI) to UBE2W catalytic Cys91 SG -- is the AF3 attack
#      pose (4.86 A in the input model, 4.6 A median across 300 AF3 samples) a
#      minimum the system stays in?
#   2. Cxa RMSD and substrate-enzyme contacts -- does the complex hold, i.e. is
#      any relaxation of (1) local rather than the substrate coming off?
#   3. Per-residue RMSF, especially UBE2W's C-terminal tail (141-151), where AF3
#      pLDDT falls 91 -> 41.
#
# Residue indices are RESOLVED, not assumed: the topology concatenates SUMO2
# (96 res) then UBE2W (151 res), so UBE2W Cys91 is topology residue 187 -- but
# this script asserts that residue 187 is a CYS in molecule 2 and that residue
# 12 is LYQ before measuring anything.
set -eo pipefail
module load amber/24

test -f sys.parm7 || { echo "sys.parm7 missing" >&2; exit 1; }
REPS=""
for r in 1 2 3; do [ -f prod_$r.nc ] && REPS="$REPS $r"; done
test -n "$REPS" || { echo "no prod_*.nc found" >&2; exit 1; }
echo "replicates present:$REPS"

# --- verify the residue identities the masks rely on -------------------------
cat > check.in <<'CPT'
parm sys.parm7
resinfo :12
resinfo :187
CPT
cpptraj -i check.in > check.log 2>&1
grep -A3 '#Res' check.log | awk '/LYQ|CYS/{print "  resolved: " $0}'
# The identity line sits THREE lines after the `resinfo :N` echo (the echo, a
# "1 residues selected." line, the #Res header, then the data), so -A2 misses it
# and the assertion fails on a correct topology. Parse the data line by its
# residue number instead of by proximity to the command echo.
assert_res () {   # assert_res <resnum> <expected name>
  local got
  got=$(awk -v n="$1" '$1==n && NF>=6 {print $2; exit}' check.log)
  [ "$got" = "$2" ] || { echo "topology residue $1 is '$got', expected $2" >&2; exit 1; }
}
assert_res 12 LYQ
assert_res 187 CYS

# --- per-replicate ------------------------------------------------------------
for r in $REPS; do
  cat > an_$r.in <<CPT
parm sys.parm7
trajin prod_$r.nc
autoimage
distance NI_SG :12@NI :187@SG out dist_rep$r.dat
rms CA_rmsd first @CA out rmsd_rep$r.dat
atomicfluct out rmsf_rep$r.dat @CA byres
nativecontacts :1-96 :97-247 distance 4.0 out contacts_rep$r.dat
# is the catalytic Cys really the nearest one, or a coincidence? measure the
# distance to every cysteine SG in the complex
distance c135 :12@NI :231@SG out cys135_rep$r.dat
distance c151 :12@NI :247@SG out cys151_rep$r.dat
run
CPT
  cpptraj -i an_$r.in > cpptraj_$r.log 2>&1 \
    || { echo "cpptraj failed on rep $r"; tail -20 cpptraj_$r.log; exit 1; }
  N=$(($(wc -l < dist_rep$r.dat) - 1))
  echo "  rep $r: $N frames = $(python3 -c "print(f'{$N*0.05:.1f}')") ns"
done

# --- pooled summary ----------------------------------------------------------
python3 - <<'PY'
import glob, numpy as np
FRAME_NS = 0.05          # ntwx=25000 x dt=0.002 ps, verified against prod_*.in
AF3_INPUT, AF3_MEDIAN, ATTACK = 4.86, 4.6, 4.0

rows, finals = [], []
for f in sorted(glob.glob("dist_rep*.dat")):
    d = np.loadtxt(f, skiprows=1)[:, 1]
    r = f.replace("dist_rep", "").replace(".dat", "")
    q = d[int(len(d) * 0.75):]                 # final quarter
    finals.append(q)
    rows.append((r, len(d), len(d) * FRAME_NS, d[0], d.mean(), q.mean(), q.std(),
                 d.min(), (d <= ATTACK).mean(), (d <= 8).mean()))

print("\nLYQ NI to UBE2W Cys91 SG")
print(f"{'rep':>4} {'ns':>6} {'start':>7} {'mean':>7} {'final-q':>9} {'min':>6} "
      f"{'<=4A':>7} {'<=8A':>7}")
for r, n, ns, s, m, fq, fsd, mn, f4, f8 in rows:
    print(f"{r:>4} {ns:6.1f} {s:7.2f} {m:7.2f} {fq:6.2f}+-{fsd:.2f} {mn:6.2f} "
          f"{f4:7.1%} {f8:7.1%}")
allq = np.concatenate(finals)
per_rep = np.array([q.mean() for q in finals])
print(f"\npooled final-quarter: {allq.mean():.2f} A "
      f"(across-replicate mean {per_rep.mean():.2f} +- {per_rep.std(ddof=1) if len(per_rep)>1 else 0:.2f})")
print(f"AF3 input model {AF3_INPUT:.2f} A; AF3 median across 300 samples {AF3_MEDIAN:.1f} A")
print(f"verdict: {'pose RETAINED' if allq.mean() < 6 else 'pose RELAXED AWAY' if allq.mean() > 10 else 'INTERMEDIATE'}")

# specificity: is Cys91 still the nearest cysteine?
for tag, lab in (("cys135", "Cys135"), ("cys151", "Cys151")):
    fs = sorted(glob.glob(f"{tag}_rep*.dat"))
    if fs:
        v = np.concatenate([np.loadtxt(f, skiprows=1)[:, 1] for f in fs])
        print(f"  NI to {lab}: mean {v.mean():.1f} A  min {v.min():.1f} A")

print("\nRMSF (Calpha, A)")
for f in sorted(glob.glob("rmsf_rep*.dat")):
    a = np.loadtxt(f, skiprows=1)
    sub = a[a[:, 0] <= 96, 1]
    core = a[(a[:, 0] > 96) & (a[:, 0] <= 236), 1]
    tail = a[a[:, 0] > 236, 1]
    print(f"  {f.split('_')[1][:-4]:>5}: SUMO2 {sub.mean():.2f}  "
          f"UBE2W core {core.mean():.2f}  UBE2W tail {tail.mean():.2f}  "
          f"(ratio {tail.mean()/core.mean():.1f}x)")

print("\nCalpha RMSD from start (A)")
for f in sorted(glob.glob("rmsd_rep*.dat")):
    a = np.loadtxt(f, skiprows=1)[:, 1]
    print(f"  {f.split('_')[1][:-4]:>5}: mean {a.mean():.2f}  final {a[-1]:.2f}  max {a.max():.2f}")
PY

tar czf md_analysis.tar.gz dist_rep*.dat rmsd_rep*.dat rmsf_rep*.dat \
    contacts_rep*.dat cys*_rep*.dat cpptraj_*.log prod_*.in
ls -lh md_analysis.tar.gz
