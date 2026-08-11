#!/bin/bash
# Test MD run: SUMO2(K11LisoK) + uncharged UBE2W, Amber 24 / pmemd.cuda on Euler.
#
# Protocol follows Dominykas Spelveris's setup, with Amber used end-to-end
# instead of converting to GROMACS -- Euler's only GROMACS module is 2021.4
# (CPU/MPI build) whereas Amber 24 ships pmemd.cuda, so staying in Amber both
# avoids a ParmEd round-trip that could break the custom LYQ linkage and gets
# GPU throughput.
#
#   ff19SB protein, OPC water, 1.5 nm padding, 0.15 M NaCl
#   staged minimisation -> 0.5 ns restrained NVT -> stepwise NPT with restraints
#   released -> production at 310 K / 1 bar, 2 fs, LINCS-equivalent (SHAKE) on
#   H bonds, PME, 1.0 nm cutoffs, Langevin thermostat, Monte Carlo barostat
#
# DELIBERATE DEVIATIONS, both because this is a feel-for-it test not the real run:
#   * 20 ns production per replicate, 3 replicates (60 ns total) rather than
#     3 x 100 ns. Enough to see whether the pose holds and what the C-terminal
#     tail does; not enough for converged MM/GBSA.
#   * Langevin + Monte Carlo barostat rather than velocity-rescale +
#     Parrinello-Rahman. These are the well-tested pmemd.cuda equivalents;
#     Parrinello-Rahman is not available in pmemd.
#
# Usage: sbatch this via the compute harness. Needs the LYQ parameters
# (lyq.prep, lyq.frcmod) and sumo2_lisok_ube2w_amber.pdb in the workdir.
set -eo pipefail
module load amber/24

NS_PROD=${NS_PROD:-20}
NREP=${NREP:-3}
STEPS_PROD=$(python3 -c "print(int($NS_PROD*1e6/2))")   # 2 fs steps

test -f lyq.prep && test -f lyq.frcmod || { echo "LYQ params missing" >&2; exit 1; }
test -f sumo2_lisok_ube2w_amber.pdb || { echo "input PDB missing" >&2; exit 1; }

echo "=== tleap: solvate and parameterise ==="
cat > build.leap <<'LEAP'
source leaprc.protein.ff19SB
source leaprc.water.opc
loadamberparams frcmod.ionsjc_tip4pew
loadamberparams lyq.frcmod
loadamberprep lyq.prep

sys = loadpdb sumo2_lisok_ube2w_amber.pdb
check sys

# vacuum topology, for reference and for MM/GBSA later
saveamberparm sys sys_vac.parm7 sys_vac.rst7

# 1.5 nm padding of OPC water, then neutralise, then 0.15 M NaCl
solvatebox sys OPCBOX 15.0
addions sys Na+ 0
addions sys Cl- 0
saveamberparm sys sys.parm7 sys.rst7
charge sys
quit
LEAP
tleap -f build.leap > tleap_build.log 2>&1 || true
test -f sys.parm7 || { echo "TLEAP FAILED:"; grep -iE 'error|fatal|missing|could not' tleap_build.log | head -20; exit 1; }
grep -iE 'unperturbed charge' tleap_build.log | tail -1

# 0.15 M NaCl on top of neutralisation: n_pairs = 0.15 mol/L * V(L) * N_A
NWAT=$(grep -c 'WAT' sys.pdb 2>/dev/null || python3 - <<'PY'
import re
t=open('tleap_build.log').read()
m=re.findall(r'Added (\d+) residues', t)
print(m[-1] if m else 0)
PY
)
python3 - > ions.txt <<PY
# OPC water is ~55.5 mol/L; n_NaCl = 0.15/55.5 * n_water
nw = $NWAT
print(max(1, round(0.15/55.5*nw)))
PY
NPAIR=$(cat ions.txt)
echo "  water residues ~$NWAT -> adding $NPAIR NaCl pairs for 0.15 M"

cat > salt.leap <<LEAP
source leaprc.protein.ff19SB
source leaprc.water.opc
loadamberparams lyq.frcmod
loadamberprep lyq.prep
sys = loadpdb sumo2_lisok_ube2w_amber.pdb
solvatebox sys OPCBOX 15.0
addions sys Na+ 0
addions sys Cl- 0
addionsrand sys Na+ $NPAIR Cl- $NPAIR
saveamberparm sys sys.parm7 sys.rst7
charge sys
quit
LEAP
tleap -f salt.leap > tleap_salt.log 2>&1 || true
test -f sys.parm7 || { echo "salt step FAILED"; tail -20 tleap_salt.log; exit 1; }
echo "  system atoms: $(python3 -c "
import re; t=open('sys.parm7').read()
i=t.index('%FLAG POINTERS'); print(int(t[i:i+400].split()[6]))" 2>/dev/null || echo '?')"

echo "=== staged minimisation ==="
cat > min1.in <<'IN'
minimise with restrained solute
 &cntrl
  imin=1, maxcyc=5000, ncyc=2500, ntb=1, cut=10.0,
  ntr=1, restraint_wt=10.0, restraintmask='!:WAT,Na+,Cl- & !@H=',
 /
IN
cat > min2.in <<'IN'
unrestrained minimise
 &cntrl
  imin=1, maxcyc=10000, ncyc=5000, ntb=1, cut=10.0, ntr=0,
 /
IN
pmemd.cuda -O -i min1.in -p sys.parm7 -c sys.rst7 -o min1.out -r min1.rst7 -ref sys.rst7
pmemd.cuda -O -i min2.in -p sys.parm7 -c min1.rst7 -o min2.out -r min2.rst7
echo "  minimisation done"

echo "=== 0.5 ns restrained NVT heating to 310 K ==="
cat > heat.in <<'IN'
NVT heat 0 -> 310 K, solute restrained
 &cntrl
  imin=0, irest=0, ntx=1, nstlim=250000, dt=0.002,
  ntc=2, ntf=2, ntb=1, ntp=0, cut=10.0,
  ntt=3, gamma_ln=2.0, temp0=310.0, tempi=10.0, ig=-1,
  ntr=1, restraint_wt=5.0, restraintmask='!:WAT,Na+,Cl- & !@H=',
  ntpr=5000, ntwx=25000, ntwr=50000,
  nmropt=1,
 /
 &wt type='TEMP0', istep1=0, istep2=200000, value1=10.0, value2=310.0 /
 &wt type='END' /
IN
pmemd.cuda -O -i heat.in -p sys.parm7 -c min2.rst7 -o heat.out -r heat.rst7 \
           -x heat.nc -ref min2.rst7
echo "  heating done"

echo "=== stepwise NPT equilibration, restraints released ==="
PREV=heat.rst7
for WT in 5.0 2.0 1.0 0.5 0.1; do
  cat > eq_$WT.in <<IN
NPT equilibration, restraint_wt=$WT
 &cntrl
  imin=0, irest=1, ntx=5, nstlim=250000, dt=0.002,
  ntc=2, ntf=2, ntb=2, ntp=1, barostat=2, pres0=1.0, taup=2.0,
  cut=10.0, ntt=3, gamma_ln=2.0, temp0=310.0, ig=-1,
  ntr=1, restraint_wt=$WT, restraintmask='!:WAT,Na+,Cl- & !@H=',
  ntpr=5000, ntwx=25000, ntwr=50000,
 /
IN
  pmemd.cuda -O -i eq_$WT.in -p sys.parm7 -c $PREV -o eq_$WT.out \
             -r eq_$WT.rst7 -x eq_$WT.nc -ref $PREV
  PREV=eq_$WT.rst7
  echo "  released to $WT kcal/mol/A^2"
done
cat > eq_free.in <<'IN'
NPT equilibration, unrestrained
 &cntrl
  imin=0, irest=1, ntx=5, nstlim=250000, dt=0.002,
  ntc=2, ntf=2, ntb=2, ntp=1, barostat=2, pres0=1.0, taup=2.0,
  cut=10.0, ntt=3, gamma_ln=2.0, temp0=310.0, ig=-1, ntr=0,
  ntpr=5000, ntwx=25000, ntwr=50000,
 /
IN
pmemd.cuda -O -i eq_free.in -p sys.parm7 -c $PREV -o eq_free.out \
           -r eq_free.rst7 -x eq_free.nc
echo "  equilibration complete"

echo "=== production: $NREP x $NS_PROD ns ==="
for R in $(seq 1 $NREP); do
  cat > prod_$R.in <<IN
production, replicate $R
 &cntrl
  imin=0, irest=0, ntx=1, nstlim=$STEPS_PROD, dt=0.002,
  ntc=2, ntf=2, ntb=2, ntp=1, barostat=2, pres0=1.0, taup=2.0,
  cut=10.0, ntt=3, gamma_ln=2.0, temp0=310.0, tempi=310.0, ig=-1,
  ntr=0, ntpr=25000, ntwx=25000, ntwr=250000, ioutfm=1,
 /
IN
  # ntx=1 with ig=-1 gives independent initial velocities per replicate
  pmemd.cuda -O -i prod_$R.in -p sys.parm7 -c eq_free.rst7 -o prod_$R.out \
             -r prod_$R.rst7 -x prod_$R.nc
  echo "  replicate $R done: $(grep -c 'ns/day' prod_$R.out || true)"
  grep -A2 'Average timings' prod_$R.out | tail -3 || true
done

echo "=== bundle ==="
tar czf md_test.tar.gz sys.parm7 sys_vac.parm7 prod_*.nc prod_*.out eq_free.rst7 \
    *.log min2.out heat.out eq_free.out
ls -lh md_test.tar.gz
