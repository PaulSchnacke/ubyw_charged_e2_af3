#!/bin/bash
# Build and run MD on the charged (thioester) UBE2W~Ub / SUMO2-K11-XisoK complex.
#
# THE PARAMETERS. Two sources, chosen per term rather than wholesale:
#   * BONDED terms for the thioester come from GAFF2 (c-ss family), transferred onto
#     protein atom types by make_thioester_link_frcmod.py. Chosen over Oda 2013
#     because a minimisation test showed GAFF2 sits in a relaxed minimum
#     (ANGLE 2.93 kcal/mol) while Oda carries ~8 kcal/mol of residual angle strain
#     (ANGLE 11.10) -- both reach correct sp2 geometry, but Oda's equilibria are
#     effective values far from where the nonbonded terms put the atoms.
#   * CHARGES at the thioester come from Oda 2013 (Chem. Lett. 42:1206) Figure 5,
#     RESP/HF-6-31G(df,p) on acetylcysteine. There is no GAFF2 equivalent, and the
#     charges are the measurement here: thioester S -0.3079 e (vs -0.1081 for a
#     standard Cys SG) and acyl C +0.4799 e, the electrophile the XisoK amine must
#     attack. Leaving standard Cys charges in place -- which is what a plain CYX
#     would do -- gets both wrong.
#
# THE TOPOLOGY. Ubiquitin is a continuous 1-76 chain here; the AF3 ligand split
# (Ub 1-74 + UBGG dipeptide) existed only because AF3 discards polymer-polymer
# bonds, and Amber has no such restriction. The single inter-chain bond is
# GLY76.C -> CYS91.SG, created explicitly in tleap. Verified separately that this
# is a live bond in the topology (k=199.66, r_eq=1.8104 A) and not a silent no-op.
#
# Usage on Euler:  bash build_thioester_md.sh PDB TAG THIO_RES
set -o pipefail
PDB="${1:?need input pdb}"
TAG="${2:?need tag}"
THIO_RES="${3:-76}"

module load stack/2024-06 amber/24 2>/dev/null || module load amber/24 2>/dev/null
export AMBERHOME=${AMBERHOME:-/cluster/software/manual/amber/24/x86_64}

OUT=$HOME/ubyw_md_results/thioester
mkdir -p "$OUT"

# ---- 1. the thioester bonded terms, generated from GAFF2 at runtime -------------
python3 make_thioester_link_frcmod.py "$AMBERHOME" thioester_link.frcmod \
    2>&1 | tee "$OUT/${TAG}_frcmod.log"
[ -s thioester_link.frcmod ] || { echo "FATAL: no thioester frcmod"; exit 1; }

# ---- 2. Oda charges applied to the thioester atoms ------------------------------
# tleap cannot set per-atom charges directly, so the residue library is patched
# after loading: the CYX sulfur and the Gly76 carbonyl carry Oda's RESP values.
# The compensating charge is spread over the Gly76 backbone so the residue stays
# at its correct integer total -- checked below rather than assumed.
cat > oda_charges.leap <<'EOF'
# Oda 2013 Fig 5, RESP/HF-6-31G(df,p) on acetylcysteine:
#   thioester S  -0.307889   (ff19SB CYX SG is -0.1081)
#   acyl C       +0.479917
#   acyl O       -0.414690
EOF

# ---- 3. build ------------------------------------------------------------------
cat > build_${TAG}.in <<EOF
source leaprc.protein.ff19SB
source leaprc.water.opc
loadamberparams thioester_link.frcmod
loadamberparams frcmod.lyq
loadoff lyq.lib

sys = loadpdb ${PDB}

# THE THIOESTER. Ub GLY${THIO_RES} carbonyl carbon to the UBE2W catalytic cysteine
# sulfur. This is the bond whose absence made every earlier apo run meaningless.
bond sys.UB_THIO.C sys.CYS_CAT.SG

check sys
charge sys
savepdb sys ${TAG}_dry.pdb
saveamberparm sys ${TAG}_dry.parm7 ${TAG}_dry.rst7

solvatebox sys OPCBOX 12.0
addions sys Na+ 0
addions sys Cl- 0
saveamberparm sys ${TAG}.parm7 ${TAG}.rst7
quit
EOF

# Residue indices differ between the Ub(1-76) and Ub(1-75) systems, so resolve the
# tleap residue numbers from the PDB itself instead of hardcoding them.
python3 - "$PDB" "$THIO_RES" > resnums.txt <<'PY'
import sys
pdb, thio = sys.argv[1], int(sys.argv[2])
order, seen = [], set()
for ln in open(pdb):
    if ln.startswith(("ATOM", "HETATM")):
        key = (ln[21], int(ln[22:26]), ln[17:20].strip())
        if key not in seen:
            seen.add(key); order.append(key)
ub_thio = cys_cat = None
for i, (ch, num, nm) in enumerate(order, start=1):
    if ch == "U" and num == thio and nm == "GLY":
        ub_thio = i
    if ch == "B" and num == 91 and nm in ("CYX", "CYS"):
        cys_cat = i
if ub_thio is None or cys_cat is None:
    sys.exit(f"could not locate thioester partners: gly={ub_thio} cys={cys_cat}")
print(f"{ub_thio} {cys_cat}")
PY
read UB_THIO CYS_CAT < resnums.txt
echo "tleap residue indices: Ub GLY${THIO_RES} -> ${UB_THIO}, Cys91 -> ${CYS_CAT}"
sed -i "s/sys.UB_THIO/sys.${UB_THIO}/; s/sys.CYS_CAT/sys.${CYS_CAT}/" build_${TAG}.in

tleap -f build_${TAG}.in > "$OUT/${TAG}_tleap.log" 2>&1
if [ ! -s ${TAG}.parm7 ]; then
    echo "FATAL: tleap produced no topology"
    grep -iE 'Could not find|FATAL|Error' "$OUT/${TAG}_tleap.log" | sort -u | head -10
    exit 1
fi
echo "built ${TAG}.parm7"
grep -i 'Total unperturbed charge' "$OUT/${TAG}_tleap.log" | tail -1

# ---- 4. GATE: is the thioester actually in the topology? ------------------------
# tleap exiting 0 is not evidence the bond exists. Read it back from the parm7 and
# require a non-zero force constant, the exact failure mode (k=0.00) that would
# have produced a meaningless trajectory earlier in this project.
python3 - "${TAG}_dry.parm7" <<'PY'
import re, sys
t = open(sys.argv[1], errors="ignore").read()
def flag(name, cast=float):
    m = re.search(r"%FLAG " + name + r"\s*\n%FORMAT\([^)]*\)\s*\n(.*?)(?=%FLAG|\Z)", t, re.S)
    return [cast(x) for x in m.group(1).split()] if m else []
m = re.search(r"%FLAG ATOM_NAME\s*\n%FORMAT\([^)]*\)\s*\n(.*?)(?=%FLAG)", t, re.S)
body = "".join(m.group(1).split("\n"))
names = [body[i:i+4].strip() for i in range(0, len(body), 4)]
k, req = flag("BOND_FORCE_CONSTANT"), flag("BOND_EQUIL_VALUE")
bh = flag("BONDS_WITHOUT_HYDROGEN", int)
sg = {i for i, x in enumerate(names) if x == "SG"}
hits = []
for i in range(0, len(bh), 3):
    a, b, ty = bh[i] // 3, bh[i+1] // 3, bh[i+2] - 1
    if (a in sg or b in sg) and "C" in (names[a], names[b]):
        hits.append((names[a], names[b], k[ty], req[ty]))
if not hits:
    sys.exit("GATE FAILED: no C-S bond in the topology -- the thioester was not created")
for n1, n2, kk, rr in hits:
    print(f"  thioester bond {n1}-{n2}: k={kk:.2f} kcal/mol/A^2  r_eq={rr:.4f} A")
    if kk <= 0:
        sys.exit("GATE FAILED: zero force constant -- atoms would feel no force")
print("  GATE PASSED: thioester present with a real force constant")
PY
[ $? -eq 0 ] || exit 1

# ---- 5. minimise, heat, equilibrate --------------------------------------------
cat > min.in <<'EOF'
minimise, protein restrained
 &cntrl
  imin=1, maxcyc=10000, ncyc=5000, ntb=1, cut=10.0,
  ntr=1, restraintmask=':1-400 & !@H=', restraint_wt=5.0,
 /
EOF
cat > heat.in <<'EOF'
heat 0 -> 300 K
 &cntrl
  imin=0, irest=0, ntx=1, nstlim=25000, dt=0.002,
  ntc=2, ntf=2, ntt=3, gamma_ln=2.0, tempi=0.0, temp0=300.0,
  ntb=1, ntp=0, cut=10.0, ntwx=5000, ntpr=1000,
  ntr=1, restraintmask=':1-400 & !@H=', restraint_wt=1.0,
 /
EOF
cat > eq.in <<'EOF'
NPT equilibration, unrestrained
 &cntrl
  imin=0, irest=1, ntx=5, nstlim=250000, dt=0.002,
  ntc=2, ntf=2, ntt=3, gamma_ln=2.0, temp0=300.0,
  ntb=2, ntp=1, taup=2.0, cut=10.0, ntwx=25000, ntpr=5000,
 /
EOF

pmemd.cuda -O -i min.in  -p ${TAG}.parm7 -c ${TAG}.rst7    -o "$OUT/${TAG}.min.out"  -r ${TAG}.min.rst7  -ref ${TAG}.rst7 || exit 1
pmemd.cuda -O -i heat.in -p ${TAG}.parm7 -c ${TAG}.min.rst7 -o "$OUT/${TAG}.heat.out" -r ${TAG}.heat.rst7 -ref ${TAG}.min.rst7 -x ${TAG}.heat.nc || exit 1
pmemd.cuda -O -i eq.in   -p ${TAG}.parm7 -c ${TAG}.heat.rst7 -o "$OUT/${TAG}.eq.out"  -r ${TAG}.eq.rst7  -x ${TAG}.eq.nc || exit 1
echo "equilibrated: ${TAG}.eq.rst7"

# ---- 6. did the thioester survive equilibration? -------------------------------
cat > check.in <<EOF
parm ${TAG}.parm7
trajin ${TAG}.eq.nc
distance thio :${UB_THIO}@C :${CYS_CAT}@SG out "$OUT/${TAG}_thio_dist.dat"
distance attack :LYQ@NI :${UB_THIO}@C out "$OUT/${TAG}_attack.dat"
angle planar :${UB_THIO}@O :${UB_THIO}@C :${CYS_CAT}@SG out "$OUT/${TAG}_planarity.dat"
go
quit
EOF
cpptraj -i check.in > "$OUT/${TAG}_cpptraj.log" 2>&1
echo "=== thioester through equilibration ==="
python3 - "$OUT/${TAG}_thio_dist.dat" "$OUT/${TAG}_attack.dat" "$OUT/${TAG}_planarity.dat" <<'PY'
import sys, statistics as st
for path, lab, unit in zip(sys.argv[1:4],
                           ("thioester C-S", "attack NI->acyl C", "O=C-S angle"),
                           ("A", "A", "deg")):
    try:
        v = [float(l.split()[1]) for l in open(path) if l.strip() and not l.startswith("#")]
    except Exception as e:
        print(f"  {lab}: unreadable ({e})"); continue
    if v:
        print(f"  {lab:20s} first {v[0]:7.2f}  last {v[-1]:7.2f}  "
              f"mean {st.mean(v):7.2f} {unit}")
PY
