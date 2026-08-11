#!/bin/bash
# Parameterise LYQ (N-epsilon-leucyl-lysine = LisoK) as an ff19SB-compatible
# residue, using the standard AmberTools route for a non-standard amino acid.
#
# Why this route rather than hand-editing charges: antechamber assigns AM1-BCC
# charges to the CAPPED model compound, then prepgen excises the ACE/NME caps
# against a mainchain definition and redistributes the cap charge so the residue
# comes out at INTEGER net charge. That is what makes the residue safe to splice
# into a protein chain -- a hand-balanced charge set is where "neutral linkage"
# quietly becomes -0.03 e per copy.
#
# The amide half of LYQ (NZ, HZ, CH, OH) is chemically identical to ALY
# (N-epsilon-acetyllysine), which Amber 24 ships in ff19SB_modAA. So this is not
# a from-scratch parameterisation of unknown quality -- it is AM1-BCC on a
# compound whose reactive half has a published reference to check against.
# check_lyq_charges.py does that comparison.
#
# Usage: bash prep_lyq.sh   (needs module load amber/24)
set -eo pipefail

test -f lyq_model.pdb || { echo "run build_lyq_model.py first" >&2; exit 1; }
command -v antechamber >/dev/null || { echo "module load amber/24 first" >&2; exit 1; }

echo "=== 1. AM1-BCC charges on the capped model compound ==="
# -nc 0: the capped compound is neutral (neutral leucyl alpha-amine).
# -at amber: use Amber atom types, not GAFF -- this residue joins a protein.
# Output BOTH formats. mol2 is human-readable and what check_lyq_charges reads
# for the AM1-BCC values; the .ac file is what prepgen requires -- `prepgen -i`
# takes an ac file, NOT a mol2. Feeding it a mol2 makes prepgen fail to parse
# any atoms and emit a prep file containing only dummy atoms, WITHOUT a nonzero
# exit status, which tleap then accepts while noting only "LYQ: no atoms".
# That silent chain is what cost three submissions here.
antechamber -i lyq_model.pdb -fi pdb -o lyq_capped.mol2 -fo mol2 \
            -c bcc -nc 0 -at amber -rn LYQ -pf y -s 2 > antechamber.log 2>&1 \
  || { echo "antechamber FAILED -- see antechamber.log"; tail -25 antechamber.log; exit 1; }
antechamber -i lyq_capped.mol2 -fi mol2 -o lyq_capped.ac -fo ac \
            -at amber -rn LYQ -pf y >> antechamber.log 2>&1 \
  || { echo "mol2 -> ac conversion FAILED"; tail -20 antechamber.log; exit 1; }
test -s lyq_capped.ac || { echo "lyq_capped.ac is empty" >&2; exit 1; }
echo "  wrote lyq_capped.mol2 and lyq_capped.ac ($(grep -c '^ATOM' lyq_capped.ac) atoms)"

echo "=== 2. excise the caps, forcing integer residue charge ==="
# Mainchain definition. GOTCHA, and it fails SILENTLY: MAIN_CHAIN must list
# every backbone atom on the path from HEAD to TAIL, and OMIT_NAME must name
# every CAP atom to delete. A previous attempt used `MAIN_CHAIN CA` with
# `OMIT_NAME H1`, which told prepgen the residue was ONE atom -- it duly deleted
# everything and emitted a prep file containing only dummy atoms. tleap then
# built the test peptide with exit code 0 and merely noted "LYQ: no atoms", so
# nothing failed loudly. Hence check_lyq_charges.py, which counts atoms.
#
# For ACE-LYQ-NME the path is N -> CA -> C, and the caps are the ACE atoms
# (CH3, HH31-33, C, O of residue 1) and the NME atoms (N, H, CH3, HH31-33 of
# residue 3). prepgen matches OMIT_NAME by atom name, and antechamber has
# flattened all three residues into one LYQ residue with duplicate names, so
# the caps are identified by the names antechamber assigned them.
{
  echo "HEAD_NAME N"
  echo "TAIL_NAME C"
  # MAIN_CHAIN takes ONE LINE PER ATOM and must list the whole backbone path
  # BETWEEN head and tail. `MAIN_CHAIN CA` on its own makes prepgen report
  # "Number of mainchain atoms (including head and tail atom): 1" -- it does not
  # parse a space-separated list, and it does not infer the path. For an amino
  # acid the path from N to C is just CA, so one line is genuinely correct here;
  # the count of 1 was the symptom that the earlier run's residue had been
  # emptied by OMIT_NAME, not that MAIN_CHAIN itself was wrong.
  echo "MAIN_CHAIN CA"
  # build_lyq_model.py emits the cap atom names and has asserted that none of
  # them collides with a residue atom name
  sed 's/^/OMIT_NAME /' lyq_model_caps.txt
  echo "PRE_HEAD_TYPE C"
  echo "POST_TAIL_TYPE N"
  echo "CHARGE 0.0"
} > lyq.mc
echo "  mainchain definition:"; sed 's/^/    /' lyq.mc
prepgen -i lyq_capped.ac -o lyq.prep -m lyq.mc -rn LYQ -rf lyq.res \
        > prepgen.log 2>&1 \
  || { echo "prepgen FAILED -- see prepgen.log"; tail -25 prepgen.log; exit 1; }
# prepgen exits 0 even when it has parsed no atoms, so assert on the CONTENT.
# A prep file for this residue must contain 17 heavy atoms; an empty one contains
# only three DUMM lines.
NREAL=$(awk 'NF>=11 && $1 ~ /^[0-9]+$/ && $2 != "DUMM"' lyq.prep | wc -l)
if [ "$NREAL" -lt 17 ]; then
  echo "prepgen produced only $NREAL real atoms (expected >=17) -- the residue" >&2
  echo "was emptied. Check lyq.mc against the atom names in lyq_capped.ac." >&2
  echo "--- lyq.prep ---" >&2; cat lyq.prep >&2
  echo "--- prepgen.log ---" >&2; cat prepgen.log >&2
  exit 1
fi
echo "  wrote lyq.prep ($NREAL real atoms)"

echo "=== 3. missing bonded parameters ==="
# parmchk2 on the PREP file fails with "Atom type of DU does not exist in
# PARMCHK.DAT" -- the three DUMM placeholder atoms prepgen writes have no type.
# Run it on the AC file instead, which has no dummies, and check against
# parm19.dat since this residue joins an ff19SB protein (the GAFF default would
# report every protein-standard term as missing).
parmchk2 -i lyq_capped.ac -f ac -o lyq.frcmod -a Y \
         -p "$AMBER_EULER_ROOT/dat/leap/parm/parm19.dat" > parmchk2.log 2>&1 \
  || parmchk2 -i lyq_capped.ac -f ac -o lyq.frcmod > parmchk2.log 2>&1 \
  || { echo "parmchk2 FAILED"; cat parmchk2.log; exit 1; }
test -s lyq.frcmod || { echo "lyq.frcmod is empty" >&2; exit 1; }
echo "  wrote lyq.frcmod"
echo "  parameters flagged as guessed (ATTN):"
grep -c ATTN lyq.frcmod || echo "    none"

echo "=== 4. build a standalone tripeptide to prove tleap accepts LYQ ==="
cat > test_lyq.leap <<'LEAP'
source leaprc.protein.ff19SB
source leaprc.water.opc
loadamberparams lyq.frcmod
loadamberprep lyq.prep
# ACE-ALA-LYQ-ALA-NME: LYQ must bond into a chain on both sides
x = sequence { ACE ALA LYQ ALA NME }
saveamberparm x test_lyq.parm7 test_lyq.rst7
charge x
quit
LEAP
tleap -f test_lyq.leap > tleap_test.log 2>&1 || true
if [ -f test_lyq.parm7 ]; then
  echo "  tleap built the test peptide"
  grep -iE 'Total unperturbed charge|unperturbed' tleap_test.log | tail -2
else
  echo "  tleap FAILED to build the test peptide:"
  grep -iE 'error|fatal|could not|missing' tleap_test.log | head -15
  exit 1
fi
echo
echo "PASS: LYQ is a usable ff19SB residue. Next: check_lyq_charges.py"
