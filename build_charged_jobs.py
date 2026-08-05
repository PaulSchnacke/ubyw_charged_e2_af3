#!/usr/bin/env python3
"""Build AF3 jobs for UBE2W CHARGED with ubiquitin, co-folded with a XisoK substrate.

This applies the ligand trick (see docs/LIGAND_TRICK_ADAPTED.md) to the question
rounds 2-5 could not answer: those all measured the LisoK neo-N-terminus against
a BARE catalytic cysteine, which turned out to be permissive -- AF3 docks the
modification into an empty active site whether or not the site is reactive
(AUC 0.367 over 14 sites), and MD showed the pose is not even a force-field
minimum (7.4 +- 0.3 A, 0/1200 frames <= 4 A).

Loading ubiquitin onto Cys91 changes the question from "can the amine reach a
cysteine" to "can it reach the carbonyl it must actually attack".

THE CHEMISTRY, and it changes which atoms matter
------------------------------------------------
An E2~Ub thioester is ubiquitin's C-terminal Gly76 carboxyl joined to the
catalytic cysteine's SG through a THIOESTER:

    Ub Gly76 C(=O)-S-CG(Cys91)          <- the loaded enzyme

Aminolysis then transfers ubiquitin to the substrate amine:

    substrate-NH2 + Ub-CO-S-Cys  ->  substrate-NH-CO-Ub + Cys-SH

So with Ub loaded:
  * the ELECTROPHILE is Ub Gly76 C, not Cys91 SG
  * the nucleophile is unchanged: the LisoK leucyl alpha-amine (ligand N01)
  * the QC distance to measure becomes N01 -> Ub Gly76 C

The friend's examples bond Ub Gly76 C to a lysine NZ, i.e. the PRODUCT
(isopeptide) state. Ours bonds Ub Gly76 C to Cys91 SG, i.e. the PRE-TRANSFER
thioester state. Same mechanism, different point along it -- and the
pre-transfer state is the one that can discriminate, because the substrate amine
is still free.

TWO STATES ARE BUILT, deliberately
----------------------------------
  charged   Ub thioester on Cys91, substrate amine free.  THE EXPERIMENT.
  product   Ub isopeptide-bonded to the substrate's LisoK amine, Cys91 free.
            A positive control: if AF3 cannot even build a plausible product
            complex, a null in the charged state is uninterpretable.

AF3 CAVEAT, stated plainly
--------------------------
AF3 has no thioester chemistry. A `bondedAtomPairs` entry between Gly76 C and
Cys91 SG tells it "these two atoms are bonded" and it will place them ~1.3-1.8 A
apart, but the geometry around that bond is not parameterised as a thioester and
the C=O will not necessarily be planar or correctly oriented. This is a
geometric proxy for a loaded E2, not a chemically faithful one. It is still a
much better proxy than a bare cysteine.

Usage:
    python build_charged_jobs.py OUTDIR [--ccd LisoK_userCCD.cif] [--seeds 1 2 3 4 5]
"""
import argparse
import copy
import json
import os
import sys

# ---------------------------------------------------------------- sequences
# Human ubiquitin, residues 1-76 (UBB P0CG47 / UBC P0CG48 / RPS27A P62979 all
# encode the identical 76-mer). Ends ...LRLRGG: Gly76 is the C-terminal residue
# whose carboxyl carbon forms the thioester.
UBIQUITIN = ("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKE"
             "STLHLVLRLRGG")

# Human UBE2W (UniProt Q96B02), 151 aa. Catalytic Cys91.
UBE2W = ("MASMQKRLQKELLALQNDPPPGMTLNEKSVQNSITQWIVDMEGAPGTLYEGEKFQLLFKFSSRY"
             "PFDSPQVMFTGENIPVHPHVYSNGHICLSILTEDWSPALSVQSVCLSIISMLSSCKEKRRPPDN"
             "SFYVRTCNKNPKKTKWWYHDDTC")

# SUMO2 construct as used in the paper: an N-terminal Pro is added because
# "As UBE2W is known to ubiquitylate SUMO2's native N-terminus we introduced an
# N-terminal proline residue to prohibit its modification." (Results, p3).
# That Pro shifts numbering by +1, so native K11 is construct residue 12.
SUMO2_NATIVE = ("MADEKPKEGVKTENNDHINLKVAGQDGSVVQFKIKRHTPLSKLMKAYCERQGLSMRQIRFRFDG"
                "QPINETDTPAQLEMEDEDTIDVFQQQTGGVY")
SUMO2_CONSTRUCT = "P" + SUMO2_NATIVE

UB_CTERM_RES = 76          # Gly76
UB_CTERM_ATOM = "C"        # carboxyl carbon -> thioester / isopeptide carbon
UBE2W_CAT_CYS = 91
CYS_SG = "SG"
LIG_COMP = "LIG-1"         # every XisoK user-CCD in this project uses this id
LIG_LINK = "C01"           # acyl carbonyl C, bonds to the substrate lysine NZ
LIG_NUCLEOPHILE = "N01"    # free alpha-amine: THE nucleophile

SITES = {"k11": 12, "k21": 22}    # construct numbering (native K11 -> 12)


def check_sequences():
    """Fail before writing anything if the hardcoded sequences are wrong."""
    problems = []
    if len(UBIQUITIN) != 76:
        problems.append(f"ubiquitin is {len(UBIQUITIN)} aa, expected 76")
    if not UBIQUITIN.endswith("GG"):
        problems.append("ubiquitin does not end in Gly-Gly")
    if UBIQUITIN[UB_CTERM_RES - 1] != "G":
        problems.append(f"ubiquitin residue {UB_CTERM_RES} is "
                        f"{UBIQUITIN[UB_CTERM_RES-1]}, expected G")
    if len(UBE2W) != 151:
        problems.append(f"UBE2W is {len(UBE2W)} aa, expected 151")
    if UBE2W[UBE2W_CAT_CYS - 1] != "C":
        problems.append(f"UBE2W residue {UBE2W_CAT_CYS} is "
                        f"{UBE2W[UBE2W_CAT_CYS-1]}, expected C (catalytic Cys)")
    if SUMO2_CONSTRUCT[0] != "P":
        problems.append("SUMO2 construct does not start with the blocking Pro")
    for name, pos in SITES.items():
        if SUMO2_CONSTRUCT[pos - 1] != "K":
            problems.append(f"SUMO2 construct residue {pos} ({name}) is "
                            f"{SUMO2_CONSTRUCT[pos-1]}, expected K")
    if problems:
        sys.exit("sequence checks failed:\n  " + "\n  ".join(problems))
    return True


def build(name, site_pos, ccd_text, seeds, state):
    """Assemble one AF3 payload.

    state='charged': Ub thioester on Cys91; substrate amine free (the experiment)
    state='product': Ub isopeptide on the substrate amine; Cys91 free (control)
    """
    if state not in ("charged", "product"):
        raise ValueError(state)

    sequences = [
        {"protein": {"id": "A", "sequence": SUMO2_CONSTRUCT}},   # substrate
        {"protein": {"id": "B", "sequence": UBE2W}},             # enzyme
        {"protein": {"id": "U", "sequence": UBIQUITIN}},         # the "ligand"
        {"ligand": {"id": "L", "ccdCodes": [LIG_COMP]}},         # XisoK acyl
    ]

    # bond 1, both states: the XisoK acyl group onto the substrate lysine.
    # This is the modification itself and is present regardless.
    bonds = [[["A", site_pos, "NZ"], ["L", 1, LIG_LINK]]]

    if state == "charged":
        # Ub Gly76 carboxyl C -> Cys91 SG. The thioester. Substrate amine free.
        bonds.append([["U", UB_CTERM_RES, UB_CTERM_ATOM],
                      ["B", UBE2W_CAT_CYS, CYS_SG]])
    else:
        # Ub Gly76 carboxyl C -> the LisoK free amine. The product isopeptide.
        bonds.append([["U", UB_CTERM_RES, UB_CTERM_ATOM],
                      ["L", 1, LIG_NUCLEOPHILE]])

    return {
        "name": name,
        "sequences": sequences,
        "modelSeeds": list(seeds),
        "dialect": "alphafold3",
        "version": 1,
        "bondedAtomPairs": bonds,
        "userCCD": ccd_text,
    }


def validate(payload, state, site_pos):
    """Assert the payload says what we think before a collaborator burns GPU."""
    p = []
    ids = [next(iter(s.values()))["id"] for s in payload["sequences"]]
    if len(ids) != len(set(ids)):
        p.append(f"duplicate chain ids: {ids}")

    chains = {}
    for s in payload["sequences"]:
        kind, body = next(iter(s.items()))
        chains[body["id"]] = (kind, body)

    for cid in ("A", "B", "U", "L"):
        if cid not in chains:
            p.append(f"chain {cid} missing")
    if p:
        return p

    # every bonded atom must name a chain that exists, at a residue that exists,
    # carrying the residue type the bond assumes
    for bond in payload["bondedAtomPairs"]:
        for cid, resnum, atom in bond:
            if cid not in chains:
                p.append(f"bond references absent chain {cid}")
                continue
            kind, body = chains[cid]
            if kind == "protein":
                seq = body["sequence"]
                if not 1 <= resnum <= len(seq):
                    p.append(f"{cid}:{resnum} out of range (len {len(seq)})")
                    continue
                aa = seq[resnum - 1]
                if atom == "NZ" and aa != "K":
                    p.append(f"{cid}:{resnum} is {aa}, but NZ requires Lys")
                if atom == SG_ATOM_NAME and aa != "C":
                    p.append(f"{cid}:{resnum} is {aa}, but SG requires Cys")
                if cid == "U" and resnum == UB_CTERM_RES and aa != "G":
                    p.append(f"U:{resnum} is {aa}, expected Gly76")
            else:
                if resnum != 1:
                    p.append(f"ligand {cid} residue {resnum}, expected 1")

    # exactly two bonds, and the acyl-onto-lysine bond must be one of them
    if len(payload["bondedAtomPairs"]) != 2:
        p.append(f"{len(payload['bondedAtomPairs'])} bonds, expected 2")
    mod_bond = [["A", site_pos, "NZ"], ["L", 1, LIG_LINK]]
    if mod_bond not in payload["bondedAtomPairs"]:
        p.append("the XisoK-onto-lysine bond is missing")

    # the two states must differ in exactly the intended way
    thioester = [["U", UB_CTERM_RES, UB_CTERM_ATOM],
                 ["B", UBE2W_CAT_CYS, CYS_SG]]
    isopep = [["U", UB_CTERM_RES, UB_CTERM_ATOM], ["L", 1, LIG_NUCLEOPHILE]]
    if state == "charged":
        if thioester not in payload["bondedAtomPairs"]:
            p.append("charged state lacks the Gly76-Cys91 thioester bond")
        if isopep in payload["bondedAtomPairs"]:
            p.append("charged state must NOT pre-form the product isopeptide")
    else:
        if isopep not in payload["bondedAtomPairs"]:
            p.append("product state lacks the Gly76-amine isopeptide bond")
        if thioester in payload["bondedAtomPairs"]:
            p.append("product state must NOT keep the thioester")

    if LIG_COMP not in payload.get("userCCD", ""):
        p.append(f"userCCD does not define {LIG_COMP}")
    for key in ("name", "sequences", "modelSeeds", "dialect", "version"):
        if key not in payload:
            p.append(f"missing top-level key {key}")
    return p


SG_ATOM_NAME = CYS_SG


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--ccd", default="LisoK_userCCD.cif")
    ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3, 4, 5])
    a = ap.parse_args()

    check_sequences()
    print(f"sequences verified: Ub 76 aa ending Gly76, UBE2W Cys{UBE2W_CAT_CYS}, "
          f"SUMO2 construct {len(SUMO2_CONSTRUCT)} aa with blocking Pro1")

    ccd_text = open(a.ccd).read()
    if LIG_COMP not in ccd_text:
        sys.exit(f"{a.ccd} does not define {LIG_COMP}")

    os.makedirs(a.outdir, exist_ok=True)
    manifest = []
    for site, pos in SITES.items():
        for state in ("charged", "product"):
            name = f"sumo2_{site}lisok_ube2w_ub_{state}"
            payload = build(name, pos, ccd_text, a.seeds, state)
            problems = validate(payload, state, pos)
            if problems:
                sys.exit(f"{name} failed validation:\n  " + "\n  ".join(problems))
            path = os.path.join(a.outdir, f"{name}.json")
            with open(path, "w") as fh:
                json.dump(payload, fh, indent=1)
            nres = sum(len(s["protein"]["sequence"])
                       for s in payload["sequences"] if "protein" in s)
            manifest.append(dict(name=name, file=os.path.basename(path),
                                 site=site, construct_residue=pos, state=state,
                                 n_protein_residues=nres,
                                 n_bonds=len(payload["bondedAtomPairs"]),
                                 seeds=a.seeds))
            print(f"  {name}: {nres} protein residues + 1 ligand, "
                  f"{len(payload['bondedAtomPairs'])} bonds -> {path}")

    with open(os.path.join(a.outdir, "manifest.json"), "w") as fh:
        json.dump(dict(jobs=manifest,
                       qc_atom_pairs=dict(
                           charged=f"L:1:{LIG_NUCLEOPHILE} -> U:{UB_CTERM_RES}:"
                                   f"{UB_CTERM_ATOM}  (attack on the thioester "
                                   f"carbonyl -- THE measurement)",
                           reference=f"L:1:{LIG_NUCLEOPHILE} -> B:"
                                     f"{UBE2W_CAT_CYS}:{CYS_SG}  (comparable to "
                                     f"the uncharged rounds)")), fh, indent=1)
    print(f"\nwrote {len(manifest)} jobs + manifest.json to {a.outdir}/")
    print(f"QC in the charged state: measure L:1:{LIG_NUCLEOPHILE} to "
          f"U:{UB_CTERM_RES}:{UB_CTERM_ATOM}, not to Cys91 SG")


if __name__ == "__main__":
    main()
