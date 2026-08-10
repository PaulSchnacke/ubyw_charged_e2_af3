#!/usr/bin/env python3
"""QC every model of the sweep, on the cluster, writing one compact CSV.

Runs remotely so nothing large has to be transferred: 13 jobs x 20 seeds x 5
samples is ~1300 CIFs, far past the harvest threshold. The CSV is a few hundred
kB and lands in a durable home directory, so it survives both the scratch purge
and a dropped VPN.

Chains are resolved FROM THE COVALE RECORDS, not from the input naming: AF3
renumbers a split chain (it absorbed a GLY ligand into the ubiquitin chain last
time, making it 75 residues rather than 74) and it renames ligand chains. Reading
the input spec would measure the wrong atoms.

Usage: python qc_sweep.py MODEL_ROOT OUT.csv
"""
import csv
import glob
import json
import os
import sys

import numpy as np
import gemmi

ATTACK = 4.0
UBE2W_CAT_CYS = 91


def ang(a, b, c):
    v1, v2 = a - b, c - b
    return np.degrees(np.arccos(np.clip(v1 @ v2 / np.linalg.norm(v1)
                                        / np.linalg.norm(v2), -1, 1)))


def load(cif):
    st = gemmi.read_structure(cif)
    st.setup_entities()
    A = {}
    for ch in st[0]:
        for res in ch:
            for at in res:
                A[(ch.name, res.seqid.num, res.name, at.name)] = \
                    np.array(at.pos.tolist())
    return st, A


def analyse(cif):
    st, A = load(cif)
    r = {}

    def get(pred):
        return [(k, v) for k, v in A.items() if pred(k)]

    lig = {k[3]: v for k, v in A.items() if k[2] == "LIG-1"}
    r["has_ligand"] = bool(lig)

    # declared bonds, measured as AF3 actually wrote them
    bonds = []
    for con in st.connections:
        p1, p2 = con.partner1, con.partner2
        k1 = (p1.chain_name, p1.res_id.seqid.num, p1.atom_name)
        k2 = (p2.chain_name, p2.res_id.seqid.num, p2.atom_name)
        v1 = next((v for k, v in A.items() if (k[0], k[1], k[3]) == k1), None)
        v2 = next((v for k, v in A.items() if (k[0], k[1], k[3]) == k2), None)
        if v1 is not None and v2 is not None:
            bonds.append((p1.atom_name, p2.atom_name, float(np.linalg.norm(v1 - v2))))
    r["n_bonds"] = len(bonds)
    for a, b, d in bonds:
        pair = "-".join(sorted([a, b]))
        r[f"bond_{pair}"] = round(d, 3)

    # the catalytic cysteine, and whether it is the NEAREST cysteine (specificity)
    cys = get(lambda k: k[2] == "CYS" and k[3] == "SG")
    cat = [v for k, v in cys if k[1] == UBE2W_CAT_CYS]

    # nucleophile: the XisoK free amine, when it is not bonded away
    nuc = lig.get("N01")
    if nuc is not None and cat:
        ds = sorted((float(np.linalg.norm(nuc - v)), f"{k[0]}/{k[1]}") for k, v in cys)
        r["n01_catcys_sg"] = round(float(np.linalg.norm(nuc - cat[0])), 3)
        r["n01_nearest_cys"] = ds[0][1]
        r["n01_nearest_cys_d"] = round(ds[0][0], 3)
        r["catcys_is_nearest"] = ds[0][1].endswith(f"/{UBE2W_CAT_CYS}")
        r["in_attack_range"] = r["n01_catcys_sg"] <= ATTACK

    # thioester carbon geometry, for the charged and tetrahedral species
    for gres in ("UBGG", "UBGT"):
        gl = {k[3]: v for k, v in A.items() if k[2] == gres}
        if not gl or "C2" not in gl or not cat:
            continue
        c, o, ca = gl["C2"], gl.get("O2"), gl.get("CA2")
        r["thio_species"] = gres
        r["thio_c_s"] = round(float(np.linalg.norm(c - cat[0])), 3)
        if o is not None and ca is not None:
            r["thio_planarity"] = round(float(ang(o, c, cat[0]) + ang(ca, c, cat[0])
                                              + ang(o, c, ca)), 1)
        r["thio_has_oxt"] = "OXT" in gl
        if nuc is not None:
            r["n01_thio_c"] = round(float(np.linalg.norm(nuc - c)), 3)

    # lysine control: no ligand, so measure the lysine NZ itself
    if not lig:
        nz = get(lambda k: k[0] == "A" and k[2] == "LYS" and k[3] == "NZ")
        if nz and cat:
            ds = sorted((float(np.linalg.norm(v - cat[0])), k[1]) for k, v in nz)
            r["nearest_lys_nz_d"] = round(ds[0][0], 3)
            r["nearest_lys_resi"] = ds[0][1]

    # substrate-enzyme interface size, as a docking sanity check
    sub = np.array([v for k, v in A.items() if k[0] == "A" and k[3] == "CA"])
    enz = np.array([v for k, v in A.items() if k[0] == "B" and k[3] == "CA"])
    if len(sub) and len(enz):
        dmin = np.sqrt(((sub[:, None, :] - enz[None, :, :]) ** 2).sum(-1)).min(1)
        r["iface_res_8A"] = int((dmin < 8).sum())

    r["chains"] = ";".join(f"{c.name}:{len(c)}" for c in st[0])
    return r


def main(root, out):
    rows = []
    for cif in sorted(glob.glob(os.path.join(root, "**", "*model*.cif"),
                                recursive=True)):
        base = os.path.dirname(cif)
        try:
            r = analyse(cif)
        except Exception as e:
            r = {"error": f"{type(e).__name__}: {e}"}
        rel = os.path.relpath(cif, root)
        parts = rel.split(os.sep)
        r["job"] = parts[0]
        r["sample"] = parts[1] if len(parts) > 2 else os.path.basename(base)
        r["model"] = rel
        cf = os.path.join(base, "summary_confidences.json")
        if os.path.exists(cf):
            try:
                c = json.load(open(cf))
                for k in ("iptm", "ptm", "ranking_score"):
                    r[k] = c.get(k)
            except Exception:
                pass
        rows.append(r)

    if not rows:
        sys.exit(f"no models found under {root}")
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} models -> {out}")

    # terse per-job summary so the log alone is informative if transfer fails
    byjob = {}
    for r in rows:
        byjob.setdefault(r["job"], []).append(r)
    for job, rs in sorted(byjob.items()):
        d = [r["n01_catcys_sg"] for r in rs if r.get("n01_catcys_sg") is not None]
        if d:
            print(f"  {job:34s} n={len(rs):4d}  amine-Cys91 min {min(d):6.2f} "
                  f"median {sorted(d)[len(d)//2]:6.2f}  "
                  f"<=4A {sum(1 for x in d if x <= ATTACK):3d}/{len(d)}")
        else:
            nz = [r["nearest_lys_nz_d"] for r in rs
                  if r.get("nearest_lys_nz_d") is not None]
            extra = f"  nearest Lys NZ min {min(nz):.2f}" if nz else ""
            print(f"  {job:34s} n={len(rs):4d}  (no free amine){extra}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
