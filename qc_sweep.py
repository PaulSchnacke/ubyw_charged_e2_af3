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

import math

# NO third-party imports. The cluster's bare python3 has no numpy and the AF3
# container has numpy but no gemmi, and an overnight run must not depend on
# resolving that. Everything here is stdlib: AF3 writes mmCIF with a fixed
# ATOM/HETATM column layout and an explicit _struct_conn block, both of which are
# straightforward to parse directly. This also removes the crash that killed the
# first sites attempt AFTER its AF3 work had already succeeded.


def dist(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def angle(a, b, c):
    v1 = [x - y for x, y in zip(a, b)]
    v2 = [x - y for x, y in zip(c, b)]
    n1 = math.sqrt(sum(x * x for x in v1))
    n2 = math.sqrt(sum(x * x for x in v2))
    if n1 == 0 or n2 == 0:
        return float("nan")
    cosv = sum(x * y for x, y in zip(v1, v2)) / (n1 * n2)
    return math.degrees(math.acos(max(-1.0, min(1.0, cosv))))


def read_cif(path):
    """Return (atoms, connections) from an AF3 mmCIF using stdlib only.

    atoms: {(chain, resnum, resname, atomname): (x, y, z)}
    connections: [(chain1, res1, atom1, chain2, res2, atom2)] from _struct_conn.
    Handles both the looped and single-value forms of _struct_conn, since AF3
    writes a loop when there are several bonds and bare key/value when there is
    one.
    """
    atoms, conns = {}, []
    lines = open(path).read().split("\n")
    # --- atoms: the atom_site loop has a fixed AF3 column order
    hdr = []
    in_atom_loop = False
    for ln in lines:
        if ln.startswith("_atom_site."):
            hdr.append(ln.strip().split(".")[1])
            in_atom_loop = True
            continue
        if in_atom_loop and (ln.startswith("ATOM") or ln.startswith("HETATM")):
            f = ln.split()
            if len(f) < len(hdr):
                continue
            rec = dict(zip(hdr, f))
            try:
                key = (rec["label_asym_id"], int(rec["label_seq_id"])
                       if rec["label_seq_id"] not in (".", "?")
                       else int(rec["auth_seq_id"]),
                       rec["label_comp_id"], rec["label_atom_id"])
                atoms[key] = (float(rec["Cartn_x"]), float(rec["Cartn_y"]),
                              float(rec["Cartn_z"]))
            except (KeyError, ValueError):
                continue
        elif in_atom_loop and ln.strip() in ("#", "") and atoms:
            in_atom_loop = False
    # --- connections. Two AF3 forms: bare key/value when there is ONE bond, and
    # a loop_ when there are several. Ligand partners carry label_seq_id "." so
    # the auth_seq_id must be used as the fallback -- an int() on "." was why an
    # earlier version of this parser silently reported zero bonds.
    def seqid(rec, side):
        for k in (f"ptnr{side}_label_seq_id", f"ptnr{side}_auth_seq_id"):
            v = rec.get(k, ".")
            if v not in (".", "?", ""):
                try:
                    return int(v)
                except ValueError:
                    pass
        return None

    def grab(rec):
        try:
            c1 = rec["ptnr1_label_asym_id"]; a1 = rec["ptnr1_label_atom_id"]
            c2 = rec["ptnr2_label_asym_id"]; a2 = rec["ptnr2_label_atom_id"]
        except KeyError:
            return None
        r1, r2 = seqid(rec, 1), seqid(rec, 2)
        if r1 is None or r2 is None:
            return None
        return (c1, r1, a1, c2, r2, a2)

    sc_hdr, single, in_loop = [], {}, False
    for ln in lines:
        s = ln.strip()
        if s.startswith("_struct_conn."):
            parts = s.split(None, 1)
            k = parts[0].split(".", 1)[1]
            if k not in sc_hdr:
                sc_hdr.append(k)
            if len(parts) > 1 and parts[1].strip():
                single[k] = parts[1].strip()   # bare key/value form
            else:
                in_loop = True                 # header of a loop_
            continue
        if in_loop and sc_hdr:
            if s in ("#", ""):
                in_loop = False
                continue
            # a data row of the loop: same field count as the header
            f = s.split()
            if len(f) == len(sc_hdr):
                g = grab(dict(zip(sc_hdr, f)))
                if g:
                    conns.append(g)
    if single:
        g = grab(single)
        if g:
            conns.append(g)
    return atoms, conns

ATTACK = 4.0
UBE2W_CAT_CYS = 91


def analyse(cif):
    A, conns = read_cif(cif)
    r = {}
    lig = {k[3]: v for k, v in A.items() if k[2] == "LIG-1"}
    r["has_ligand"] = bool(lig)

    # declared bonds, measured as AF3 actually wrote them (chains resolved from
    # the covale records, NOT the input spec -- AF3 renumbers split chains)
    r["n_bonds"] = len(conns)
    for c1, r1, a1, c2, r2, a2 in conns:
        v1 = next((v for k, v in A.items() if (k[0], k[1], k[3]) == (c1, r1, a1)), None)
        v2 = next((v for k, v in A.items() if (k[0], k[1], k[3]) == (c2, r2, a2)), None)
        if v1 and v2:
            r["bond_" + "-".join(sorted([a1, a2]))] = round(dist(v1, v2), 3)

    cys = [(k, v) for k, v in A.items() if k[2] == "CYS" and k[3] == "SG"]
    cat = [v for k, v in cys if k[1] == UBE2W_CAT_CYS]
    nuc = lig.get("N01")

    if nuc and cat:
        ds = sorted((dist(nuc, v), f"{k[0]}/{k[1]}") for k, v in cys)
        r["n01_catcys_sg"] = round(dist(nuc, cat[0]), 3)
        r["n01_nearest_cys"] = ds[0][1]
        r["n01_nearest_cys_d"] = round(ds[0][0], 3)
        r["catcys_is_nearest"] = ds[0][1].endswith("/%d" % UBE2W_CAT_CYS)
        r["in_attack_range"] = r["n01_catcys_sg"] <= ATTACK

    # thioester carbon geometry for the charged / tetrahedral species
    for gres in ("UBGG", "UBGT"):
        gl = {k[3]: v for k, v in A.items() if k[2] == gres}
        if not gl or "C2" not in gl or not cat:
            continue
        c, o, ca = gl["C2"], gl.get("O2"), gl.get("CA2")
        r["thio_species"] = gres
        r["thio_c_s"] = round(dist(c, cat[0]), 3)
        if o and ca:
            r["thio_planarity"] = round(angle(o, c, cat[0]) + angle(ca, c, cat[0])
                                        + angle(o, c, ca), 1)
        r["thio_has_oxt"] = "OXT" in gl
        if nuc:
            r["n01_thio_c"] = round(dist(nuc, c), 3)

    # lysine control: no ligand, so measure the lysine NZ itself
    if not lig and cat:
        nz = [(k, v) for k, v in A.items()
              if k[0] == "A" and k[2] == "LYS" and k[3] == "NZ"]
        if nz:
            ds = sorted((dist(v, cat[0]), k[1]) for k, v in nz)
            r["nearest_lys_nz_d"] = round(ds[0][0], 3)
            r["nearest_lys_resi"] = ds[0][1]

    # interface size, as a docking sanity check
    sub = [v for k, v in A.items() if k[0] == "A" and k[3] == "CA"]
    enz = [v for k, v in A.items() if k[0] == "B" and k[3] == "CA"]
    if sub and enz:
        r["iface_res_8A"] = sum(1 for s in sub if min(dist(s, e) for e in enz) < 8)

    chains = {}
    for k in A:
        chains.setdefault(k[0], set()).add(k[1])
    r["chains"] = ";".join(f"{c}:{len(v)}" for c, v in sorted(chains.items()))
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
