#!/usr/bin/env python3
"""Compare UBE2W charged with Ub(1-76) against Ub(1-75), and the non-covalent pair.

THE QUESTION. Paul engineered UBE2W to accept Ub(1-75) as well as Ub(1-76);
experimentally it failed, and the failure could sit at UBA1 (charging) or at UBE2W
(accommodating the shorter donor). If UBE2W cannot reach the shorter ubiquitin's
C-terminal carbonyl, that is a structural explanation on the UBE2W side.

WHAT IS MEASURED, and why each one:

  1. thioester distance      acyl C -> Cys91 SG. In the charged jobs this bond is
                             DECLARED, so the interesting number is whether AF3 kept
                             it (it silently drops polymer-polymer bonds) and what
                             geometry it built around it.
  2. attack distance         XisoK free amine (N01) -> acyl C. The reaction
                             coordinate, for the charged jobs only.
  3. reach (non-covalent)    Ub C-terminal carbonyl C -> Cys91 SG with NO bond
                             imposed. This is the honest version of the question:
                             can the shorter tail get there at all?
  4. clashes                 heavy-atom pairs < 2.8 A between the Ub C-terminal
                             segment and UBE2W, excluding the bonded partners. A
                             shorter tail may be forced into the enzyme.
  5. per-residue contacts    which UBE2W residues sit within 5 A of the Ub tail, per
                             variant. The direct answer to "are there particular
                             residues that do not permit the shorter monomer" --
                             residues contacted in 1-75 but not 1-76 are candidates.
  6. confidence              pLDDT for the Ub tail and Cys91 region, plus interface
                             PAE and the summary iptm/ptm. AF3 confidence is NOT a
                             reactivity predictor (established over four rounds), so
                             these are reported as model quality, not as evidence
                             about the engineering.

STDLIB ONLY. The cluster's bare python3 has neither numpy nor gemmi, and a QC step
importing them once crashed AFTER AF3 had written 101 models.

Usage: python analyse_ub75.py MODEL_DIR OUT_PREFIX
"""
import csv
import glob
import json
import math
import os
import statistics as st
import sys

CATALYTIC_CYS = 91
ENZ, UB, SUB, LIG = "B", "U", "A", "L"
CLASH = 2.8          # heavy-atom contact closer than this is a clash
CONTACT = 5.0        # residue-level contact cutoff
TAIL_FROM = 70       # "C-terminal segment" of ubiquitin


def read_cif(path):
    """Minimal mmCIF atom_site + struct_conn reader. Returns (atoms, conns)."""
    atoms, conns = [], []
    hdr, inloop = [], False
    cn_hdr, cn_loop = [], False
    for ln in open(path, errors="ignore"):
        s = ln.rstrip("\n")
        if s.startswith("_atom_site."):
            hdr.append(s.strip().split(".", 1)[1]); inloop = True; continue
        if inloop and s.startswith(("ATOM", "HETATM")):
            f = s.split()
            if len(f) < len(hdr):
                continue
            r = dict(zip(hdr, f))
            try:
                atoms.append(dict(
                    chain=r.get("auth_asym_id") or r.get("label_asym_id"),
                    resi=int(r.get("auth_seq_id") or r.get("label_seq_id")),
                    resn=r.get("label_comp_id"), name=r.get("label_atom_id"),
                    el=r.get("type_symbol"),
                    x=float(r["Cartn_x"]), y=float(r["Cartn_y"]), z=float(r["Cartn_z"]),
                    b=float(r.get("B_iso_or_equiv") or 0.0)))
            except (ValueError, KeyError):
                continue
            continue
        if inloop and atoms and not s.startswith(("ATOM", "HETATM")):
            inloop = False
        if s.startswith("_struct_conn."):
            cn_hdr.append(s.strip().split(".", 1)[1]); cn_loop = True; continue
        if cn_loop and s.startswith("covale"):
            f = s.split()
            if len(f) >= len(cn_hdr):
                conns.append(dict(zip(cn_hdr, f)))
    return atoms, conns


def d(a, b):
    return math.dist((a["x"], a["y"], a["z"]), (b["x"], b["y"], b["z"]))


def get(atoms, chain=None, resi=None, name=None, resn=None, heavy=False):
    out = atoms
    if chain is not None:
        out = [a for a in out if a["chain"] == chain]
    if resi is not None:
        out = [a for a in out if a["resi"] == resi]
    if name is not None:
        out = [a for a in out if a["name"] == name]
    if resn is not None:
        out = [a for a in out if a["resn"] == resn]
    if heavy:
        out = [a for a in out if a["el"] != "H"]
    return out


def analyse_model(path, variant):
    atoms, conns = read_cif(path)
    if not atoms:
        return None
    r = dict(model=os.path.basename(path), variant=variant, n_atoms=len(atoms),
             n_covale=len(conns))
    chains = sorted({a["chain"] for a in atoms})
    r["chains"] = "/".join(chains)

    sg = get(atoms, chain=ENZ, resi=CATALYTIC_CYS, name="SG")
    sg = sg[0] if sg else None

    # The acyl carbon. In the CHARGED jobs it is the last glycine LIGAND's C; in the
    # NON-COVALENT jobs it is the C-terminal residue of the ubiquitin protein chain.
    ub = get(atoms, chain=UB)
    ub_last = max((a["resi"] for a in ub), default=None)
    # THREE possible encodings of the ubiquitin C-terminal glycines, all seen in this
    # project's own files. Getting this wrong silently reports the WRONG atom's
    # distance: on a test file the Arg74 fallback gave 6.80 A where the real
    # thioester was 1.68 A, so each form is detected explicitly and the choice is
    # recorded in acyl_source for every row.
    #   (a) two GLY LIGAND chains (current build_jobs_v2 output): acyl C is the C of
    #       the last glycine ligand chain.
    #   (b) one UBGG dipeptide ligand (older CCD): atoms N1/CA1/C1/O1 + N2/CA2/C2/O2,
    #       and the acyl carbon is C2. If only Gly75 is present it is C1.
    #   (c) no ligand at all (non-covalent jobs): acyl C is the C of the last residue
    #       of the ubiquitin protein chain.
    gly_lig = [a for a in atoms if a["resn"] == "GLY" and a["chain"] not in (UB, ENZ, SUB)]
    ubgg = [a for a in atoms if a["resn"] in ("UBGG", "UBGT")]
    if gly_lig:
        last_ch = sorted({a["chain"] for a in gly_lig})[-1]
        acyl = get(gly_lig, chain=last_ch, name="C")
        r["acyl_source"] = f"GLY ligand {last_ch}"
        r["ub_len"] = (ub_last or 0) + len({a["chain"] for a in gly_lig})
    elif ubgg:
        names = {a["name"] for a in ubgg}
        which = "C2" if "C2" in names else "C1"
        acyl = [a for a in ubgg if a["name"] == which]
        n_gly = 2 if "C2" in names else 1
        r["acyl_source"] = f"{ubgg[0]['resn']} {which}"
        r["ub_len"] = (ub_last or 0) + n_gly
        gly_lig = ubgg          # treat as the tail for clash/contact purposes
    else:
        acyl = get(atoms, chain=UB, resi=ub_last, name="C")
        r["acyl_source"] = f"protein U{ub_last}"
        r["ub_len"] = ub_last
    acyl = acyl[0] if acyl else None
    if acyl is None:
        r["acyl_source"] += " (NOT FOUND)"

    r["thio_c_sg"] = round(d(acyl, sg), 3) if (acyl and sg) else None
    r["thio_formed"] = (r["thio_c_sg"] is not None and r["thio_c_sg"] < 2.2)

    # attack coordinate, charged jobs only (the XisoK ligand carries N01)
    n01 = get(atoms, name="N01")
    r["attack_n01_acyl"] = round(d(n01[0], acyl), 3) if (n01 and acyl) else None
    r["attack_n01_sg"] = round(d(n01[0], sg), 3) if (n01 and sg) else None

    # clashes between the ubiquitin C-terminal segment and UBE2W
    tail = [a for a in atoms
            if a["chain"] == UB and a["resi"] >= TAIL_FROM and a["el"] != "H"]
    tail += [a for a in gly_lig if a["el"] != "H"]
    enz = get(atoms, chain=ENZ, heavy=True)
    bonded = {id(acyl), id(sg)} if (acyl and sg) else set()
    clashes, contacts = 0, {}
    for t in tail:
        if id(t) in bonded:
            continue
        for e in enz:
            if id(e) in bonded:
                continue
            dd = d(t, e)
            if dd < CLASH:
                clashes += 1
            if dd < CONTACT:
                key = f"{e['resn']}{e['resi']}"
                contacts[key] = min(contacts.get(key, 99.0), round(dd, 2))
    r["n_clash_lt2.8"] = clashes
    r["n_contact_res"] = len(contacts)
    r["_contacts"] = contacts

    # confidence: pLDDT is written into the B-factor column by AF3
    r["plddt_ub_tail"] = round(st.mean([a["b"] for a in tail]), 2) if tail else None
    cys_env = [a for a in atoms if a["chain"] == ENZ
               and abs(a["resi"] - CATALYTIC_CYS) <= 3 and a["el"] != "H"]
    r["plddt_cys91_env"] = round(st.mean([a["b"] for a in cys_env]), 2) if cys_env else None
    r["plddt_all"] = round(st.mean([a["b"] for a in atoms if a["el"] != "H"]), 2)
    return r


def load_conf(model_path):
    """Match a model CIF to its summary_confidences JSON.

    Two layouts, both of which occur:
      * AF3's own tree      <job>/seed-N_sample-M/{model.cif, summary_confidences.json}
      * the flattened harvest  <job>__seed-N_sample-M__model.cif  and
                               <job>__seed-N_sample-M__summary_confidences.json
    Note the flattened form replaces "model.cif" WITHOUT an extra underscore, so
    stripping "_model.cif" and re-appending "_summary..." misses by one character --
    which silently produced empty iptm/ptm columns on the first run.
    """
    cands = [
        model_path.replace("__model.cif", "__summary_confidences.json"),
        model_path.replace("_model.cif", "_summary_confidences.json"),
        os.path.join(os.path.dirname(model_path), "summary_confidences.json"),
    ]
    for cand in cands:
        if cand != model_path and os.path.exists(cand):
            try:
                return json.load(open(cand))
            except Exception:
                return {}
    return {}


def main(model_dir, out_prefix):
    VARIANTS = ("ube2w_ub76_charged", "ube2w_ub75_charged",
                "ube2w_ub76_free", "ube2w_ub75_free")
    rows, contact_rows = [], []
    for cif in sorted(glob.glob(os.path.join(model_dir, "**", "*.cif"), recursive=True)):
        variant = next((v for v in VARIANTS if v in os.path.basename(cif)
                        or v in cif), None)
        if variant is None:
            continue
        r = analyse_model(cif, variant)
        if r is None:
            continue
        conf = load_conf(cif)
        r["iptm"] = conf.get("iptm")
        r["ptm"] = conf.get("ptm")
        r["ranking_score"] = conf.get("ranking_score")
        cn = r.pop("_contacts")
        for res, dist in cn.items():
            contact_rows.append(dict(variant=variant, model=r["model"],
                                     ube2w_res=res, min_dist=dist))
        rows.append(r)

    if not rows:
        sys.exit(f"no models found under {model_dir}")

    fields = [k for k in rows[0] if not k.startswith("_")]
    with open(f"{out_prefix}_models.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in fields})
    with open(f"{out_prefix}_contacts.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["variant", "model", "ube2w_res", "min_dist"])
        w.writeheader()
        w.writerows(contact_rows)

    print(f"{len(rows)} models, {len(contact_rows)} contact records\n")
    by = {}
    for r in rows:
        by.setdefault(r["variant"], []).append(r)

    def agg(rs, key):
        v = [x[key] for x in rs if x.get(key) is not None]
        return (min(v), st.median(v), len(v)) if v else (None, None, 0)

    print(f"{'variant':22s} {'n':>3s} {'Ub':>3s} {'thio min':>9s} {'thio med':>9s} "
          f"{'formed':>7s} {'clash med':>9s} {'pLDDT tail':>10s} {'iptm med':>8s}")
    for v in VARIANTS:
        rs = by.get(v)
        if not rs:
            print(f"{v:22s}  -- no models --")
            continue
        tmin, tmed, _ = agg(rs, "thio_c_sg")
        formed = sum(1 for x in rs if x["thio_formed"])
        _, cmed, _ = agg(rs, "n_clash_lt2.8")
        _, pmed, _ = agg(rs, "plddt_ub_tail")
        _, imed, _ = agg(rs, "iptm")
        ublen = rs[0]["ub_len"]
        print(f"{v:22s} {len(rs):3d} {ublen:3d} "
              f"{(f'{tmin:9.2f}' if tmin is not None else '       --')} "
              f"{(f'{tmed:9.2f}' if tmed is not None else '       --')} "
              f"{formed:4d}/{len(rs):<3d}"
              f"{(f'{cmed:9.1f}' if cmed is not None else '       --')} "
              f"{(f'{pmed:10.1f}' if pmed is not None else '        --')} "
              f"{(f'{imed:8.3f}' if imed is not None else '      --')}")

    # THE ANSWER TO THE ENGINEERING QUESTION: UBE2W residues contacted by the SHORT
    # tail but not the long one, and vice versa. Restricted to residues seen in at
    # least a quarter of that variant's models so single-model noise is excluded.
    print("\nUBE2W residues differentially contacted by the ubiquitin C-terminus:")
    for pair in (("ube2w_ub76_charged", "ube2w_ub75_charged"),
                 ("ube2w_ub76_free", "ube2w_ub75_free")):
        long_v, short_v = pair
        n_long = len(by.get(long_v, []))
        n_short = len(by.get(short_v, []))
        if not (n_long and n_short):
            continue
        def freq(v, n):
            c = {}
            for row in contact_rows:
                if row["variant"] == v:
                    c[row["ube2w_res"]] = c.get(row["ube2w_res"], 0) + 1
            return {k: n_ / n for k, n_ in c.items() if n_ / n >= 0.25}
        fl, fs = freq(long_v, n_long), freq(short_v, n_short)
        only_short = sorted(set(fs) - set(fl), key=lambda k: -fs[k])
        only_long = sorted(set(fl) - set(fs), key=lambda k: -fl[k])
        print(f"  {long_v} (n={n_long}) vs {short_v} (n={n_short}):")
        print(f"    contacted ONLY by Ub-short: {only_short[:12] or 'none'}")
        print(f"    contacted ONLY by Ub-long : {only_long[:12] or 'none'}")
    print("\nCaveat: AF3 confidence has been anti-predictive for reactivity across four")
    print("rounds in this project, so iptm/pLDDT here describe model quality only.")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
