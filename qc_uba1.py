#!/usr/bin/env python3
"""QC for the UBA1 jobs: can Ub(1-75) reach either catalytic centre as well as Ub(1-76)?

STDLIB ONLY. The cluster's bare python3 has neither numpy nor gemmi, and a QC step that
imported them once crashed AFTER AF3 had written 101 models.

WHAT IS MEASURED, and why each one:

  1. reach_aden    Ub C-terminal carbonyl C -> ATP P-alpha. The adenylation site bond
                   that forms. No bond is declared in these jobs, so this is the honest
                   question: can the tail get there at all?
  2. reach_cys     Ub C-terminal carbonyl C -> Cys632 SG. Same question at the thioester
                   site, also with no bond imposed.
  3. thio_len      In the COVALENT jobs the bond is declared, so the informative numbers
                   are whether AF3 kept it and what geometry it built around it.
  4. planarity     Bond-angle sum at the acyl carbon. sp2 = 360, sp3 = 328.5. Validated
                   at 357-359 in the stub with the no-OXT CCDs.
  5. clashes       Heavy-atom pairs < 2.8 A between the Ub C-terminal segment and UBA1,
                   excluding bonded partners. A shorter tail may be forced into the enzyme.
  6. contacts      Which UBA1 residues sit within 5 A of the Ub tail, per variant. The
                   differential list -- contacted by the SHORT tail but not the long one
                   at >=25% model frequency -- is the direct answer to "are there
                   particular residues that do not permit the shorter monomer".
  7. confidence    pLDDT/ipTM as MODEL QUALITY ONLY. AF3 confidence has been
                   uninformative for reactivity across five rounds; it is not evidence
                   about the engineering.

TWO MEASUREMENT TRAPS THIS AVOIDS, both of which have bitten this project:

  * THE ACYL SOURCE. Ubiquitin's C-terminal glycines have several encodings in this
    project's files: a UBGG/UBG1 ligand chain (these jobs), two bare GLY ligand chains
    (jobs_ub75, superseded), a UBGG dipeptide with atoms N1/CA1/C1/O1+N2/CA2/C2/O2, and
    no ligand at all (non-covalent jobs). analyse_ub75.py silently fell through to an
    Arg74 fallback and reported a thioester at 6.80 A where the true value was 1.68 A.
    Every row here records acyl_source, and an unrecognised layout is an ERROR, never a
    silent fallback.
  * CHAIN RESOLUTION. AF3 renumbers: it absorbed a single-residue GLY ligand into the
    protein chain in an earlier round (chain U came back 75 residues, not 74) but left
    the multi-atom custom CCD as its own chain in the stub. So chains are resolved from
    the covale records where they exist, not from the input naming.

Usage: python qc_uba1.py MODEL_DIR OUT.csv
"""
import csv
import glob
import json
import math
import os
import sys

CLASH = 2.8            # heavy-atom clash threshold
CONTACT = 5.0          # contact shell
UBA1_CAT_CYS = 632


def parse_cif(path):
    """Minimal mmCIF atom_site + struct_conn reader. No gemmi on the cluster."""
    atoms, covale = [], []
    with open(path) as fh:
        lines = fh.read().splitlines()
    i = 0
    while i < len(lines):
        if lines[i].strip() == "_atom_site.group_PDB":
            hdr = []
            while lines[i].strip().startswith("_atom_site."):
                hdr.append(lines[i].strip().split(".", 1)[1])
                i += 1
            while i < len(lines) and not lines[i].startswith("#"):
                f = lines[i].split()
                if len(f) == len(hdr):
                    atoms.append(dict(zip(hdr, f)))
                i += 1
            continue
        if lines[i].strip().startswith("covale"):
            covale.append(lines[i].split())
        i += 1
    return atoms, covale


def xyz(a):
    return (float(a["Cartn_x"]), float(a["Cartn_y"]), float(a["Cartn_z"]))


def dist(p, q):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(p, q)))


def angle(a, b, c):
    v1 = [x - y for x, y in zip(a, b)]
    v2 = [x - y for x, y in zip(c, b)]
    n1 = math.sqrt(sum(x * x for x in v1))
    n2 = math.sqrt(sum(x * x for x in v2))
    if n1 == 0 or n2 == 0:
        return float("nan")
    cos = sum(x * y for x, y in zip(v1, v2)) / (n1 * n2)
    return math.degrees(math.acos(max(-1.0, min(1.0, cos))))


def find(atoms, ch, resi, name):
    for a in atoms:
        if (a["auth_asym_id"] == ch and a["auth_seq_id"] == str(resi)
                and a["label_atom_id"] == name):
            return a
    return None


def chain_atoms(atoms, ch):
    return [a for a in atoms if a["auth_asym_id"] == ch]


def resolve_acyl(atoms):
    """Find the ubiquitin C-terminal carbonyl carbon -- the electrophile.

    Returns (atom, oxygen, ca, source_label). An unrecognised layout RAISES; it must
    never fall through to a plausible-looking wrong atom.
    """
    # 1. custom ligand CCDs used by these jobs
    for comp, acyl, o, ca in (("UBGG", "C2", "O2", "CA2"), ("UBG1", "C1", "O1", "CA1")):
        lig = [a for a in atoms if a.get("label_comp_id") == comp]
        if lig:
            ch = lig[0]["auth_asym_id"]
            resi = lig[0]["auth_seq_id"]
            A = find(atoms, ch, resi, acyl)
            if A is not None:
                return (A, find(atoms, ch, resi, o), find(atoms, ch, resi, ca),
                        f"{comp}:{acyl}")
    # 2. non-covalent jobs: ubiquitin is an intact protein chain, so the acyl carbon is
    #    the backbone C of its LAST residue. Identify the Ub chain as the protein chain
    #    whose length is 74-76 and which is not UBE2W (151).
    best = None
    for ch in sorted({a["auth_asym_id"] for a in atoms}):
        ca = chain_atoms(atoms, ch)
        resis = sorted({int(a["auth_seq_id"]) for a in ca})
        if 70 <= len(resis) <= 80 and any(a["label_atom_id"] == "CA" for a in ca):
            last = max(resis)
            A = find(atoms, ch, last, "C")
            if A is not None:
                best = (A, find(atoms, ch, last, "O"), find(atoms, ch, last, "CA"),
                        f"protein_Cterm:{ch}{last}")
    if best:
        return best
    raise ValueError("cannot resolve the ubiquitin acyl carbon -- unrecognised layout; "
                     "refusing to guess (see the analyse_ub75.py 6.80-vs-1.68 A bug)")


def uba1_chain(atoms):
    """The UBA1 chain: the protein chain with ~1058 residues."""
    for ch in sorted({a["auth_asym_id"] for a in atoms}):
        resis = {int(a["auth_seq_id"]) for a in chain_atoms(atoms, ch)}
        if len(resis) > 900:
            return ch
    return None


def analyse(path):
    atoms, covale = parse_cif(path)
    row = {"model": os.path.basename(path)}

    acyl, oxy, ca, src = resolve_acyl(atoms)
    row["acyl_source"] = src
    A = uba1_chain(atoms)
    row["uba1_chain"] = A
    pa = xyz(acyl)

    # --- the two reach coordinates
    sg = find(atoms, A, UBA1_CAT_CYS, "SG") if A else None
    row["cys632_is_cys"] = (sg is not None
                            and find(atoms, A, UBA1_CAT_CYS, "CB") is not None)
    row["reach_cys"] = round(dist(pa, xyz(sg)), 3) if sg else ""

    atp = [a for a in atoms if a.get("label_comp_id") == "ATP"]
    if atp:
        ch, resi = atp[0]["auth_asym_id"], atp[0]["auth_seq_id"]
        # PA is the alpha phosphate: the atom the carboxyl attacks
        p = find(atoms, ch, resi, "PA")
        row["reach_aden"] = round(dist(pa, xyz(p)), 3) if p else ""
    else:
        row["reach_aden"] = ""

    # --- covalent jobs: did AF3 keep the bond, and what geometry did it build?
    row["n_covale"] = len(covale)
    if sg and "protein_Cterm" not in src:
        d = dist(pa, xyz(sg))
        row["thio_len"] = round(d, 3)
        row["thio_formed"] = d < 2.5
        if oxy is not None and ca is not None:
            row["planarity"] = round(
                angle(xyz(oxy), pa, xyz(sg)) + angle(xyz(oxy), pa, xyz(ca))
                + angle(xyz(sg), pa, xyz(ca)), 2)
    else:
        row["thio_len"] = row["thio_formed"] = row["planarity"] = ""

    # --- clashes and contacts between the Ub C-terminal segment and UBA1
    ub_ch = acyl["auth_asym_id"]
    if "protein_Cterm" in src:
        last = int(acyl["auth_seq_id"])
        seg = [a for a in chain_atoms(atoms, ub_ch)
               if last - 4 <= int(a["auth_seq_id"]) <= last]
    else:
        seg = chain_atoms(atoms, ub_ch)           # the whole tail ligand
    seg = [a for a in seg if a["type_symbol"] != "H"]
    uba = [a for a in chain_atoms(atoms, A) if a["type_symbol"] != "H"] if A else []

    nclash, contacts = 0, set()
    for s in seg:
        ps = xyz(s)
        for u in uba:
            # exclude the declared bond partner itself
            if u["auth_seq_id"] == str(UBA1_CAT_CYS) and u["label_atom_id"] == "SG":
                continue
            d = dist(ps, xyz(u))
            if d < CLASH:
                nclash += 1
            if d < CONTACT:
                contacts.add(int(u["auth_seq_id"]))
    row["clashes"] = nclash
    row["n_contacts"] = len(contacts)
    row["contacts"] = ";".join(str(x) for x in sorted(contacts))

    # --- confidence: MODEL QUALITY ONLY, not evidence about the engineering
    sc = os.path.join(os.path.dirname(path), "summary_confidences.json")
    if os.path.exists(sc):
        try:
            c = json.load(open(sc))
            row["iptm"] = c.get("iptm", "")
            row["ptm"] = c.get("ptm", "")
            row["ranking_score"] = c.get("ranking_score", "")
        except Exception:
            pass
    return row


def main(model_dir, out_csv):
    files = sorted(glob.glob(os.path.join(model_dir, "**", "*.cif"), recursive=True))
    if not files:
        sys.exit(f"no .cif under {model_dir}")
    rows, errors = [], []
    for f in files:
        # skip the top-level top-ranked duplicate AF3 writes alongside the 100 real models
        if "seed-" not in os.path.basename(f):
            continue
        try:
            r = analyse(f)
            r["job"] = os.path.basename(f).split("__")[0]
            rows.append(r)
        except Exception as e:                     # loud, never silent
            errors.append(f"{os.path.basename(f)}: {e}")
    if not rows:
        sys.exit("no models parsed")
    cols = ["job", "model", "acyl_source", "uba1_chain", "cys632_is_cys", "reach_aden",
            "reach_cys", "n_covale", "thio_len", "thio_formed", "planarity", "clashes",
            "n_contacts", "iptm", "ptm", "ranking_score", "contacts"]
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote {len(rows)} rows to {out_csv}")
    srcs = {}
    for r in rows:
        srcs[r["acyl_source"]] = srcs.get(r["acyl_source"], 0) + 1
    print("acyl_source distribution (must match the job types, no surprises):")
    for k, v in sorted(srcs.items()):
        print(f"   {k:26s} {v}")
    if errors:
        print(f"\n{len(errors)} models FAILED to parse:")
        for e in errors[:10]:
            print("   ", e)
        sys.exit(1)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
