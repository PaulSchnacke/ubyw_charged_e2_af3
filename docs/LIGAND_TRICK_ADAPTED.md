# The ligand trick, and how these jobs adapt it

Credit for the pattern goes to Paul's collaborator, whose note and notebook
described it. This records what was taken unchanged and what had to change, so
the collaborator can see at a glance whether we've broken their convention.

## What the trick actually is

Despite the name, it does **not** use AF3 `ligand` entries for ubiquitin. From
the original note:

> the trick is not to emit AlphaFold ligand entries (`{"ligand": ...}`) in
> `sequences`. Instead, the workflow builds ubiquitin topologies as protein
> chains, and emits covalent links via `bondedAtomPairs`.

That is the insight worth having. Ubiquitin is 76 residues of well-conserved
sequence with thousands of PDB entries and a deep MSA. Declaring it a protein
chain gets all of that; declaring it a ligand throws it away. The covalent
attachment is then just a `bondedAtomPairs` entry, which AF3 honours regardless of
whether the partners are protein or ligand.

## Kept unchanged

* ubiquitin as a **protein chain**, residues 1–76, ending `...LRLRGG`
* the attachment expressed through **`bondedAtomPairs`**, not a ligand entry
* the electrophilic atom is ubiquitin's **C-terminal Gly76 `C`** — the carboxyl
  carbon
* the payload shape: `name`, `sequences`, `modelSeeds`, `dialect`,
  `bondedAtomPairs`, `version`

Verified against the supplied fixtures: our ubiquitin chain is byte-identical to
the sequence in `protein_bound_ubiquitin_monomer.json`, and our bond entries use
the same `[[chain, resnum, atom], [chain, resnum, atom]]` nesting.

## Changed, and why

### 1. The acceptor atom: `NZ` → `SG`

The supplied examples bond Gly76 `C` to a lysine **`NZ`**:

```json
[["A", 76, "C"], ["B", 48, "NZ"]]
```

`NZ` is the ε-amine of a lysine side chain, so that models the **product** —
ubiquitin already transferred onto a substrate lysine.

Our charged jobs bond Gly76 `C` to **`SG`** of UBE2W Cys91:

```json
[["U", 76, "C"], ["B", 91, "SG"]]
```

`SG` is the cysteine sulfur, so that models the **thioester on the enzyme**, i.e.
the state *before* transfer. This is the point of the whole exercise: in the
pre-transfer state the substrate amine is still free, so whether it can reach the
thioester carbon is a real question with a measurable answer. In the product state
the bond is already declared and nothing is being tested.

We do also build product-state jobs, using the friend's `NZ`-style pattern
adapted to our ligand's amine (`L:1:N01` rather than a protein `NZ`, since our
acceptor amine lives on the XisoK ligand) — but as a positive control, not as the
experiment.

### 2. A third chain, and a real ligand entry alongside

The examples have two or three protein chains. Ours have three proteins **plus**
one genuine AF3 `ligand` entry:

| chain | what |
|---|---|
| `A` | SUMO2 construct (96 aa, N-terminal Pro) — the substrate |
| `B` | UBE2W (151 aa) — the enzyme |
| `U` | ubiquitin (76 aa) — via the trick |
| `L` | `LIG-1`, the LisoK acyl group — a real `userCCD` ligand |

So both patterns coexist in one payload. The original note flagged exactly this
as out of scope:

> If Paul needs true small-molecule ligand modeling in AF3 (`ligand` entries with
> `ccdCodes` or `smiles`), that is a different pattern than what this repository
> currently uses for ubiquitin jobs.

It turns out they compose without conflict: ubiquitin as a protein chain,
the XisoK modification as a `userCCD` ligand, and `bondedAtomPairs` linking both.
Two bonds per job, one of each kind.

### 3. Two bonds instead of one

Every job carries:

1. `A:<site>:NZ` ↔ `L:1:C01` — the XisoK modification onto the substrate lysine.
   Note this one *does* use `NZ`, in the friend's sense: the acyl group attaches
   to a lysine ε-amine. It is present in both states.
2. either the thioester (`U:76:C` ↔ `B:91:SG`) **or** the product isopeptide
   (`U:76:C` ↔ `L:1:N01`), never both.

`build_charged_jobs.py` asserts that the two states differ in exactly this way —
that the charged state does not pre-form the product bond, and the product state
does not retain the thioester. Getting that wrong would silently answer a
different question.

## If the collaborator wants to regenerate these through their own pipeline

The four JSONs are self-contained and runnable as-is, so this is optional. But if
it is preferable to emit them through `iterate_through_ubiquitin` and
`build_alphafold_json`, the only things needed from our side are:

* the three sequences (in `seqs.json`, fetched from UniProt: Q96B02, P61956,
  P0CG47 residues 1–76)
* the `userCCD` text for `LIG-1` (in `LisoK_userCCD.cif`)
* the bond list above, with construct numbering (native K11 → residue 12)

The one thing that must survive any regeneration: **the acceptor for the
ubiquitin bond is `SG` of Cys91, not an `NZ`.** Everything else is convention.
