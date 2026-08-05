# Superseded — do not send these

The two `*_charged.json` files here express the thioester as a bond between two
PROTEIN chains (`U:76:C` → `B:91:SG`). AF3 3.0.1 **silently discards** that bond
(`structure_cleaning.py`: "Reducing number of bonds ... 1 are polymer-polymer
bonds"), exits 0, and writes a model with those atoms 26.25 Å apart.

Use `../jobs/` instead, where ubiquitin is split into a 1–74 protein chain plus
two `GLY` ligands so no bond is polymer–polymer. See
`../docs/AF3_POLYMER_BOND_LIMIT.md`.

The two `*_product.json` files here were fine — kept only for the record; the
current ones differ solely in seed count.
