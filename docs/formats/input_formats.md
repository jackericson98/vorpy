# Input Formats

| Format | Extension |
|---|---|
| Protein Data Bank | `.pdb` |
| Crystallographic Information File | `.cif` |
| MDL MOL | `.mol` |
| Tripos MOL2 | `.mol2` |
| GROMACS structure | `.gro` |

Not every format stores the same metadata. Differences can include atom naming, residue naming, chain information, element fields, bonds, box information, and molecule identifiers.

**TODO:** Document precisely which fields VorPy reads from each format and how missing metadata is handled.
