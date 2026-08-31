# Quick Start

## Launch the GUI
```bash
python vorpy
```

## Open the File Browser
```bash
python vorpy browse
```

## Analyze a Structure
```bash
python vorpy path/to/structure.pdb
```

Supported formats: PDB, CIF, MOL, MOL2, and GRO.

## Change Network Type
```bash
python vorpy example.pdb -s nt pow
```

Primary network types: `aw`, `pow`, and `prm`.

## Restrict the Analysis
VorPy can operate on atom, residue, chain, and molecule indices; general molecular classes such as protein, ligand, DNA, and RNA; and increasingly specific structural identities such as `chain A` and `THY 124`.

## Choose Exports
```bash
python vorpy example.pdb -e small
python vorpy example.pdb -e small and shell
```
