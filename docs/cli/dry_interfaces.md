# Interfaces with optional solvent

Run the normal interface path with two molecular groups:

```console
python vorpy 2kai -g c 0 and c 1 -g c 2 -i
```

Solvent is optional. Water and ions use the shared residue classification in
`vorpy/src/boundary.py`; HETATM alone does not imply solvent. PDB author chain
labels do not determine solvent ownership. Virtual boundary generators are
numerical geometry and do not become physical solvent residues.

The bundled `2KAI.pdb` contains 2,236 protein atoms and 10 crystallographic water
oxygens. The command above preserves those input records. A water-free copy of
the same structure exercises the genuinely dry path with the same chain groups.

The default Large interface export includes the full physical interface in
`surfs.off` and a separate `alpha_interface` directory. Open
`alpha_interface/selected_interface.pml` in PyMOL to inspect the alpha-selected
physical interface. Smaller info-only presets describe the same analysis;
`-e surfs` writes both full and selected physical geometry explicitly.

The selected directory includes `selected_surfaces.off`, physical edges and
vertices, `mapping.csv`, `selection.json`, and the PyMOL loader. Selection is
read from the cached `InterfaceGeometryAnalysis` representation. The exporter
checks analysis pair IDs against the cached mappings and checks each physical
surface's generator pair against its solved network row. Mapping failures stay
explicit in the CSV; missing geometry is never manufactured. Export performs
zero network solves and does not independently query the alpha filtration.

AW keeps the existing experimental additive restricted-cell filtration. Power
and Primitive keep their respective weighted and ordinary alpha filtrations.
Curvature totals retain the existing completeness/certification requirements;
partial accounting remains partial whether solvent is present or absent.

When Apollonius visualization is enabled by the existing export preset or
`-e dual`, `dual/apollonius_interface.pml` displays the cached full/restricted
dual and selected physical AW surfaces from the same solved network. This
selected-interface export does not add dual exports to ordinary group solves.

Regression coverage is in `vorpy/tests/interface/test_dry_pdb_interface.py`.
Its slow test runs the real command pipeline against a water-free 2KAI copy,
checks one network solve, and compares exported physical mesh coordinates and
surface IDs with the cached selected geometry.
