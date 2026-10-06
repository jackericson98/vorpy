# Apollonius dual layers for AW interfaces

Use the normal two-group interface workflow and request `dual` (alias
`apollonius`) as an export product:

```powershell
python vorpy vorpy/data/2KAI.pdb -g c A and c B -g c I -i -s nt aw -e dual -e dir output/2KAI-dual
```

This defines group A as chains A+B and group B as chain I. The existing
`interface_alpha` setting, or its existing alpha-zero default, supplies the
query; the visualization does not choose a new threshold. Large/all presets
also include dual visualization for AW interfaces. Tiny/small/medium do not.
`Interface.export(dual=True)` exposes the same product to Python callers.

The export reads the cached dual, filtration, `geometry_analysis`, and physical
representation populated by `Interface.build()`. It performs no additional
network solve, alpha-birth calculation, or curvature integration. Missing or
inconsistent caches produce an explicit error rather than a replacement
selection. Power/Primitive visualization remains available through the existing
generic dual exporters; `-e dual` on a normal interface currently targets AW.

Under `<interface-output>/dual/`:

```text
apollonius_interface.pml / .py
molecule_A.pdb, molecule_B.pdb
apollonius_selected_edges.off, apollonius_selected_vertices.pdb
apollonius_rejected_edges.off, apollonius_rejected_vertices.pdb
apollonius_unresolved_edges.off, apollonius_unresolved_vertices.pdb
dual_interface_mapping.csv, dual_pair_status.csv
generator_identity.csv
contacts.json, interface_dual_summary.json
aw/
  dual/dual_full_edges.off, generators.pdb, dual_full_*.csv
  dual/dual_full_triangles.off, dual_full_tetrahedra_*.off
  alpha_<query>/alpha_edges.off, alpha_vertices.pdb
  alpha_<query>/interface_alpha_surfaces.off
  mapping/dual_voronoi_mapping.csv, alpha_interface_surface_mapping.csv
```

Empty mesh layers are omitted. Full-dual faces/tetrahedra are present only when
the existing dual contains them; they are hidden by default. Higher-dimensional
AW alpha births are not invented by this pair-based interface analysis.

Open `apollonius_interface.pml` in PyMOL (`@path/to/apollonius_interface.pml`).
Objects are `molecule_A`, `molecule_B`, `apollonius_full`,
`apollonius_selected`, and `aw_interface`. Optional objects include
`dual_vertices`, `rejected_dual`, `unresolved_dual`, `apollonius_faces`, and
`apollonius_tetrahedra`. The groups `apollonius_dual` and `interface_geometry`
toggle the layers. The full dual is thin/subdued, selected bicolor connections
are prominent, and physical surfaces are semi-transparent. Rejected and
unresolved layers start hidden and can be enabled separately.

The full dual is the complete **existing incidence dual of this solved
network**, including its unresolved support records. It is not a new global
tessellation. Its vertices are generator centers; its straight connections
are the repository's abstract Apollonius visualization. Curved AW surface
points and triangles come directly from the selected physical representation.
Every selected dual pair has the same generator indices as the interface
network and the same pair ID as `geometry_analysis`.

`dual_interface_mapping.csv` records stable atom IDs, group identities, native
birth values, surface IDs/areas, cached smooth H/K, and support/mapping status.
`generator_identity.csv` connects PyMOL molecule atom serials with network
generator IDs and full stable atom metadata. Missing curvature remains empty.
A selected alpha pair whose physical patch
cannot be mapped remains selected with `mapping_unresolved`; it is not silently
removed. `dual_pair_status.csv` also includes rejected and unresolved pairs.
An unsupported/missing birth is **unresolved**, never rejected. A persisted
certified lower-bound exclusion can be rejected while its exact birth remains
empty. Numerical certification does not assert interval-arithmetic proof.

To isolate a contact, copy its `dual_feature_id` from the mapping table and run
in the PyMOL Python console, for example:

```python
show_contact("1:10,13")
```

Use the IDs in your own table. This replaces three reusable contact objects:
`contact_dual_edge`, `contact_generators`, and `contact_aw_surface`. The helper
copies stored coordinates/triangles; it does not invoke VorPy or build geometry.
Re-enable the two main groups to return to the overview.

The cached 2KAI example can be regenerated with
`python scripts/export_cached_2kai_interface_dual.py` when the existing
`output/comparative_study` archive/checkpoint and curvature tables are present.
That script verifies generator identity and the persisted selection, then uses
the normal exporter. It does not run the older demonstration's Power solve.
