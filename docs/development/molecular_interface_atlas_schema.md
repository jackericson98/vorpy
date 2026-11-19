# Molecular Interface Atlas — schema and capability boundary

Status: Phase 2 cache-only exporter demonstrated. This document does not
authorize a molecular solve, frame batch, wet calculation, or geometry rebuild.

## 1. Common record envelope

Every machine-readable atlas row uses the `molecular_interface_atlas.v1`
envelope. IDs are derived from producer/archive identity and canonical feature
identity; they must not depend on an output path or row order.

Required fields:

```json
{
  "record_id": "atlas:v1:<source-network>:<method>:<representation>:<kind>:<feature-key>",
  "kind": "surface | edge | vertex",
  "structure_id": "2KAI | 1AAY",
  "state": "dry | solvated",
  "frame_id": "static | fN",
  "method": "prm | pow | aw | cazals_alpha",
  "representation": "primal | dual | alpha_selected | boundary | interface",
  "orientation": {"from": "A", "to": "B"},
  "status": "CERTIFIED | PARTIAL | UNRESOLVED | NOT_CALCULATED | NOT_SUPPORTED | STALE",
  "status_reason": "explicit reason when the value is unavailable or partial",
  "source_identity": {
    "input_role": "_og | solvated_post_md | dry_post_md | frame",
    "input_sha256": "...",
    "archive_id": "...",
    "source_network_id": "...",
    "solve_universe_id": "...",
    "generator_digest": "..."
  },
  "provenance": {"producer": "...", "settings_digest": "...", "units": {}}
}
```

`status` is scientific state, not export success. A missing or unsupported
quantity is `null` with a non-available status; zero is never used as a
placeholder. Ordered interfaces retain both group IDs and the orientation
(`A -> B` versus `B -> A`).

## 2. Primitive schemas

### Surface rows

Required: `surface_id`, `generator_ids`, `source_atom_ids`,
`source_residue_ids`, `source_chain_ids`, `surface_type`, `selection_status`,
`bounded`, `complete`, `incident_edge_ids`, `incident_vertex_ids`,
`geometry` and `area_A2`.

`surface_type` distinguishes `intra_group`, `cross_group`, `support`,
`boundary`, `alpha_selected`, `alpha_excluded`, and `unclassified`.
`geometry` contains the source mesh or an explicit `NOT_SUPPORTED` status;
generator-center dual triangles are marked abstract and are not substituted for
physical surfaces. Mean/Gaussian curvature fields carry their native units and
individual statuses.

### Edge rows

Required: `edge_id`, `endpoint_vertex_ids`, `adjacent_surface_ids`,
`generator_ids`, `analytic_type`, `bounded`, `complete`, `length_A`,
`angle_rad`, `mean_curvature_contribution_A_rad`,
`gaussian_curvature_contribution`, `incidence_status` and provenance.

An edge with one or zero resolved incident faces is retained as
`boundary`/`unbounded`/`UNRESOLVED`; it is not dropped from an aggregate and its
unknown contribution is not replaced by zero.

### Vertex rows

Required: `vertex_id`, `xyz_A`, `generating_atom_ids`,
`incident_edge_ids`, `incident_surface_ids`, `bounded`, `complete`,
`angle_defect_rad`, `gaussian_curvature_contribution`, `status` and
provenance. Undefined angle defects remain null with an explicit status.

## 3. Identity and reconciliation rules

- Primal feature keys use the producer network namespace plus the stable
  feature ID. Dual keys use sorted stable generator IDs and dimension.
- Atom, residue and chain IDs come from the audited input identity, not from
  display names or a reordered export.
- Aggregate area is reconciled to the unique primitive surface set. The report
  includes supported/total counts, unresolved IDs, duplicate IDs and the
  absolute/relative residual.
- Edge and vertex totals are reconciled only when incidence is complete and
  the population is explicitly defined. Boundary and alpha-selected
  populations are separate from full-network populations.
- Gauss–Bonnet certification requires complete incidence, stated boundary
  terms, orientation and a residual. Otherwise the total is `PARTIAL` or
  `UNRESOLVED`.
- PRM, ordinary POW, AW and Cazals-style alpha selection have separate
  `method`/`representation` identities and units. No method equivalence is
  inferred.

## 4. Export and visualization contract

The implemented cache-only bundle is:

```text
atlas_manifest.json
surfaces.csv / surfaces.jsonl
edges.csv / edges.jsonl
vertices.csv / vertices.jsonl
interface_summary.csv
validation_report.json
geometry/surfaces_<classification>.off
geometry/edges.obj
geometry/vertices.pdb
geometry/geometry_index.json
visualization.pml
```

The PyMOL launcher must load surfaces, edges and vertices as independent
objects, preserve group/interface colors, and expose full, boundary and alpha
selection objects separately. For AW, straight generator-center dual meshes
are labeled visualization-only; physical AW geometry remains the solved
network geometry.

## 5. Current capability boundary

Already present in VorPy:

- Network archives retain named balls, vertices, edges, surfaces, atom
  metadata, generator provenance and solvent membership.
- `InterfaceSurface`, `InterfaceEdge` and `InterfaceVertex` provide a useful
  cached interface record model and CSV export, including selection,
  boundedness, completeness and incidence status.
- Typed Results provide identity, quantities, topology, dependencies,
  `GeometryLineage`, unresolved statuses and archive round-trips when a typed
  interface result is attached.
- Dual/alpha visualization exports write full and selected layers, mapping
  tables, OFF/PDB geometry and PyMOL launchers without rebuilding geometry.
- Existing synthetic tests cover archive/cache-only preservation and several
  bounded/alpha/incidence failure modes.

Implemented in Phase 2 for the archived 1AAY ordinary POW network:

- `vorpy/src/atlas/cache_export.py` reads the archive without rebuilding
  geometry and refuses an existing output directory;
- `surfaces.csv`, `edges.csv`, `vertices.csv`, `contacts.csv`,
  `interface_summary.csv`, `validation_report.json`, and `atlas_manifest.json`
  are demonstrated in `output/1AAY_atlas_cache_export_v6/`;
- all 14,983 surfaces, 27,191 edges and 13,839 vertices are retained, with
  explicit partial/nonmanifold statuses and archive/source hashes;
- PyMOL layers are independently loadable, and `geometry_index.json` maps
  table IDs to surface meshes, edge lines, vertex pseudoatoms and a selected
  surface highlight layer.

Still missing for the full Atlas contract:

- complete group/solvent/ion classification and ordered interface catalogs for
  every requested 2KAI and 1AAY pair;
- primitive-to-aggregate reconciliation and method-labelled diagnostics in a
  single manifest;
- validated PRM/AW/Cazals-alpha comparison datasets;
- a frame-directory identity/checksum resolver and dry/wet paired catalog;
- proof that typed `GeometryLineage` survives every required archive path.

The next implementation seam is typed Results/GeometryLineage attachment for
legacy raw archives and separately approved method/wet-state datasets. The
Phase 2 pairwise contact rows remain adjacency classifications, not ordered
interface solutions.
