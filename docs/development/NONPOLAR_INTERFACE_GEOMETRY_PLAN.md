# Nonpolar Interface Geometry Plan

## Scope and scientific question

The first question is deliberately structural:

> Do local geometric descriptors distinguish nonpolar surface exposed at an existing molecular interface from nonpolar surface exposed away from that interface?

The initial target is `interface surface` versus `non-interface surface`. “Active area,” function, binding site, and biological activity are out of scope for the first analysis. Whole-molecule integrated Gaussian curvature is also out of scope as a local predictor; it is topology-constrained and should be reported, if at all, only as a QC quantity.

The analysis should consume a solved VorPy `System`/`Network` where possible. CSV logs are useful for reproducibility and batch analysis, but they omit mesh-level ownership and some topology needed for an authoritative interface label.

## Existing architecture and data audit

### System, atoms, and identity

`vorpy/src/system/system.py` owns the loaded molecular system. `System.balls` is the parent atom/cell table. Input parsing in `vorpy/src/inputs/pdb.py` creates atom records through `vorpy/src/objects/atom.py`; records retain:

- stable local atom/cell index (`num` and the DataFrame index)
- coordinates (`loc`) and Voronoi radius (`rad`)
- element, atom name, residue name/sequence, chain name, segment, and mass
- links to residue/chain objects and network topology lists (`surfs`, `edges`, `verts`)

The PDB reader also distinguishes water and common ion categories through `chemistry_interpreter.py`. Virtual water charge sites are excluded from the atom geometry representation. The identity fields are sufficient to produce `system_id`, `frame_id`, chain, residue, atom name, element, and atom index, subject to an explicit index convention.

### Network and topology

`vorpy/src/network/network.py` stores the solved representation:

- `net.balls`: atom/cell rows, including `loc`, `rad`, `surfs`, `edges`, `verts`, and analyzed properties.
- `net.surfs`: pairwise Voronoi surface/contact faces. Important columns include `balls`, `points`, `tris`, `mean_tri_curvs`, `gauss_tri_curvs`, `mean_curv`, `avg_mean_curv`, `gauss_curv`, `avg_gauss_curv`, `int_mean_curv`, `int_mean_curv_sq`, `int_gauss_curv`, `sa`, `vols`, `contact_area`, `overlap`, and surface-energy-related values.
- `net.edges`: three-generator edge records, length, endpoints, and (for current AW calculations) cell-relative mean/Gaussian curvature contributions.
- `net.verts`: four-generator vertices, location/radius, and cell-relative vertex Gaussian curvature contributions.

`vorpy/src/network/analyze.py` performs the cell aggregation. For AW networks it explicitly combines:

`H_cell = H_surface + H_edge`

`K_cell = K_surface + K_edge + K_vertex`

It also calculates surface area, volume, maxima and averages of local triangle curvature, `integrated mean curvature squared`, neighbor lists/counts, contact area, overlap quantities, completeness, center of mass, moments, and bounding boxes. The AW cell-relative caches (`int_mean_curv_by_ball`, `int_gauss_curv_by_ball`, and related surface/edge/vertex caches) are important because signed mean curvature depends on the cell side.

### Surface-energy calculations

`vorpy/src/calculations/surface_energy.py` integrates triangle-level `H` and `K` over a selected surface mesh. It returns area-weighted average mean/Gaussian curvature and:

- `Area`
- `Integrated Mean Curvature`
- `Integrated Mean Curvature Squared`
- `Integrated Gaussian Curvature`

Its morphometric and Helfrich helpers convert those geometric integrals into energy-like quantities. These are derived geometry/parameterized models, not independently measured biological energies. The proof of concept should retain the geometry fields and treat `Representative Surface Energy` as an optional secondary descriptor, with its formula and coefficients recorded.

### Groups and interfaces

`vorpy/src/group/group.py` defines selected groups as subsets of system ball indices and builds a `Network` for each group. `vorpy/src/interface/interface.py` creates an `Interface` between two groups (or a group and a surrounding set). It stores `group1_indices`, `group2_indices`, `interface_id`, partial groups, and an interface-owned network. The interface network receives:

- the full system coordinates/radii for geometric validity
- a restricted `group` containing both sides
- `iface_grps=(group1_indices, group2_indices)`

`vorpy/src/network/build_net.py` contains `spans_interface`, which identifies a topology element whose three generators include members from both interface sides. This is the strongest existing structural signal for an interface edge/vertex. Interface surfaces are pairwise faces whose two generators lie on opposite sides; this must be implemented as an explicit topology predicate and validated against the interface network rather than inferred from an arbitrary distance cutoff.

The interface’s surrounding-group construction currently uses a spatial candidate radius to assemble the network context. That radius is a construction envelope, not the final interface label. The label should be based on the solved cross-side surface topology.

### Contacts and ownership

`network/analyze.py` derives first-neighbor relationships from `net.surfs['balls']`, assigns each cell its surface/edge/vertex lists, and calculates `surface_contact_area` and `contact_area`. Thus a surface has explicit two-ball ownership, and an atom/cell can inherit area and curvature by summing its owned surfaces. This supports both observation levels:

1. surface rows, with exact pair ownership and interface status;
2. atom/cell rows, with sums or area-weighted aggregates over owned surfaces.

The surface row is the better primary label unit. Atom rows should be a first convenience table, not treated as independent physical samples.

## What is already exported to `aw_logs.csv`

The current writer is `vorpy/src/output/logs.py`; the filename is commonly `<network_type>_logs.csv` in the active output directory.

### Atom section

The `Atoms` section currently exports, where the cell is complete and has nonzero area:

- identity: index, name, residue, residue sequence, chain, mass, coordinates, radius;
- volume, VDW volume, surface area, completeness;
- maximum and average mean curvature;
- maximum and average Gaussian curvature;
- integrated mean curvature, integrated mean curvature squared, integrated Gaussian curvature;
- representative surface-energy-like value;
- sphericity/isoperimetric quotient, inner-ball flag;
- neighbor count/nearest neighbor and distance summaries;
- spike distances, overlap/contact/non-overlap quantities, COM, inertia, bounding box, neighbors;
- explicit face/edge/vertex/total curvature breakdown columns.

This is enough for a useful atom-level geometry table, but the log does not encode the underlying surface IDs or triangle ownership. The atom `Surface Area` is total exposed network-cell area, not necessarily “nonpolar area” and not necessarily a uniquely labeled interface area.

### Surface section

The `Surfaces` section exports pair ownership (`Ball 1`, `Ball 2`), area, mean/Gaussian curvature and averages, all three principal curvature integrals, representative surface energy, volume contributions, contact area, overlap, and a Gaussian face field. This is the most useful current log section for an exploratory surface-level analysis.

### Edge and vertex sections

Edges export their three generator balls, length, three cell-relative integrated mean terms, three cell-relative integrated Gaussian terms, and six directed face-incidence Gaussian terms. Vertices export four generator balls, location/radius, and integrated Gaussian curvature. These are useful for diagnostics and cell aggregation, but edge/vertex rows are not surface observations and should not be mixed into the first classifier without a clear area measure.

### Build/group metadata

Build metadata includes input location, network type, resolutions, timing, and VorPy version. Group metadata includes volume, surface area, mass, density, COM/inertia, integrated curvature totals, Euler/Gauss-Bonnet fields, and boundary completeness/closure/manifold/orientability diagnostics.

`vorpy/src/analyze/tools/compare/read_logs2.py` now reads the modern header-driven fields, including curvature breakdowns, per-face/per-ball edge values, vertex values, and boundary diagnostics. It remains a good interchange reader, but it cannot recover mesh ownership or the original `Interface` object from a log alone.

## Missing or not yet reliable data

The following should not be silently fabricated in the first table:

- a canonical `system_id` and `frame_id` are not currently written in every log;
- a canonical surface-side/interface label is not exported as a column;
- an atom’s “interface area” and “interface fraction” are not currently stored as explicit fields;
- `nonpolar` is not a first-class VorPy atom property;
- neighborhood-averaged curvature is not stored and requires a documented aggregation over the surface/atom adjacency graph;
- triangle-level geometry is in-memory in `net.surfs`, but the standard CSV does not export triangle IDs/coordinates;
- surface `contact_area` describes VorPy contact/overlap geometry and must not automatically be equated with cross-group interface area;
- a surface may be owned by two generators, while an interface classification is a relation between two group sides. The relation must be evaluated from `iface_grps` and the surface’s `balls`.

## Proposed observation unit

Use a two-table design, generated by a separate analysis module after a network is solved.

### Primary: surface observation

One row per valid `net.surfs` surface:

```text
system_id, frame_id, network_id, surface_id
ball_1, ball_2
chain_1, residue_1, residue_seq_1, atom_name_1, element_1
chain_2, residue_2, residue_seq_2, atom_name_2, element_2
surface_area
mean_curvature, average_mean_curvature
gaussian_curvature, average_gaussian_curvature
integrated_mean_curvature, integrated_mean_curvature_squared
integrated_gaussian_curvature
representative_surface_energy
volume_contribution_1, volume_contribution_2
contact_area, overlap
side_1, side_2, spans_interface, interface_label
nonpolar_1, nonpolar_2, nonpolar_surface_label
```

The first positive class should be a surface whose two generators are in opposite interface groups. The first non-interface comparison set should be exposed surfaces within the same group, with the same network and completeness criteria. Surfaces involving water/ions or unsupported chemical categories should be excluded or reported in a separate stratum, not silently pooled.

### Secondary: atom/cell observation

One row per exposed, complete `net.balls` cell, with an area-weighted sum over owned surfaces and the already-calculated cell fields:

```text
system_id, frame_id, network_id, atom_index
chain, residue_name, residue_number, atom_name, element, radius
volume, surface_area, exposed_surface_area
max_mean_curvature, average_mean_curvature
max_gaussian_curvature, average_gaussian_curvature
integrated_mean_curvature, integrated_mean_curvature_squared
integrated_gaussian_curvature
integrated_mean_curvature_face, integrated_mean_curvature_edge
integrated_gaussian_curvature_face, integrated_gaussian_curvature_edge
integrated_gaussian_curvature_vertex
neighbor_count, neighbors
interface_area, interface_fraction, interface_label
nonpolar_class
```

An atom is `interface_label=1` only if its interface-associated surface area exceeds a declared positive threshold (initially `> 0`, with sensitivity analyses for a minimum area). The atom table must carry `surface_count` and `interface_surface_count`; otherwise atom-level labels conceal how many faces support the label.

Patch rows should wait until these two tables and their topology tests are stable. A future patch can be a connected component of nonpolar surface triangles or surface IDs, with area-weighted curvature summaries and a clear graph definition.

## Nonpolar definition

VorPy currently has element, residue, atom-name, and broad residue-category information, but no authoritative `polar_class` field. Element-only classification is inadequate for scientific claims: carbonyl carbon, aromatic carbon, heteroatom-containing residues, protonation variants, and ligand atom types need chemical context. The existing exploratory script `F4_atomic_curvature/SurfaceEnergyCurvature.py` uses `C/H` versus `O/N/S/P`, but this is a transparent pair-energy heuristic, not a validated atom-level hydrophobicity model.

For the first proof of concept, define and publish two labels:

1. **Strict element baseline:** nonpolar atom if element is `C` or `H`; exclude water, ions, unknown elements, and virtual sites.
2. **Residue-aware sensitivity label:** use an explicit, versioned table for atom/residue types; classify ambiguous heteroatom-containing residues separately rather than forcing them into nonpolar.

The strict element label is acceptable as a baseline only. It should not be described as a physical solvent-exposed hydrophobicity classification. A later force-field or solvent-accessible-atom typing layer can replace it, while retaining the original label and rule version in the output. For a surface pair, initially define `nonpolar_surface_label` as both generators passing the strict rule; also report mixed pairs as a separate category.

Exclude waters, ions, virtual sites, and incomplete cells from the primary comparison. Hydrogens should be retained only if the input/typing policy treats them consistently; otherwise report heavy-atom-only results as the main analysis and hydrogen-inclusive results as sensitivity analysis.

## Interface definition

For a two-group interface network:

```python
def spans_interface(surface_balls, iface_grps):
    side_a, side_b = map(set, iface_grps)
    balls = set(surface_balls)
    return bool(balls & side_a) and bool(balls & side_b)
```

For a two-generator surface, the stronger label is exactly `ball_1 in side_a and ball_2 in side_b` or the reverse. The existing `iface_grps`/`spans_interface` machinery already provides the side membership; no new distance cutoff is needed. For a group-versus-surrounding interface, use the same side sets after the surrounding set has been materialized. Keep the construction-envelope distance (`get_surrounding_indices`) out of the final label.

At atom level, sum surface area over cross-side surfaces for each generator/cell and define:

`interface_fraction = interface_area / exposed_surface_area`.

Store the raw area and fraction, not just a binary label. This permits threshold and class-balance sensitivity checks. For protein-ligand analyses, group membership must be declared from the input selection and recorded in metadata; “protein” and “ligand” should not be inferred solely from an arbitrary residue-name fallback.

## First feature set and exploratory analysis

Start with nonredundant, local or cell-level features:

- exposed surface area and volume;
- mean curvature and Gaussian curvature averages;
- maximum mean/Gaussian curvature;
- integrated mean curvature;
- integrated mean-curvature-squared;
- integrated Gaussian curvature only as a local cell/surface term;
- neighbor count and neighborhood mean/Gaussian curvature after adjacency is defined;
- interface area/fraction only as the label-derived structural quantity, never as a predictor in a leakage-prone model.

Use signed and absolute mean curvature deliberately. AW signed mean curvature is cell-relative; report the sign convention and include `|H|` or squared terms where orientation-invariant comparisons are desired. Preserve units: area Å², volume Å³, `H` Å⁻¹, `K` Å⁻², integrated `H` Å, integrated `H²` dimensionless under Å units, and integrated `K` dimensionless.

For each descriptor, report sample count, missingness, median/IQR, robust effect size (rank-biserial or Cliff’s delta), bootstrap confidence interval at the system level, and both Pearson and Spearman correlations where correlation is relevant. Plot distributions and interface/non-interface scatter plots. PCA is descriptive only and must be fit within the training partition if later used for prediction.

The sampling unit for uncertainty is the molecular system/complex, not the atom. Atom rows from one structure are correlated. Report per-system summaries alongside pooled rows and use clustered or hierarchical statistics for the first comparison.

## Simple classification after exploration

Use interpretable regularized logistic regression only after labels/features pass QC:

- **Model A:** chemistry/nonpolar baseline (element, residue/atom category, radius, area/volume controls);
- **Model B:** geometry alone;
- **Model C:** chemistry plus geometry.

Scale numeric features within each training fold. Remove or group near-duplicates (for example, area-weighted sums that algebraically define another field). Compare held-out system-level AUROC, AUPRC, balanced accuracy, calibration, and confidence intervals. The key result is whether Model C improves over Model A on unseen complexes, not whether geometry-only separates atoms from the same protein.

Hold out complete proteins/complexes, preferably by sequence/structural family where possible. Never perform an atom-level random split across one molecule. Frames from one trajectory belong to the same system group and must remain in the same split; frame stability is a later longitudinal analysis.

## Output format

Create a separate analysis export, for example:

```text
nonpolar_interface_geometry.parquet
nonpolar_interface_geometry.csv
```

Parquet is preferred for nested neighbor metadata and numeric types; CSV is a readable interchange option. Include a sidecar JSON/YAML metadata record containing VorPy version, input path/hash, system/frame IDs, network type/settings, index convention, nonpolar rule/version, excluded categories, interface group IDs, area threshold, completeness rule, units, and analysis timestamp. Do not overload `aw_logs.csv` with classifier-specific labels.

## Implementation phases

### Phase 0 — audit and fixtures

This document is the Phase 0 result. Add small solved-network fixtures with known two-group topology, one nonpolar/nonpolar surface, one mixed surface, an incomplete cell, and a water/ion exclusion case.

### Phase 1 — read-only table builder

Add a module under `vorpy/src/analyze/nonpolar_interface/` that accepts a solved `Network` plus system/interface metadata and returns surface and atom DataFrames. It should not modify `System`, `Network`, or solver tables. Implement identity joins, area/curvature extraction, side membership, nonpolar rule versioning, completeness filtering, and QC counts.

### Phase 2 — exploratory reports

Add scripts/notebooks that consume the analysis table and produce descriptive summaries, distributions, effect sizes, correlation matrices, and PCA. Keep plots and statistical code outside the solver package where practical.

### Phase 3 — held-out baseline models

Implement grouped cross-validation and Models A–C in a separate analysis package/module. Store predictions and fold assignments with system IDs.

### Phase 4 — frames and patches

Add `frame_id` from the trajectory runner, repeated-measures summaries, and stability metrics. Only then define connected patch construction over a documented surface adjacency graph.

## Specific files likely to change later

Initial analysis work should add new files rather than alter the solver:

- `vorpy/src/analyze/nonpolar_interface/table.py` — surface/atom observation tables;
- `vorpy/src/analyze/nonpolar_interface/classification.py` — explicit nonpolar rules and interface predicates;
- `vorpy/src/analyze/nonpolar_interface/export.py` — Parquet/CSV plus metadata;
- `vorpy/src/analyze/nonpolar_interface/explore.py` — summaries and plots;
- `vorpy/tests/analyze/test_nonpolar_interface_table.py` — topology/label fixtures.

Only if profiling shows a real need should `vorpy/src/output/logs.py` gain optional `system_id`, `frame_id`, or explicit interface columns. `read_logs2.py` should remain backward-compatible and need not become the canonical source for topology labels.

## Tests and QC requirements

- surface ownership and area sums agree with `net.balls` aggregation;
- cross-side surface predicate is symmetric under swapping group sides;
- surfaces with same-side generators are non-interface;
- atom interface area equals the sum of cross-side owned surface areas;
- incomplete cells are excluded by default and counted in QC output;
- water, ion, virtual-site, and unknown categories are excluded or explicitly stratified;
- strict `C/H` rule is deterministic and versioned;
- signed AW curvature is preserved, while orientation-invariant derived features are reproducible;
- no whole-molecule Gaussian integral enters the local feature table;
- repeated frame IDs and atom IDs remain stable and split-safe;
- logs and direct-network table builders agree on shared scalar fields within rounding tolerance.

## Smallest scientifically useful first implementation

Implement one read-only surface/atom table builder for an already-solved two-group interface network. Use strict `C/H` heavy-atom nonpolar baseline, exclude waters/ions/incomplete cells, label cross-side `net.surfs` by `iface_grps`, aggregate `interface_area` and `interface_fraction` to atoms, and export a table plus metadata. Run descriptive system-clustered comparisons of area, local `H`, `H²`, `K`, and curvature integrals before adding patches, machine learning, new solver fields, or changes to the standard logs.

## Phase 1 status

The first read-only table builder now lives in `vorpy/src/analyze/nonpolar_interface/`. `build_geometry_tables(network, ...)` returns `surfaces`, `atoms`, and `metadata` DataFrames/dictionary without mutating the solved network. It uses the versioned strict `C/H` rule, excludes recognized water/ion residues from the nonpolar class, labels cross-side surfaces from `iface_grps`, and aggregates cross-side area and fraction to atom rows. `export_geometry_tables(...)` writes separate surface/atom CSVs and a metadata JSON (with opt-in Parquet when a pandas engine is installed). `summarize_surface_features(...)` provides descriptive count/mean/median/IQR summaries grouped by chemistry and interface label. Large/all exports now include these tables under each group or interface directory in `nonpolar_interface_geometry/`; ordinary group networks retain missing interface labels, while interface networks use their explicit cross-side groups. Tests cover same-side versus cross-side labels, missing interface metadata, exclusions, invalid rule names, export round-trips, preset integration, and descriptive summaries. Plotting and real multi-complex aggregation remain the next phase so the table schema can be reviewed against real solved interface networks first.
