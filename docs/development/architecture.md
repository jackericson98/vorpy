# Architecture

VorPy separates molecular input, network construction, geometry calculations,
analysis, and output. The morphometric extension follows the same separation:

| Layer | Responsibility | Principal modules |
|---|---|---|
| Input/system | Molecular coordinates, radii, identities, groups | `inputs/`, `system/`, `group/` |
| Network topology | Cells, surfaces, edges, vertices, adjacency | `network/build_net.py`, `network/network.py` |
| Surface construction | Pairwise boundary meshes | `network/build_surf.py`, `network/build_surfs.py` |
| Surface mathematics | Implicit functions and local/integrated curvature | `calculations/surf.py`, `calculations/curvature.py`, `calculations/surface_geometry.py` |
| Edge mathematics | Parametric edge geometry and intrinsic descriptors | `calculations/edge_geometry.py` |
| Boundary ownership | Group/external classification, outward normals, junction angles | `network/boundary_geometry.py` |
| Analysis | Atom and group aggregation | `network/analyze.py`, `group/sort.py` |
| Output | Logs, summaries, meshes, coloring | `output/`, `group/export.py`, `interface/export.py` |

Geometry calculation and thermodynamic calibration remain separate. VorPy
produces geometric and chemistry-resolved descriptors; a later calibration
layer determines how those descriptors enter a free-energy model.

See [Morphometric Geometry Implementation](morphometric_geometry.md) for the
current method-to-code map, compatibility rules, validation gates, and
production status.
