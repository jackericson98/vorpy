# Morphometric Geometry Implementation

This is the living connection between VorPy's scientific methodology and its
implementation. Update this page in the same commit whenever a formula, field,
algorithm, validation rule, or output contract changes.

Status labels:

- **Production**: used in normal VorPy calculations and outputs.
- **Validated scaffold**: implemented and unit tested but deliberately not
  connected to production descriptors.
- **Planned**: design exists, but implementation or scientific validation is
  incomplete.

## Method-to-code map

| Method or quantity | Definition | Implementation | Tests | Status |
|---|---|---|---|---|
| AW pairwise surface | `F_ij(x)=0` | `calculations/surf.py::calc_surf_func` | `tests/calculations/test_surf.py` | Production |
| Local mean curvature | Analytic implicit-surface `H` | `calculations/curvature.py::mean_curvature` | `tests/calculations/test_surface_geometry.py` | Production formula; direct analytic coverage |
| Local Gaussian curvature | Analytic implicit-surface `K` | `calculations/curvature.py::gaussian_curvature` | `tests/calculations/test_surface_geometry.py` | Production formula; direct analytic coverage |
| Surface integration | Triangle-centroid area quadrature | `calculations/curvature.py::calc_surf_tri_curvs_both` | Molecular regressions; convergence tests still needed | Production |
| Common surface API | Value, gradient, normal, `H`, `K` | `calculations/surface_geometry.py` | `tests/calculations/test_surface_geometry.py` | Validated scaffold |
| Generic edge API | Point, derivatives, speed, curvature, quadrature | `calculations/edge_geometry.py::EdgeGeometry` | `tests/calculations/test_edge_geometry.py` | Validated scaffold |
| Straight edge | Exact line segment | `calculations/edge_geometry.py::LineEdgeGeometry` | `tests/calculations/test_edge_geometry.py` | Validated scaffold |
| AW curved edge branch | `x(rho)=x0+u rho +/- n sqrt(q)` | `calculations/edge_geometry.py::AdditivelyWeightedTrisectorBranch` | `tests/calculations/test_edge_geometry.py` | Validated scaffold; not connected to network edges |
| Complete AW conic chart | Angle, hyperbolic-angle, or normal-coordinate parameterization across turning points | `calculations/edge_geometry.py::AdditivelyWeightedTrisectorConic` | `tests/calculations/test_edge_geometry.py` | Validated diagnostic fallback |
| AW endpoint matcher | Recover `rho`, select branch, validate endpoints | `calculations/edge_geometry.py::match_aw_trisector_to_endpoints` | `tests/calculations/test_edge_geometry.py` | Validated diagnostic; production network not changed |
| AW sample validator | Maximum and RMS reconstruction residual | `calculations/edge_geometry.py::validate_edge_samples` | `tests/calculations/test_edge_geometry.py` | Validated diagnostic |
| AW network-edge audit | Per-edge analytic match, sampled residual, and length comparison | `network/edge_geometry_diagnostics.py::diagnose_aw_edge_geometry` and `Network.diagnose_aw_edges` | `tests/network/test_build_net.py` | Validated diagnostic; read-only and opt-in |
| AW audit CLI | `--diagnose-edges` after build and before export | `command/vpy_cmnd.py::Command.run_edge_diagnostics` | `tests/command/test_vpy_cmnd.py` | User-visible diagnostic; no production-field changes |
| Edge metadata adapter | Stable view of `net.edges["vals"]` | `calculations/edge_geometry.py::normalize_edge_construction` | `tests/calculations/test_edge_geometry.py` | Validated scaffold |
| Canonical edge tangent | Endpoints ordered by vertex index | `calculations/edge_geometry.py::canonical_edge_endpoints` | `tests/calculations/test_edge_geometry.py` | Validated scaffold |
| Barrier-edge ownership | One/two group generators versus external generators | `network/boundary_geometry.py::get_boundary_surfaces` | `tests/network/test_boundary_geometry.py` | Validated scaffold |
| Outward surface normal | Group generator toward external generator | `network/boundary_geometry.py::outward_surface_normal` | `tests/network/test_boundary_geometry.py` | Validated scaffold |
| Signed junction angle | `atan2(T dot (n1 cross n2), n1 dot n2)` | `network/boundary_geometry.py::signed_surface_intersection_angle` | `tests/network/test_boundary_geometry.py` | Validated scaffold |
| Edge mean-curvature term | `1/2 integral theta ds` | Not connected | Analytic wedge and molecular validation required | Planned |
| Intrinsic edge curvature | `integral kappa ds` and `integral kappa^2 ds` | Generic quadrature in `EdgeGeometry` | Circle, line, and rigid-transform tests | Validated scaffold; optional descriptor only |
| Gaussian edge term | Gauss--Bonnet-derived contribution | Not implemented | Analytic reference geometries required | Planned |
| Gaussian vertex term | Corner/normal-space contribution | Not implemented | Solid-angle and closed-boundary tests required | Planned |
| Representative energy | `2 integral H^2 dA` | `calculations/surface_energy.py` and `network/build_surfs.py` | Existing log/mesh regressions | Production; surface-only reference quantity |
| Morphometric free energy | `pV + gamma A + kappa M + kappa_bar G` | Calibration layer not implemented | FreeSolv train/validation protocol required | Planned |

## Current data flow

```text
calc_surf_func
  -> build_surf
  -> calc_surf_tri_curvs_both
  -> net.surfs curvature fields
  -> network.analyze atom sums
  -> group/interface boundary sums
  -> logs, summaries, coloring, and exports
```

The canonical `int_mean_curv` and `int_gauss_curv` fields currently follow
this surface-only path. None of the edge scaffolding changes current output.

The optional network-edge audit follows a separate read-only path:

```text
net.edges + net.verts + net.balls
  -> Network.diagnose_aw_edges()
  -> analytic branch or straight-line matching
  -> per-edge residual and length records
  -> aggregate diagnostic summary
```

It does not modify `net.edges` or feed analysis, logs, exports, or the
Workbench.

The matcher attempts the already validated clearance branch first. Only
`cross_branch` and `turning_point` cases proceed to the sampled-path-selected
nonsingular conic. Reports distinguish these recoveries as
`matched_nonsingular` so coverage changes remain auditable.

The CLI exposes this path with `--diagnose-edges`. The report is printed after
network construction and before export. It remains observational: invoking the
flag does not change the network or any stored scientific quantity.

## Compatibility rules

1. Preserve `int_mean_curv`, `int_mean_curv_sq`, `int_gauss_curv`, and
   `surf_energy` until a documented versioned transition.
2. Develop complete-boundary quantities under explicit component names.
3. Keep `int_mean_curv_sq` and representative `surf_energy` separate from the
   classical edge correction to integrated mean curvature.
4. Do not introduce a fitted coefficient for intrinsic edge curvature as a
   classical morphometric term without empirical model comparison.
5. Do not implement Gaussian edge or vertex terms heuristically.

## Validation gates before production connection

Analytic AW edges must:

1. Reproduce both stored endpoint vertices.
2. Agree with existing sampled `edge["points"]` within a documented tolerance.
3. Handle generator and endpoint reordering deterministically.
4. Handle conic turning points by explicit splitting or nonsingular
   parameterization.
5. Preserve length and curvature under rigid translations and rotations.

Complete integrated mean curvature must then pass:

1. Plane and straight-edge limits.
2. Constant-angle wedge tests.
3. Convex and concave sign tests.
4. Closed analytic geometries with known results.
5. Resolution-convergence tests.
6. Molecular regression tests confirming unchanged legacy fields.

## Documentation rule for future stages

Every implementation stage must update, together:

- the relevant theory page;
- this method-to-code map;
- docstrings and units for new public fields;
- unit or regression tests;
- output-format documentation when a field becomes user-visible;
- release notes when an existing public meaning changes.
