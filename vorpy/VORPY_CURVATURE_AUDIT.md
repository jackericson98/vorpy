# VorPy curvature integration audit

## Outcome

The current canonical fields are integrations over triangulated pairwise
surface patches only. No edge or vertex contribution is included in
`int_mean_curv` or `int_gauss_curv`.

This package adds a common implicit-surface API without redirecting any current
calculation through it. Existing numerical outputs therefore remain unchanged.

## Current data flow

1. `src/calculations/surf.py::calc_surf_func()` creates the 14-value implicit
   quadratic coefficient vector for an AW pairwise surface.
2. `src/network/build_surfs.py::build_surfs()` orders the two generators by
   radius and calls `src/network/build_surf.py::build_surf()`.
3. `build_surf()` triangulates the pairwise patch. For a curved AW patch it
   calls `calc_surf_tri_curvs_both()`; Power, Primitive, and equal-radius AW
   patches are marked flat and assigned zero curvature.
4. `src/calculations/curvature.py::calc_surf_tri_curvs_both()` evaluates mean
   curvature `H` and Gaussian curvature `K` at each triangle centroid and
   performs the following mesh quadratures:

   - `int_mean_curv = sum(H_t * A_t)` in angstroms
   - `int_mean_curv_sq = sum(H_t**2 * A_t)` dimensionless
   - `int_gauss_curv = sum(K_t * A_t)` dimensionless

5. `build_surfs()` stores triangle values and patch integrals in `net.surfs`.
   It sets `surf_energy = 2 * int_mean_curv_sq`.
6. `src/network/analyze.py::analyze()` sums every incident surface integral for
   each atom and stores the same canonical field names in `net.balls`.
7. `src/group/sort.py::get_info()` sums only first-layer boundary surfaces for
   a group. The raw canonical fields remain in the native implicit-function
   orientation. A separate `oriented_int_mean_curv` is produced by comparing
   `grad(F)` with the group-inside to group-outside direction.
8. The canonical values are consumed by group/interface summaries, water
   interface summaries, CSV logs, mesh coloring, CLI/GUI scheme selection, and
   log readers.

## Definitions and sign conventions

- Mean curvature uses
  `(trace(Hess(F))*|grad(F)|^2 - grad(F)^T Hess(F) grad(F)) /
  (2*|grad(F)|^3)`.
- Gaussian curvature uses
  `grad(F)^T adj(Hess(F)) grad(F) / |grad(F)|^4`.
- Surface triangle area is always nonnegative.
- Raw mean-curvature sign depends on the native `grad(F)` orientation.
- Generator sorting in `build_surfs()` makes the raw convention deterministic
  for ordinary surfaces, but it is not a molecular-group outward convention.
- Gaussian curvature, squared mean curvature, and representative surface
  energy are invariant under normal reversal.
- Group-oriented mean curvature exists separately and is the correct starting
  convention for signed boundary-edge angles.

## Sampling dependence

All three canonical integrals are centroid quadratures over the generated
surface mesh. Their accuracy depends on `surf_res`, perimeter construction,
triangulation, projection accuracy, and treatment of degenerate triangles.
Flat surfaces bypass this quadrature. No convergence test across surface
resolutions currently protects these descriptors.

## Backward-compatibility surface

The following names must retain their current meanings until the complete
boundary result is explicitly versioned and validated:

- `mean_curv`, `avg_mean_curv`, `mean_tri_curvs`
- `gauss_curv`, `avg_gauss_curv`, `gauss_tri_curvs`
- `int_mean_curv`, `int_mean_curv_sq`, `int_gauss_curv`
- `surf_energy`
- CSV headings for integrated curvature and representative surface energy
- CLI/GUI/coloring aliases for the four canonical integrated fields

During development, new contributions should use explicit component names such
as `int_mean_curv_surface` and `int_mean_curv_edge`. The existing canonical
field should not become their total until analytic validation and an announced
compatibility transition.

## Existing test coverage

The supplied tests cover low-level surface, edge, network, log, and mesh
behavior, but contain no direct analytic or convergence tests for the current
mean/Gaussian integrations. `test_mesh.py` only checks that `surf_energy` can be
exported as a mesh property.

Added in this package:

- analytic sphere checks for mean and Gaussian curvature
- flat-plane behavior
- normal reversal behavior
- equivalence between the new quadratic API and the current scalar functions
- explicit failure for an undefined normal

## Performance-sensitive paths

`calc_surf_tri_curvs_both()` is called once for every curved surface and loops
over every surface triangle. It deliberately combines mean curvature, Gaussian
curvature, squared mean curvature, and triangle area in one pass. The new API
must not replace that optimized loop with per-triangle Python object dispatch.
It is intended initially for analytic edge quadrature and validation.

Full edge geometry should be computed only for boundary/barrier edges selected
for a group, not every internal network edge.

## Implemented foundation

`src/calculations/surface_geometry.py` now supplies:

- abstract `SurfaceGeometry`
- `QuadraticSurfaceGeometry`
- `PlaneSurfaceGeometry`
- implicit value, gradient, oriented unit normal, mean curvature, and Gaussian
  curvature methods

These are exported from `src/calculations/__init__.py`. No existing caller was
changed, so this addition does not alter logs or calculated values.

The second foundation stage also supplies:

- abstract `EdgeGeometry`
- exact `LineEdgeGeometry`
- validated `ParametricEdgeGeometry` adapter
- arc length, intrinsic curvature, and squared-curvature quadrature
- deterministic `1_group_2_external` / `2_group_1_external` classification
- exact selection of the two group/external surface views
- group-to-external surface-normal orientation
- robust signed surface-intersection angle using `atan2`

These additions are likewise not connected to canonical output fields yet.

## Next integration step

The edge sampler obtains rich construction data from `calc_edge_dir()` and
stores it in the existing `net.edges["vals"]` column alongside sampled points
and polyline length. Curved edges retain case-specific direction, center, and
projection metadata; straight edges currently retain only center and radius.
Before edge-curvature formulas are added:

1. Normalize the retained `vals` metadata without changing existing `points`,
   `vals`, or `length` columns.
2. Extend the new `EdgeGeometry` protocol with conic parameterizations that reproduce current
   endpoints for every AW edge case and represents straight Power/Primitive
   edges exactly.
3. Add endpoint, generator-order, rigid-transform, and arc-length tests.
4. Add deterministic group-relative boundary classification using the existing
   `edge.balls` and `edge.surfs` topology.
5. Only then implement signed dihedral quadrature and store it first as
   `int_mean_curv_edge` alongside—not inside—the canonical field.

Gaussian edge and vertex terms remain deferred until a validated
Gauss–Bonnet/normal-space formulation is selected.
