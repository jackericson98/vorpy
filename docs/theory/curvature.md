# Curvature

VorPy describes each pairwise molecular boundary with an implicit surface
`F(x, y, z) = 0`.

For additively weighted (AW) networks, the pairwise boundary is generally a
quadratic surface. Power, Primitive, and equal-radius AW boundaries are flat.

## Local surface curvature

At a point on an implicit quadratic surface, VorPy evaluates mean curvature
`H` and Gaussian curvature `K` from the gradient and Hessian of `F`:

```text
H = [tr(Hess(F)) |grad(F)|^2 - grad(F)^T Hess(F) grad(F)]
    / [2 |grad(F)|^3]

K = grad(F)^T adj(Hess(F)) grad(F) / |grad(F)|^4.
```

The formulas are analytic. Their integration is numerical: VorPy evaluates
`H` and `K` at each surface-triangle centroid and weights them by that
triangle's area.

Current surface-patch descriptors are

```text
int_mean_curv     = sum(H_t A_t)       [angstrom]
int_mean_curv_sq  = sum(H_t^2 A_t)     [dimensionless]
int_gauss_curv    = sum(K_t A_t)       [dimensionless].
```

Here `t` indexes surface triangles. Results can depend on surface resolution,
perimeter construction, triangulation, projection accuracy, and degenerate
triangle handling. Resolution-convergence testing is therefore required for
quantitative use.

## Orientation

Reversing the surface normal reverses the sign of `H` but does not change `K`
or `H^2`. The raw per-surface `int_mean_curv` follows the native
implicit-function orientation; it is not automatically outward relative to an
arbitrary molecular group.

Group analysis separately determines an outward direction from the group-side
generator toward the external generator and stores
`oriented_int_mean_curv`. This convention is the basis for future signed edge
contributions.

## Complete piecewise-smooth boundary

A molecular boundary is piecewise smooth. Its curvature is distributed over
smooth surface patches, their junction edges, and their corner vertices. The
target complete descriptors are

```text
M_total = M_surface + M_edge
G_total = G_surface + G_edge + G_vertex.
```

For a boundary edge `e`, the candidate classical mean-curvature contribution
is

```text
M_edge,e = 1/2 integral_e theta(s) ds.
```

Here `theta` is the signed angle between the two outward-facing surface normals
about a canonically oriented edge tangent.

## Current compatibility status

The public fields `int_mean_curv` and `int_gauss_curv` currently contain
surface-patch quadrature only. Edge and vertex terms are not yet included.
During development, new contributions will use explicit component names:

- `int_mean_curv_surface`
- `int_mean_curv_edge`
- `int_mean_curv_total`
- `int_gauss_curv_surface`
- `int_gauss_curv_edge`
- `int_gauss_curv_vertex`
- `int_gauss_curv_total`

The canonical fields will not be redefined until complete-boundary formulas
have passed analytic, invariance, convergence, and molecular regression tests.
`int_mean_curv_sq` remains a smooth-surface bending descriptor and is not
silently augmented with edge curvature.

See [Morphometric Geometry Implementation](../development/morphometric_geometry.md)
for the equation-to-code map and validation status.
