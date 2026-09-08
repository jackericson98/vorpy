# Additively Weighted Voronoi

For a generator centered at `c_i` with radius `r_i`, the additively weighted
distance is

```text
d_AW(x, i) = ||x - c_i|| - r_i.
```

A point belongs to the generator with the smallest surface distance. Pairwise
boundaries satisfy

```text
||x - c_i|| - r_i = ||x - c_j|| - r_j
```

and are generally quadratic surfaces. An equal-radius pair reduces to a flat
bisector.

## Three-generator edges

An AW edge is the locus with equal clearance `rho` from three generators:

```text
||x - c_0|| - r_0
  = ||x - c_1|| - r_1
  = ||x - c_2|| - r_2
  = rho.
```

Subtracting the squared equations produces two constraints linear in `x` and
`rho`. Solving those constraints and substituting back into one sphere equation
gives an analytic planar conic branch:

```text
x(rho) = x0 + u rho +/- n sqrt(q(rho)),
```

where `n` is normal to the generator plane and `q` is quadratic. The sign
selects one conic branch. The coefficient of `rho^2` distinguishes ellipse,
parabola, and hyperbola forms.

The current analytic API intentionally rejects an interval that reaches a root
of `q`. The `rho` parameter is singular at such a turning point, even when the
geometric curve itself is smooth. Production connection therefore requires
splitting the edge or introducing a nonsingular angle/hyperbolic
parameterization.

## Matching stored network edges

`match_aw_trisector_to_endpoints` recovers the clearance at both stored
vertices, selects the conic branch from the signed generator-plane normal
coordinate, and checks both endpoint reconstruction and equal-clearance
residuals. The returned curve is affinely parameterized on `[0, 1]`, from the
lower vertex index to the higher vertex index, even when clearance decreases in
that direction.

The matcher reports explicit diagnostic statuses:

- `matched`: both vertices lie on one analytic branch within tolerance;
- `cross_branch`: the vertices select opposite conic branches;
- `turning_point`: clearance cannot parameterize the complete interval;
- `residual_too_large`: reconstruction or equal-clearance validation failed;
- `invalid`: the analytic trisector could not be constructed.

`validate_edge_samples` then projects each legacy sample through its recovered
clearance and reports maximum and RMS reconstruction errors. These diagnostics
do not yet replace `edge.points` or alter production descriptors.

At network level, `diagnose_aw_edge_geometry` applies that validation to the
constructed `net.edges` table. Equal-radius AW triples are checked against their
exact straight-line limit; unequal-radius triples are checked against the
analytic conic. It reports endpoint, clearance, sample-path, and length
residuals per edge, plus aggregate match coverage. Unsupported topology and
individual numerical failures remain explicit records so one problematic edge
does not abort the audit.

Implementation: `calculations/edge_geometry.py::AdditivelyWeightedTrisectorBranch`.
Validation: `tests/calculations/test_edge_geometry.py` checks equal-clearance
residuals, derivatives, intrinsic curvature integration, rigid transforms, and
turning-point rejection. The same files implement and test endpoint matching,
canonical orientation, branch rejection, and sampled-path residuals.
Network-table integration is validated in
`tests/network/test_build_net.py` without mutating the input DataFrames.

VorPy network name: `aw`.
