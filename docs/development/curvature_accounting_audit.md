# Curvature Accounting Audit

This audit separates pointwise curvature from integrated boundary measures.
It records the current implementation state; it does not promote a partial
geometric measurement to a complete morphometric measurement.

## Conventions, orientation, and units

VorPy uses

```text
H = (k1 + k2) / 2                 [Å^-1]
K = k1 k2                         [Å^-2]
C = integral(H dA)                [Å]
X = integral(K dA)                [1]
```

For a quadratic implicit surface `F = 0`, the native normal is
`grad(F) / |grad(F)|`.  With that normal the analytic implementation evaluates

```text
H = (tr(Hess(F)) |grad(F)|^2 - grad(F)^T Hess(F) grad(F))
    / (2 |grad(F)|^3)
K = grad(F)^T adj(Hess(F)) grad(F) / |grad(F)|^4.
```

Reversing the normal reverses the sign of `H` and leaves `K` unchanged.  The
legacy AW surface path retains this native implicit-function orientation;
`oriented_int_mean_curv` is a separately derived group-relative quantity.  The
molecular-contact path uses an outward molecular normal.  Values from those
two orientation conventions must not be added without an explicit conversion.

A regular plane has measured `H = K = 0`.  A singular implicit field has no
normal and no curvature measurement.  `mean_curvature`, `gaussian_curvature`,
and their combined triangle quadrature now raise `ValueError` for that case;
they never substitute zero for an unavailable result.

## Implemented measurement paths

| Path | Pointwise inputs | Accumulation | What it measures | Completeness |
| --- | --- | --- | --- | --- |
| Legacy network surface | Analytic `H`, `K` at each triangle centroid | `sum(H_t A_t)`, `sum(H_t^2 A_t)`, `sum(K_t A_t)` | Smooth triangulated pairwise patches | Does not include physical seam, boundary, or corner terms |
| Canonical molecular contact surface | Exact spherical-patch area and physical refined seams | Smooth, seam, boundary, and junction components below | Selected union-of-balls contact boundary | Certified only after incidence/orientation and Gauss--Bonnet audits pass |

The legacy fields `int_mean_curv`, `int_mean_curv_sq`, and `int_gauss_curv`
are centroid quadratures, not pointwise curvature and not automatically
complete-boundary invariants.  `surf_energy = 2 * int_mean_curv_sq` is an
orientation-invariant reference bending score over smooth patches only.  It is
not a free energy and is not augmented with edge terms.

For the canonical molecular-contact surface, spherical patches use
`H = 1 / r` and `K = 1 / r^2` in the outward molecular orientation:

```text
C_smooth = sum_p A_p / r_p
C_seam   = 1/2 sum_e beta_e L_e
C_total  = C_smooth + C_seam.
```

`beta_e` is the signed dihedral angle of the two outward patch normals and
`L_e` is the physical refined-seam length.  An exterior free boundary has no
automatic scalar mean-curvature term, so its typed
`boundary_contribution` is `None / NOT_SUPPORTED`, rather than zero.

The molecular Gaussian audit keeps intrinsic and completed quantities
distinct:

```text
X_intrinsic = integral(K dA) + sum_internal_seams integral(k_g ds)
              + sum_interior_junctions defect

X_completed = X_intrinsic + sum_exterior_boundary integral(k_g ds)
              + sum_boundary_corners turning
            = 2 pi chi
```

The seam `k_g` terms are intrinsic geodesic curvature.  They are not the
extrinsic dihedral terms in `C_seam`.  Auxiliary topology cuts are excluded
from both physical curvature sums.  The typed
`integrated_gaussian_curvature` currently names `X_intrinsic`; its scope is
`intrinsic_surface`.  The Gauss--Bonnet-completed value is emitted separately
as `gauss_bonnet_completion`.

## Status rules

`CERTIFIED` requires the selected geometry, orientation, required incidence,
and relevant accounting audit to be complete.  `PARTIAL` may retain a known
component with its scope, but cannot certify the total.  `UNRESOLVED` means a
known obstruction prevents a value; `NOT_CALCULATED` means no applicable
producer ran; `NOT_SUPPORTED` documents a deliberately unavailable term (the
free-boundary scalar mean term); and `STALE` means a once-calculated value no
longer has valid provenance.  None of these states is represented by a zero
measurement.

The molecular-contact conversion already transfers
`total_area_A2` to `MolecularContactSurfaceResult.metrics["area"]` with units
`Å²`, the construction method, scope `molecular_contact_surface`, provenance,
and its independent `area_status`.  The conversion regression covers
`CERTIFIED`, `PARTIAL`, `STALE`, `UNRESOLVED`, `NOT_CALCULATED`, and
`NOT_SUPPORTED`; no Results-contract change is required for area population.

## Synthetic validation map

The existing lightweight fixtures in `geometry.validation_shapes` use the
same `H = (k1 + k2)/2` convention.

| Shape | Expected `C` | Expected `X` | Interpretation |
| --- | --- | --- | --- |
| Sphere of radius `R` | `4 pi R` | `4 pi` | Smooth closed genus-zero surface |
| Rectangular box `a x b x c` | `pi(a+b+c)` | `4 pi` | Smooth faces contribute zero; edge and corner completion are required |
| Cube of side `a` | `3 pi a` | `4 pi` | Box special case |
| Spherocylinder, radius `r`, cylinder length `L` | `pi L + 4 pi r` | `4 pi` | Smooth cylindrical and spherical terms |
| Torus, major radius `R` | `2 pi^2 R` | `0` | Geometric deformation with genus one |
| Deep pocket | no certified analytic `C` | `4 pi` | A recess changes geometry, not topology, while the boundary remains closed genus zero |
| Closed spherical internal cavity | `4 pi(R_outer - r_cavity)` for the spherical shell convention | `8 pi` | A new closed boundary component is a topology change |

The deep-pocket mean-curvature reference is deliberately `NaN`, never zero.
The cavity fixture checks component and orientation behavior separately.  The
full-sphere, closed equal-sphere union, boundary-hole, multitorus, deep-pocket,
and cavity regressions exercise the applicable topology and Gauss--Bonnet
paths; a general production sharp-polyhedron completion is still not present
in the legacy network path.

## Morphometric free-energy framework

The proposed quantity is a separate calibrated layer:

```text
Delta G = p V + gamma A + kappa C + kappa_bar X.
```

`V`, `A`, `C`, and `X` are geometric measurements.  In Å-based coordinates,
the supplied coefficients must carry energy/Å^3, energy/Å^2, energy/Å, and
energy units respectively.  A future evaluator must accept coefficient values,
units, fit dataset/version, temperature/reference state, and provenance as
external inputs; it must not fit or invent them from geometry.

It may calculate a free-energy result only when all four selected geometric
terms are compatible, complete for the declared model, and `CERTIFIED`.  If a
selected term is partial, unresolved, unsupported, not calculated, or stale,
the free-energy value must remain null with the corresponding status.  The
current `surf_energy` does not satisfy these requirements and must continue to
be reported as a representative surface-energy/bending descriptor rather than
a validated free-energy predictor.
