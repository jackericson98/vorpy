# Independent open-interface curvature validation

This work does not integrate with or redesign Phase 1. Production solvers,
filtrations, curvature analyzers, InterfaceGeometryAnalysis and exports remain
unchanged. The independent module is `vorpy.src.analyze.open_interface_curvature`.

## Definitions and mean curvature

Assume a compact orientable manifold with boundary, assembled from regular
smooth disk patches with matching seam positions and arc lengths. A topological
audit must establish these conditions before using the formulas. Use
`H=(k1+k2)/2`; synthetic examples use `Dn` so an outward sphere has positive H.
Static results retain the existing Group 1 -> Group 2 sign convention.

Round an internal crease using a thin cylindrical strip of radius epsilon.
Its transverse principal curvature is signed `1/epsilon`, its area element is
`epsilon d(angle) ds`, and longitudinal bounded terms vanish as epsilon tends
to zero. Its limiting scalar integral is therefore `1/2 integral beta ds`.
Vertex smoothing caps have area O(epsilon squared) and H O(1/epsilon), hence
no finite point mass. For this scalar normal-turn/smoothing convention:

`C_H = sum_faces integral H dA + 1/2 sum_internal_seams integral beta ds`.

A free boundary supplies no jump between two selected tangent planes and has
no additional scalar H mass. First variation does contain a boundary conormal
force. It is not an omitted scalar H contribution. Rounding a solid closure or
the boundary of a thickened sheet defines a different object. Also distinguish
scalar angle-times-length curvature from vector first-variation curvature,
whose crease magnitude involves `2 sin(beta/2)`.
These distinct scalar and vector constructions are discussed in
[Sullivan, sections 4.3–4.4](https://arxiv.org/html/0710.4497).

## Intrinsic curvature and Gauss–Bonnet

Orient each face boundary with the face on the left and set
`kg = (dT/ds) dot (n cross T)`. Reversal changes both n and the induced T;
kg is invariant. For an internal seam let
`q_e = kg_left_face + kg_right_face`, each evaluated in its own induced
boundary direction. This is the intrinsic line curvature measure. It can be
nonzero on curved glued patches. The extrinsic dihedral angle is not q_e.

At interior vertices let `d_v = 2pi - sum incident face-sector angles`.
At boundary vertices let `t_v = pi - sum incident face-sector angles`.
The complete identity is

`smooth_K + sum_internal integral q_e ds + sum_interior d_v`
`+ sum_free_boundary integral kg ds + sum_boundary t_v = 2pi chi`.

The first three terms are intrinsic Gaussian curvature of the surface.
The last two are exterior boundary terms. The Gauss–Bonnet left-hand side is
not interchangeable with the surface's integrated Gaussian curvature.
For background on the boundary theorem see
[UNC differential geometry notes, section 7](https://lindagreen.web.unc.edu/wp-content/uploads/sites/5262/2020/12/Section7_GaussBonnetTheorem.pdf).
For curvature measures on curved seams see
[Sullivan, section 4.1](https://arxiv.org/html/0710.4497).

An explicit cell-count derivation prevents corner double counting. Summing
individual disk-face Gauss–Bonnet identities gives 2pi F. Their corner sum is
`sum_v (n_v*pi - sum_sector_angles)`. Replacing it with assembled d_v/t_v
changes that sum by `2pi V_i + pi V_b - pi(2E_i+E_b)`.
The boundary is a union of cycles, so V_b=E_b. The change equals 2pi(V-E).
Consequently the assembled identity has right side 2pi(V-E+F), not 2pi F.
Do not sum face exterior corners and assembled defects together.

For planar faces all straight-edge kg terms vanish, including folded seams.
Only interior vertex defects contribute intrinsic K. A disk may have K=0
and boundary turning 2pi; folding a straight seam changes H but neither term.

## Topology and numerical validation

Reconstruct ordered physical face cycles. Require edge incidence one or two,
connected vertex links that are a path at the boundary or a cycle internally,
regular boundary cycles, and compatible orientation constraints. Then count
V-E+F and test each component using `chi=2-2g-b`. Every face must be a disk;
faces with holes require their own cell decomposition. Nonmanifold or pinched
complexes do not receive genus or curvature certification.

Power polygon coordinates map uniquely to cached solved dual vertex identities
(the matching coordinates have zero distance in this fixture). Every physical
edge is checked against the shared dual triangle of its two endpoints. The
mean-curvature audit independently reconstructs each edge's normal turn from
generator-oriented normals, adapting the sign to the existing production
convention, and compares all 835 contributions with the saved analyzer output.
The convex planar polygon adapter is not a general concave-face angle routine.

The known-answer cases are a polygonal disk, hemisphere, spherical latitude
cap, triangulated flat disk and folded disk. Sphere checks independently
integrate the parameterization area element with Gauss quadrature; the fold
check reconstructs the signed angle from normals instead of passing its
expected value into the accounting routine. Orientation reversal, a disk fan,
annulus, pinched vertex, three-face edge and unresolved/nonfinite mean terms
are also tested. Synthetic checks precede the static audit.

Certification here means complete supported incidence/contributions with
numerical consistency at the stated tolerance, not an interval-arithmetic
proof of floating-point error. A small global residual alone cannot certify
each contribution. The AW audit also exports per-face residual diagnostics
and a quadrature-order sweep. No residual is tuned away or assigned to a
missing term. A supported subcomplex has new boundaries and must not be used
as a substitute for the entire selected interface.

## Static matched-probe 2KAI findings

| Quantity | Power | AW |
|---|---:|---:|
| Selected pairs | 350 | 339 |
| Area (angstrom squared) | 898.118579 | 787.092578 |
| Full face-complex V, E, F | 605, 955, 350 | 584, 924, 339 |
| Euler characteristic | 0 | -1 |
| Components / boundary loops | 2 / 4 | 2 / 5 |
| Genus of each orientable component | 0 | 0 |
| Smooth H (angstrom) | 0 | 1.621963 partial |
| Internal-crease H (angstrom) | -64.450496 | -55.772124 partial |
| Free-boundary H | 0 | 0 |
| Complete H | -64.450496 | unavailable |

Power's components have (chi,b) equal to (-1,3) and (1,1). AW's full
face-complex components have (-2,4) and (1,1). The AW supported subset has
334 faces, V=571, E=906, chi=-1, one component and three boundary loops.
These counts concern physical face cells, not display triangles.

Power smooth K, internal-seam K and free-boundary geodesic curvature are zero.
Its interior vertex defects sum to **5.188337486733097**; its boundary corner
turns sum to **-5.188337486733230**. Thus its intrinsic integrated Gaussian
curvature is 5.188337486733097, whereas its full Gauss–Bonnet left-hand side
is -1.3234e-13, with target 2pi chi=0. Complete numerical accounting passes.
The old target used 605-955-350=-700, giving -4398.229715; that target was
incorrect. A sum of all boundary/interior corner terms near zero was a
Gauss–Bonnet sum, not an intrinsic Gaussian total of zero.

AW's full H is unresolved: five smooth faces and 14 internal mean-curvature
edges are missing. Its supported H contributions sum to -54.1501609113.
The Gaussian experiment separately evaluates the **334-face subcomplex**:

| Contribution at quadrature order 256 | Value |
|---|---:|
| Smooth K | 0.315254503411232 |
| Internal seam face-side geodesic curvature | -0.444077542436740 |
| Interior vertex defects | 4.776016756751914 |
| Exterior boundary geodesic curvature | 0.175166071574653 |
| Boundary corner turns | -11.106210591726630 |
| Gauss–Bonnet left-hand side | -6.283850802425571 |
| 2pi chi | -6.283185307179586 |
| Residual | -0.000665495245984 |

Residuals at orders 64, 128 and 256 are approximately -0.000773397,
-0.000686816 and -0.000665495. They exceed the fixed tolerance
`1e-7*max(1,abs(2pi chi))`. The supported-subcomplex experiment therefore
does not certify Gaussian accounting. Its artificial subset boundary terms
are not the full selected interface's boundary terms. Full AW intrinsic K,
full Gaussian boundary decomposition and full residual remain unavailable.
No omitted term is replaced by zero, no topology target is used to manufacture
missing curvature, and no quadrature or alpha threshold is fitted to the fixture.

All 334 supported disk-face diagnostics have computable terms. Their individual
Gauss–Bonnet residuals sum to -0.000665495245958, matching the assembled
residual. The largest absolute face residual is 0.000069047275080 on surface
10498. Thus the remaining discrepancy is already visible in local face
accounting; it is not explained by an incorrect assembled Euler target. This
does not identify which smooth/edge/angle numerical contribution is responsible.
The existing analyzers are left unchanged.

Final focused validation: **20 independent mathematical tests**; **68 tests**
in the combined mathematical, Phase 1, water, Power/AW and frozen Cazals suite.
Independent Power normal turns match all 835 existing internal-edge contributions
within 1.4433e-15 angstrom. All static checks reuse saved geometry; zero solves.

## Production proposal (not implemented)

Keep Phase 1 and its public result schema frozen. Add an optional independent
audit consumer that takes the same selected physical faces and existing
curvature outputs, plus solved face/edge/vertex incidence. It would gate H
certification on complete internal seams and manifold incidence; exposed
boundaries would not require a fictitious mean-curvature contribution.

Gaussian auditing must expose intrinsic K and Gauss–Bonnet boundary terms
separately, with per-feature support and residual/error diagnostics. Reuse
existing AW geodesic-curvature and vertex-angle functions through adapters.
For now AW remains unresolved. No production certification flags should change
on the basis of this static experiment, and no missing patch should be inferred
from the topology target. This proposal requires a separate integration task.

## Reproduction

From the repository root (the local pytest timeout plugin is unavailable):

```powershell
python -m pytest -o addopts='' vorpy/tests/analyze/test_open_interface_curvature.py -q
python -m pytest -o addopts='' vorpy/tests/analyze/test_open_interface_curvature.py vorpy/tests/analyze/test_power_interface_curvature.py vorpy/tests/analyze/test_aw_interface_curvature.py vorpy/tests/analyze/test_cazals_2kai_curvature_forensics.py vorpy/tests/analyze/test_cazals_validation.py vorpy/tests/interface/test_geometry_analysis.py vorpy/tests/interface/test_interface_water.py -q
python -m vorpy.src.analyze.open_interface_curvature --output-dir output/open_interface_validation/2KAI
```

The audit reads the saved matched-probe pair selection, Power polygon cache,
existing curvature component CSVs and solved `full_aw.vpy`. It performs zero
network solves. Results are written separately under the specified directory.
