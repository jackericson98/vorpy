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

## G2b production seam contribution

The AW interface path now promotes only a mathematically validated regular
analytic seam to an edge mean-curvature contribution.  It reuses the exact
line, regular clearance-branch, or nonsingular conic parameterization from
the AW edge resolver, and evaluates the two retained bicolor incident
surfaces with the validated A-to-B normal convention.  For an oriented edge,

`C_e = 1/2 integral(beta ds)`

with `H=(k1+k2)/2`; the per-edge record keeps the raw `integral(beta ds)`
and the half-scaled signed contribution separately.  Exactly two retained
bicolor surfaces are required.  Unsupported Apollonius incidence,
incomplete generator cells or edge features, degenerate parameterizations,
and failed quadrature remain unresolved (`None`), never zero.  Edge
contribution completeness is reported independently from total mean-curvature
certification; the latter still requires the higher-level selected-surface
topology and smooth-term accounting audit.

## G2c interface-level integrated mean curvature

The AW result now exposes one signed accounting record after all retained edge
records are measured.  It keeps the smooth-patch integral, the signed seam
integral (with raw `integral(beta ds)` retained separately), and the partial
sum distinct.  All integrated mean-curvature values use angstroms: `H dA` and
`C_e` both have units of length because radians are dimensionless.  Provenance
is the solved AW analytic surface records and validated A-to-B edge records.

Free-boundary and corner terms are explicitly not applicable to the scalar
mean-curvature measure; boundary conormals and corner turning terms belong to
the separate intrinsic/Gauss--Bonnet audit.  Missing or unsupported smooth
patches and seams remain unresolved and are not converted to zero.  The
standalone AW aggregate therefore reports `total=None` and `certified=False`;
only the topology-aware interface aggregation may certify a total after every
required smooth and internal-seam contribution is accounted for.

## G2d analytic validation

Focused synthetic validation checks the independent geometry identities.  For
the equal-radius straight seam, `L=2 A`, `beta=pi/2 rad`, so the raw seam
integral is `pi A rad` and `C_e=pi/2 A`; reversing A and B negates the signed
values while preserving length and absolute beta.  Independent cube and
concave-ridge checks give the expected positive convex and negative concave
signs under the same right-hand orientation, and the sphere patch check uses
`integral(H dA)=4 pi R A` with `H=(k1+k2)/2`.

Regular clearance-branch and nonsingular-conic AW examples are integrated by
the production analyzer and compared with an independent 256-node
Gauss--Legendre evaluation of the analytic point/tangent parameterization.
Tests require converged quadrature, finite angstrom-length and angstrom-radian
terms, and the half-scaled seam value.  Missing incident surfaces or
surface-referenced edge rows, unsupported/degenerate geometry, and
maximum-order quadrature remain unresolved; they do not contribute zero.  The
aggregate remains uncertified because free-boundary/corner intrinsic terms and
full selected-interface topology are owned by the separate accounting gate.

## G3 topology and corner accounting

For the scalar measure with `H=(k1+k2)/2`, a certified total requires:

1. every selected smooth patch to have a finite oriented `integral(H dA)`;
2. every internal seam to have exactly two regular selected incident patches,
   a finite signed `integral(beta ds)`, and the `1/2` seam factor;
3. every selected face to have one valid physical boundary cycle, with edge
   incidence one or two, path/cycle vertex links, regular boundary cycles, and
   a compatible global orientation; and
4. no missing, nonfinite, unsupported, degenerate, nonmanifold, pinched, or
   unresolved feature in the selected incidence.

An open interface is complete when its regular one-face boundary edges are
classified and its internal seams and smooth patches are complete.  Its free
boundary contributes no scalar H mass.  Boundary conormals, boundary geodesic
curvature, and boundary corner turns belong to first variation or the
intrinsic Gauss--Bonnet identity.  Interior vertex defects and seam endpoint
turns are likewise Gaussian/topological terms, not substitutes for mean
curvature.  A closed interface has no boundary terms but still requires the
same manifold, orientability, smooth-patch, and internal-seam checks.

The minimal open-disk and closed-tetrahedron checks certify zero scalar H while
separately retaining the open disk's `2pi` Gaussian boundary-corner term.  A
pinched vertex or edge with incidence other than one or two remains
nonmanifold and cannot certify a total.  The current full AW molecular result
therefore remains blocked by incomplete smooth/seam coverage and unresolved
full-interface topology/Gauss--Bonnet residuals; no Gaussian angle defect is
added to, or used to complete, integrated mean curvature.

Keep Phase 1 and its public result schema frozen. Add an optional independent
audit consumer that takes the same selected physical faces and existing
curvature outputs, plus solved face/edge/vertex incidence. It would gate H
certification on complete internal seams and manifold incidence; exposed
boundaries would not require a fictitious mean-curvature contribution.

Gaussian auditing must expose intrinsic K and Gauss–Bonnet boundary terms
separately, with per-feature support and residual/error diagnostics. Reuse
existing AW geodesic-curvature and vertex-angle functions through adapters.
AW seam aggregation is production-supported only for the validated regular edge
classes above.  No certification is inferred from a topology target when a
required patch or seam is missing.

## G4 EDTA molecular completeness audit

A bounded audit was run on the existing `vorpy/data/EDTA.pdb` AW solve, using
generator IDs A=`0..31` (EDTA) and B=`{32}` (Mg). It reused the solved network
and did not construct an additional network. The network contains 612 stored
balls, 424 surfaces, 892 edges and 500 vertices; 32 generator cells are
complete and 580 are incomplete. The global Apollonius audit reports 1,474
missing face relations, all attributed to incomplete network geometry.

The EDTA--Mg interface has 6 candidate smooth patches and 15 candidate primal
edges. Smooth coverage is 0/6 eligible: all six stored `integral(H dA)` values
are finite, but every participating generator-cell proof is incomplete. Of
the 15 edges, 11 retain exactly two bicolor surfaces and 4 have one retained
surface (open-boundary incidence); no surface-referenced edge row is missing.
All 15 edge parameterizations independently resolve as `matched_nonsingular`,
so no unsupported analytic edge class was observed. The 11 exact-two seams all
converge at order 16. Their total length is 9.2343944860 A, raw signed
`integral(beta ds)` is -10.2380670445 A rad, and the conventional seam
contribution is `C_e=-5.1190335222 A`. Maximum reported quadrature change is
8.67e-12 A rad; A-to-B reversal changes the signed result by sign only.
Radians are dimensionless, so both `H dA` and `C_e` are lengths.

The six-patch incidence is one connected orientable regular open disk:
V/E/F=`10/15/6`, chi=`1`, four boundary edges, one boundary loop, and no
nonmanifold or pinched feature. Thus topology is not the current EDTA blocker;
construction completeness is. Curvature analysis has no unsupported or
non-converged analytic edge blocker, but cannot promote the numerically
available seams while their generator-cell proofs are incomplete. Boundary
edges receive no fictitious scalar mean-curvature term; their conormal,
geodesic and corner quantities remain separate Gaussian/Gauss--Bonnet
accounting. The EDTA total therefore remains unresolved, not zero or
CERTIFIED.

Prioritized remediation: (1) complete and validate the Mg-participating AW
generator cells and missing face relations; (2) re-audit exact-two seam
incidence and open-boundary cycles; (3) run the existing analytic seam and
orientation regressions, then the topology/Gauss--Bonnet residual audit; and
(4) certify only when all smooth patches, internal seams, boundary
classification, and selected-interface manifold/orientation checks pass.
The focused AW, orientation, edge-geometry, open-topology, interface-analysis,
and EDTA Apollonius regressions passed; no broad suite was run.

## G5 local AW interface completeness

Global generator-cell completeness is stronger than local patch completeness.
The former proves every participating AW cell and all higher-order closure
relations in the network. A scoped patch audit need only prove the selected
physical faces, their selected edge incidence, and the analytic geometry used
by those faces and seams. It must still report global closure as unresolved
and must never promote a local result to a global interface certificate.

The minimum local evidence is: finite bounded selected surface records with
valid one-cycle boundaries; every selected pair face present; every internal
edge referenced by exactly two selected faces and every exposed boundary edge
referenced by exactly one; finite endpoint/clearance residuals from the
canonical analytic edge resolver; converged signed seam quadrature; and a
manifold, orientable selected topology with valid vertex links. Missing
relations are partitioned by scope: a missing selected pair, seam coface, or
referenced edge is a local blocker, while a missing non-interface coface is an
out-of-scope global blocker that must remain visible.

The focused EDTA diagnostic confirms the distinction. None of the 1,474 global
missing-face relations names one of the six selected EDTA--Mg pair faces or
one of the 11 exact-two seam triples. Four relations touch candidate boundary
triangles `(5,32,330)`, `(9,32,330)`, `(11,32,330)`, and `(15,32,330)`; each
only lacks the non-interface pair `(32,330)`. The selected patch remains the
validated orientable open disk from G4. This supports a scoped status such as
`LOCAL_COMPLETE_GLOBAL_UNRESOLVED` when the local evidence passes, but the
current production status remains unresolved because all six patches still
fail the global generator-cell gate.

The regression is
`test_edta_interface_local_faces_are_distinct_from_global_closure_gaps`.

## G6 molecule-only AW completeness

F10 is now the authoritative producer for this diagnostic. Its default
`MOLECULE_ONLY` solve excludes the physical-solvent generator universe before
AW construction; group selection does not perform a second filtering solve.
Using the existing EDTA input with `boundary_mode="none"` produced 32 retained
generators, 76 surfaces, 174 edges and 101 vertices. Two generator cells are
complete and 30 are incomplete; 484 missing face relations remain. The input's
Mg row is not in this producer's molecule-only generator membership, so an
EDTA--Mg claim is not inferred here. If Mg must be a molecular interface
partner, its F10 membership/provenance classification is an upstream blocker,
not a curvature value to fabricate.

For a bounded local diagnostic, the 11 bicolor surfaces incident to retained
generator 0 form a closed orientable sphere (`V/E/F=18/27/11`, chi=2), with 27
exact-two selected edge incidences and no boundary. The selected pair faces
are present in the incidence table, while missing cofaces occur between other
retained molecular generators. Thus solvent exclusion is demonstrated by
membership evidence; the remaining missing geometry is finite/incomplete
molecular-cell geometry, not a solvent contribution.

The existing conservative AW local certificate was also run on this same
network without rebuilding it. Its extracted view contains the expected 11
surfaces, but the shared dual support component reaches all 101 vertices, 174
edges, 76 surfaces and 32 generators. It therefore returns
`proven=false`, `incomplete_relevant_competitor`; the selected patch is not
promoted to a local geometric certificate merely because unrelated generators
are outside the requested group. Narrowing that support component would be a
new mathematical completeness rule and remains a GAUSS blocker.

Curvature coverage is analytically complete for the selected edge geometry,
but not globally support-complete. Fifteen edges resolve directly as
`matched_nonsingular` conics and twelve as regular `matched_curve` clearance
branches. The latter are independently rematched to nonsingular conics before
quadrature when the clearance parameter is poorly conditioned; all 27 converge
under order 8--64 quadrature. Their length sum is 26.3211533485 A, raw signed
`integral(beta ds)` is -31.7816569569 A rad, and the H-convention seam sum is
-15.8908284784 A (`C_e=0.5*integral(beta ds)`). The largest edge quadrature
change is 2.13e-10 A rad and A/B reversal negates each signed result. Five
edges still fail the global Apollonius support test because their non-interface
cofaces are missing; analytic validity does not erase that topology blocker.

Independent AW Gauss--Bonnet re-integration resolves the selected faces,
incidences and vertex angles, but is not certified: after validated conic
rematching, the residual is -8.8835714837e-4 at orders 16 and 32 (also stable
at order 64). The quadrature blocker is therefore removed, while the persistent
nonzero residual remains an unresolved geometry/accounting discrepancy. Local
topology is diagnostic evidence only; neither angle defects nor Gaussian terms
substitute for missing mean-curvature contributions. Global curvature
certification remains fail-closed, with H=(k1+k2)/2 and
`C_e=0.5*integral(beta ds)` unchanged.

The focused regressions are
`test_edta_interface_local_faces_are_distinct_from_global_closure_gaps` and
`test_edta_molecule_only_interface_is_local_complete_but_global_unresolved`.

## G9 analytic curvature and Gaussian certification

The clearance-parameterized `matched_curve` class is a mathematically regular
AW branch whenever its stored interval avoids a conic turning point. Near a
turning point its clearance parameter can make Gauss--Legendre convergence
poor even though the curve is exact. GAUSS now rematches such a branch to a
nonsingular conic only after validating endpoints and all stored samples;
failed rematching or quadrature remains unresolved. Resolver status and the
signed convention are unchanged. The EDTA local audit consequently integrates
all 27 selected regular edges (15 direct conics and 12 rematched branches),
with maximum edge quadrature change `2.13e-10 A rad`.

Independent reference regressions now cover a radius-2 sphere
(`integrated K=4 pi`, `integrated H=8 pi`) and a side-2 convex cube. The cube
has `V/E/F=8/12/6`, smooth K=0, zero face-boundary geodesic K, vertex defect
`4 pi`, and piecewise-smooth mean total `6 pi` from twelve length-2 edges with
exterior turn `pi/2`. This keeps Gaussian vertex/edge accounting distinct from
the extrinsic mean-curvature seam sum.

The EDTA Gaussian audit is still not certified: selected-face, edge-incidence
and vertex-angle evaluation has no unresolved primitive rows, but the
Gauss--Bonnet residual is stably `-8.8835714837e-4` at orders 16, 32 and 64,
above the certification tolerance. Conic rematching has removed the prior
clearance-branch convergence defect; the remaining blockers are global
Apollonius support for five cofaces, incomplete smooth-face/cell provenance,
and the persistent geometry/accounting residual. At order 32 the explicit
decomposition is smooth-face `K=0.12690830325`, seam edge/geodesic
`K=1.99521641743`, interior vertex defect `K=10.44335753654`, and zero free
boundary/corner terms; the sum is `12.56548225721` versus `2*pi=12.56637061436`.
This is an accounting diagnostic, not a certified total: no total is promoted
to CERTIFIED.

## G10 Gauss--Bonnet residual root-cause audit

For the selected closed orientable patch, `chi=2`, so Gauss--Bonnet requires
`2*pi*chi=4*pi`, not `2*pi`: the latter is the disk/hemisphere value
(`chi=1`). The reference is independent of A/B orientation; reversing the
retained cell group leaves smooth, seam, vertex, total, and residual terms
unchanged.

The order-32 EDTA decomposition has smooth-face `K=0.12690830325`, internal
face-edge geodesic `K=1.99521641743`, interior vertex defect
`K=10.44335753654`, and no boundary terms. All 11 face records, 27 edges
(two selected incidences each), and 18 vertices are complete. Face residuals
localize the discrepancy, with the largest current terms on surfaces 1, 4,
and 7 (`-2.9537e-4`, `-2.8739e-4`, and `-2.9397e-4`); no global correction is
applied.

The residual is surface-resolution dependent while analytic seam and vertex
geometry remain fixed: `-8.88357e-4` at resolution `0.2`, `-2.38917e-4` at
`0.1`, and `-5.97456e-5` at `0.05`, with order-16/32 agreement at each tested
resolution. The current smooth-face value is accumulated from Gaussian
curvature at triangle centroids times flat triangle area, so this is the
expected finite-surface-mesh approximation error, not quadrature failure or a
sign/incidence defect.

As an independent topology check, merging the selected face triangulations by
coordinate gives 1,160 triangles, 582 vertices, 1,740 edges, no boundary or
nonmanifold edges, and discrete vertex defect `12.566370614359007`, agreeing
with `4*pi` to `1.7e-13`. This validates the sampled mesh's closed topology
but does not replace exact AW smooth-face provenance. The five missing global
Apollonius cofaces do not alter selected pair-face incidence; they remain
global completeness blockers only. Certification therefore remains
fail-closed until smooth-face integration/provenance is upgraded or bounded
error is formally accounted for.

## G11 accurate smooth-face Gaussian integration

The selected AW smooth faces are regular pairwise quadratic implicit surfaces.
Their existing metadata is sufficient for a curved integral: the 14-coefficient
implicit equation, the stored lifted mesh vertices and triangles, the projection
plane point, and the generator-center direction used to select the quadratic
branch. No producer metadata change was needed.

GAUSS now provides a curved-surface path for the AW audit. Each stored triangle
is projected to its existing plane, Duffy-mapped tensor Gauss--Legendre nodes
are lifted back to the validated quadratic branch, and the exact implicit
Gaussian curvature is multiplied by the implicit graph Jacobian. Stored
vertices are revalidated against the same branch; missing, singular, or
inconsistent data raises an unresolved-face condition. The previous centroid
curvature times flat-triangle-area calculation remains available as
`smooth_face_method="triangulated"` for baseline comparisons. The curved path
is the default, with order 4 in the EDTA audit and explicit higher orders for
convergence studies.

The spherical-cap regression converges in quadrature order 8 to the analytic
value `2*pi*(1-cos(theta))`; order 16 changes by less than `1e-6` for the
128-sided cap domain (the remaining `2e-3` test allowance is the polygonal
boundary approximation). Existing sphere and cube regressions retain their
analytic references: a closed sphere is `4*pi` Gaussian curvature and the
piecewise-flat cube obtains `4*pi` from vertex defects, with zero smooth-face
curvature.

On molecule-only EDTA, the curved method reduces the order-16 seam-quadrature
residual at surface resolution `0.2` from the centroid/flat baseline
`-8.883571e-4` to approximately `-1.364640e-5`; at resolution `0.1`, order-4
curved integration gives `-3.591374e-6` versus the baseline
`-2.389175e-4`. The order-8/16 EDTA smooth integral is unchanged to the shown
precision in the independent single-face study, and whole-patch seam
quadrature is stable across its tested orders, so the remaining error is finite
surface-domain resolution rather than quadrature non-convergence. The curved route evaluates four times as many
nodes per triangle at order 8 as order 4; this is the bounded cost/accuracy
tradeoff retained for explicit studies, while the audit default stays order 4.

These improvements are numerical only. Seam/geodesic and vertex terms are
unchanged, A/B reversal remains invariant for the Gaussian decomposition, and
the EDTA result is still `UNRESOLVED`: global Apollonius completeness and
smooth-face provenance are incomplete, and the residual is not zero. A small
residual or converged quadrature does not promote the total to `CERTIFIED`.

## G12 finite AW surface-domain error study

The bounded EDTA study holds the selected seam/vertex decomposition and seam
quadrature at order 16, then varies only the curved smooth-face quadrature. The
reported residual is therefore a Gauss--Bonnet residual, not a fabricated
curvature correction:

| surface resolution | smooth order 2 | smooth order 4 | independent order 8 |
| ---: | ---: | ---: | ---: |
| 0.2 | `-1.05642018e-5` | `-1.36463602e-5` | `-1.36464044e-5` |
| 0.1 | `-3.58353000e-6` | `-3.59137383e-6` | `-3.59137393e-6` |

Orders 4 and 8 differ by `4.42e-11` at resolution 0.2 and `9.77e-14` at
resolution 0.1. This is numerical quadrature convergence. The order-4
residual decreases by a factor `3.80` when the resolution is halved, giving
observed order `p=1.926`; the remaining residual is consequently dominated by
finite projected-domain/surface-mesh approximation rather than quadrature.
The empirical estimate is useful, but it is not a rigorous error bound: the
current meshes are not provenance-backed nested refinements, and exact AW face
boundary domains and global Apollonius cofaces remain incomplete.

For a tensor rule, the smooth integrator evaluates `q^2` nodes per projected
triangle: 4, 16, 64, and 256 nodes at orders 2, 4, 8, and 16. Thus order 8
costs four times order 4 for the smooth term. The q4/q8 two-resolution study
took approximately 70.6 seconds on the development workstation; the focused
regression uses q2/q4 to stay bounded while the independent q8 values remain
documented evidence.

The existing sphere-cap order-8/16 convergence and cube `4*pi` vertex reference
remain the analytic checks; neither provides a rigorous EDTA face-domain bound.
The next numerical improvement should be local adaptive subdivision of the
existing projected triangles with per-face q4/q8 differences and an explicit
unresolved status when the estimated domain error is above threshold. This
would improve diagnostics without changing exact quadric definitions or
weakening certification. Global completeness, smooth-face provenance, and a
validated error budget are still required before certification.

## G13 adaptive AW domain integration feasibility

The current AW representation supports a local prototype: project each stored
triangle into its existing plane, subdivide that planar triangle, and lift the
new vertices with the same validated quadratic branch. The outer polygon and
its boundary vertices remain unchanged, so this experiment changes quadrature
resolution only; it cannot remove finite boundary/domain approximation error.

On the validated eight-sided spherical-cap fixture, the fixed-domain values
were:

| mesh | q4 | q8 | q16 reference |
| --- | ---: | ---: | ---: |
| fixed, 8 triangles | `2.7070761993` | `2.7071271912` | `2.7071272523` |
| one projected subdivision, 32 triangles | `2.7071336864` | `2.7071272524` | — |

Relative to the fixed q16 value, fixed q4 has error `5.1053e-5`, while one
subdivision with q4 has error `6.4341e-6`, an eight-fold improvement. However,
the two alternatives have the same node count: fixed q8 uses `8*8^2=512`
nodes and subdivided q4 uses `32*4^2=512` nodes. Fixed q8 is much more accurate
on this smooth face (`6.11e-8` error), so subdivision is not a default
replacement for increasing quadrature order here. Subdivided q8 reaches about
`1e-10` against the reference at four times the node count again.

The q4/q8 difference is a useful local quadrature indicator in this fixture:
on the fixed mesh it tracks the q4-to-q16 error within about 0.2%, and on the
subdivided mesh it tracks the q4 error to the displayed precision. It still
does not estimate the unchanged polygonal/domain error. A future local
adaptive implementation should therefore stop only when every face's q4/q8
delta meets the caller's numerical error budget, all branch/regularity checks
pass, and an independent nested-domain or boundary-resolution comparison is
available. Otherwise the contribution remains unresolved. This prototype does
not justify a global adaptive engine or any certification change.

## G14 curvature integration validation and cost audit

The bounded EDTA comparison keeps the seam and vertex terms identical and
changes only the smooth-face method. Timings are wall-clock measurements from
the focused regression on the development workstation. The triangulated time
is the full fixed-term Gauss--Bonnet call (including seam/vertex assembly); the
curved time is the smooth-face q4 integration only, so it represents the
incremental curved cost when the other terms are reused.

| resolution | triangulated residual | triangulated call | curved q4 residual | curved smooth q4 |
| ---: | ---: | ---: | ---: | ---: |
| 0.2 | `-8.88357148e-4` | `0.255 s` | `-1.36463602e-5` | `3.018 s` |
| 0.1 | `-2.38917450e-4` | `0.243 s` | `-3.59137383e-6` | `10.773 s` |

The curved result is substantially closer to the Gauss--Bonnet reference at
both resolutions. Its cost grows with the number of projected triangles and
with the tensor node count `q^2`; the per-node quadratic branch root solve,
implicit gradient/Jacobian evaluation, and Gaussian-curvature calculation are
the dominant operations. No seam or vertex computation was changed.

The validated q4/q8 stopping data are `4.42e-11` at resolution 0.2 and
`9.77e-14` at 0.1, against residual magnitudes of `1.36e-5` and `3.59e-6`.
Thus q4 is numerically converged for these EDTA faces, while the G13 coarse
sphere fixture still demonstrates that the same indicator detects a
quadrature-limited face. The recommended production default remains curved q4,
with q8 as a diagnostic/escalation order when a per-face q4/q8 difference
exceeds the caller's numerical budget. This is a numerical stopping rule, not
a certification rule.

The sphere and cube regressions continue to validate analytic `K`/Gauss--Bonnet
accounting, and the EDTA tests continue to verify A/B orientation invariance,
complete selected incidence rows, and fail-closed global status. EDTA remains
uncertified because global Apollonius completeness, smooth-face provenance,
and a rigorous finite-domain error bound are still absent.

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
