# AW geometry producer audit and envelope

## Scope and units

`vorpy.src.network.aw_geometry_envelope.extract_aw_geometry_envelope` is a
read-only producer adapter for a completed additively weighted (AW) network.
Coordinates, radii, and analytic edge parameters are in the `Network` input
coordinate unit (Å for normal molecular inputs).  It calculates no area,
integrated mean curvature, integrated Gaussian curvature, Euler
characteristic, or boundary/corner measurement.  Those remain GAUSS/ATLAS
responsibilities.

The envelope is deliberately separate from the Results model.  It is an
internal FORGE hand-off: a downstream consumer may use its stable IDs and
incidence, but must attach its own quantities, scientific status, and Results
provenance through the ATLAS-owned contract.

## Producer audit

| Requirement | Existing producer state | Envelope state |
| --- | --- | --- |
| Atoms and cells | `balls` stores the generator centre/radius and `complete` flag. | Atom and one-to-one cell records have deterministic IDs and full incidence. |
| Pairwise surfaces | `surfs` stores two generator indices, perimeter topology, mesh, and usually `func` (14 AW quadratic coefficients).  `surface_geometry` already evaluates that exact implicit surface. | Coefficients are preserved; if `func` is absent they are reconstructed exactly from the two generators.  No mesh fit is used. |
| Spherical surface patches | Not an AW primitive: pairwise AW boundaries are planes/quadrics, not atom-carried spherical patches. | Every AW face says `spherical_patch_status=NOT_SUPPORTED`; it is never relabelled as a sphere. |
| Boundary edges/arcs | The canonical edge resolver supplies an exact line, AW trisector branch, or nonsingular AW conic when it can prove one. | Reconstructible parameters are emitted for those exact cases.  Resolver failures are `UNRESOLVED` with no polyline fallback. |
| Junction vertices | `verts` stores defining balls, location, clearance radius, and doublet information. | Stable junction IDs and an AW-clearance residual check are emitted. Invalid/nonregular records are `UNRESOLVED`. |
| Connected components | Build diagnostics can report counts but do not expose stable component records. | Stable components are computed from explicit surface/edge/vertex incidence and carry boundary-incidence diagnostics. |
| Surface ownership | A surface has an unordered generator pair but no first-class directed views. | Each face has two owner records. A normal means `aw_clearance_gradient(owner_cell, neighbor_cell)`, so its direction is explicit and point-dependent. |
| Orientation and normals | `aw_outward_normal`, surface-mean-curvature orientation, and canonical perimeter orientation already exist. | The envelope references the existing exact clearance-gradient convention rather than sampling a normal or choosing an arbitrary pair order. |
| Provenance | `SolveUniverseKey` and opt-in solver provenance exist in experimental reuse. | The solve-key digest, producer policy, units, and presence of solver provenance are carried without creating Results records. |
| Completeness | `balls.complete`, reuse certificates, and opt-in provenance exist. Production vertex/topology provenance deliberately says it is not a proof. | Completeness is `PARTIAL` or `UNRESOLVED` until all relevant proof exists; it never supplies zero for an unavailable primitive. |

## Status policy

`CERTIFIED` on an envelope primitive means only that this adapter has a
validated exact geometric representation of that primitive. It does **not**
certify a curvature or topology measurement.

- `PARTIAL` means some exact records are usable but global cells, stable input
  identities, boundary incidence, or solver completeness are incomplete.
- `UNRESOLVED` means the producer has a row but cannot preserve an exact
  representation (for example, an unsupported edge topology); `analytic` is
  `None` rather than zero or a sampled approximation.
- `NOT_CALCULATED` is used for an unbuilt network with no topology primitives.
- `NOT_SUPPORTED` is used specifically for the inapplicable spherical-patch
  interpretation of an AW pairwise face.
- `STALE` is not assigned by FORGE. It must be assigned by the cache/result
  owner when an otherwise valid envelope is reused outside its solve identity.

## Reuse audit

The experimental reuse path remains disabled by default. It already has
deterministic solve-universe keys, immutable extraction/materialization,
direct-solver comparison, local component checks, environment rejection, and
an explicit `EXACT_REUSE` decision only after a matching direct interface.

Its production provenance intentionally sets vertex, edge, and surface
completeness to unproven. Therefore a production envelope is not a reuse
certificate. A direct comparison can still demonstrate geometric/topological
equivalence, but a production enablement proposal requires solver-instrumented
completeness proof on the selected universe plus retained direct-solve
regressions.

## Validation evidence (2026-10-09)

- Synthetic envelope regressions cover deterministic IDs, atom/cell/surface/
  edge/vertex incidence, two-sided ownership, exact generator reconstruction,
  and a fail-closed unresolved edge.
- The existing AW reuse suite passes its 24 synthetic and bounded-production
  tests. It confirms immutable extraction, direct-solver mismatch reporting,
  environment rejection, physical-equivalence diagnostics, and disabled
  production completeness certification.
- The focused envelope, reuse, analytic surface, and analytic edge suites pass
  together (106 tests).
- The real cached EDTA interface regression is intentionally skipped when its
  two documented historical export artifacts are not present. The artifacts
  were not regenerated during this audit; therefore this change has no new
  real-EDTA equivalence claim.
- A standalone bounded timing probe stopped at surface construction because
  this runtime promotes Numba's reflected-list deprecation warning to an
  exception outside pytest. The same synthetic production path passes through
  pytest's repository warning policy; this is an environment/tooling item to
  resolve before using standalone timing scripts as acceptance evidence.

## GAUSS/ATLAS hand-off

GAUSS can consume, per component, the exact quadratic surface coefficients,
exact edge parameters, junction coordinates, owner-normal convention, and
closed/open boundary diagnostics. ATLAS should retain ownership of the typed
quantities and statuses for:

- area;
- smooth, edge, boundary, and corner portions of integrated mean curvature;
- smooth, boundary, and corner portions of integrated Gaussian curvature;
- Euler characteristic and genus.

No Results contract was modified by this producer change.

## Bounded performance plan

1. Keep reuse opt-in and first profile the existing timings on synthetic and
   cached EDTA fixtures. The current timers already separate vertex search,
   topology, edges, surfaces, curvature, and analysis.

   A seven-generator synthetic full solve measured 3.277 s with 2.700 s
   (82%) in surface construction; its recorded geometric predicates were
   0.671 s and mesh generation 0.333 s. This is one cold bounded observation,
   not a benchmark. Prioritize surface predicate/mesh profiling before any
   replacement kernel. The corresponding one-face interface solve took
   0.429 s, including 0.416 s in surface construction.
2. Use immutable extraction/materialization only when direct equivalence and
   completeness evidence are present; avoid rebuilding mesh/perimeter data for
   identical solve universes.
3. Profile DataFrame row iteration and repeated perimeter/edge resolution
   before changing kernels. Cache exact resolved edges only with solve-key
   provenance and immutable inputs.
4. Parallelize only disconnected components or independently validated
   surfaces after deterministic ordering and regression comparisons are in
   place.
5. Consider compiled kernels only after a bounded profile shows a material
   bottleneck and geometry/topology regression tests protect the result.
