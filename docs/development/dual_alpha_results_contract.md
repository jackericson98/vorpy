# Dual and alpha Results contract

This is an additive Results contract. It does not change the eight-column
`interfaces/logs.csv` schema or make an abstract dual complex a physical
surface.

## Contract inventory

| Contract | Current populated source | Explicitly unavailable or producer-owned state |
| --- | --- | --- |
| `InterfaceResult` for `voronoi` | `adapt_interface_analysis()` snapshots the physical AW/Power cache, including area, topology, mean curvature, intrinsic Gaussian curvature, alpha-selection counts, coverage, and cache provenance. | Missing source quantities stay null with `NOT_CALCULATED`, `NOT_SUPPORTED`, or `UNRESOLVED`; no default zero is created. |
| `MolecularContactSurfaceResult` | Molecular contact selection, patch provenance, topology, and geometry counts have typed structures. | Molecular geometry measurements depend on its dedicated producer and are not copied from the physical partition. |
| `DualContactComplexResult` | `adapt_dual_contact_complex()` traverses only dimension-one dual incidences. It records generator IDs, primal-feature references, bounded/complete/supported flags, and exact cached cardinalities. | It never reports physical area, mean curvature, Gaussian curvature, or a primal topology measurement. |
| `DualGeneratorComplexResult` | `adapt_dual_generator_complex()` traverses all requested cached simplices into `DualSimplexRecord` and `GeneratorComplex`. | Topology is accepted only when a producer supplies typed quantities; no Euler characteristic or curvature is inferred from counts. |
| `AlphaSelectionResult` | `adapt_alpha_selection()` snapshots alpha identity, selected generator pairs, selected system-pair references, selected dual IDs, selection counts, and coverage. | A selected dual ID is not a physical surface; area/curvature remain absent unless a separate physical producer measures them. |

All result subclasses inherit the common typed fields for identity, group and
selection identity, environment, `Quantity` coverage, topology, mean and
Gaussian curvature, scientific status, provenance, and geometry lineage.
`DualSimplexRecord` is deliberately incidence-only; its flags describe source
support and completeness, not a geometric measure.

## Units, orientation, and status

- Dual cardinalities and alpha counts/coverage use unit `1`. The alpha value
  remains in the explicit native units carried by `ResultIdentity` (for example
  `A^2` for a power convention); it is not converted by these adapters.
- The primal physical interface orientation is still `group_A -> group_B`.
  Dual incidences and alpha-selected pairs have no curvature sign convention.
- `CERTIFIED` dual traversal requires a valid source incidence audit and only
  certified traversed simplex records. An invalid or unavailable audit, or an
  incomplete/unsupported source record, produces `PARTIAL` while preserving
  the source flags.
- `UNRESOLVED`, `NOT_CALCULATED`, and `NOT_SUPPORTED` selection quantities
  contain `null`, never `0`. In particular, the historical unresolved
  alpha-cache `coverage: 0.0` sentinel is converted to an unavailable null.
  `STALE` remains a stale snapshot rather than a fresh certification.

## Producer integration checklist

Before connecting a live producer, the responsible owners confirm these items
against a small synthetic or EDTA case:

1. FORGE supplies a read-only dual object with `simplices[dimension]`, canonical
   `simplex_id` (`"dimension:sorted_generator_ids"`), sorted network-local
   generator IDs, primal-feature references, `bounded`, `complete`,
   `supported`, and source status. Its `audit()` exposes `valid`, generator
   count, per-kind feature counts, and unmapped primal-feature count so coverage
   retains support/total counts. FORGE supplies the explicit
   side/subset IDs; ATLAS never infers side membership from geometry.
2. FORGE supplies cached alpha mappings using the existing vocabulary:
   `selected_generator_pairs`, `selected_system_pairs`,
   `selected_physical_system_pairs`, selection counts, `coverage`, `certified`,
   and `status`. Pair IDs are source identities, not physical measurements.
3. GAUSS supplies any dual topology, mean curvature, or Gaussian curvature as
   already typed `Quantity`, `Topology`, `MeanCurvature`, or
   `GaussianCurvature` blocks with method, scope, coverage, units, sign
   convention, scientific status, provenance, and geometry lineage. It must
   not reuse a primal area or curvature field for a dual record.
4. FORGE and GAUSS preserve the Result identity fields: partition (`aw` or
   `power`), network/radius configuration, groups, interface ID, side,
   orientation, alpha method/value/units, and any selection source/ID. A
   missing geometry lineage remains unresolved rather than being labelled a
   direct solve.
5. VECTOR serializes only supported quantities through the established
   eight-column vocabulary. New dual or alpha rows use the existing canonical
   representations `dual_contact_complex`, `dual_generator_complex`, and
   `alpha_selection`; adding columns or new log quantities needs separate
   authorization.

## Validation boundary

Focused Results tests exercise serialization-ready nested structures,
unavailable alpha values, alpha identity/correspondence, geometry-lineage
provenance, dual traversal, and the unchanged eight-column log header. They do
not run a molecular benchmark or calculate geometry.
