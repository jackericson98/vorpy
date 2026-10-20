# NEXUSA interface/export acceptance contract

Published 2026-10-19 by the coordinator for FORGE, PRISM, VECTOR, GAUSS and
ATLAS review. This is the shared seam for the current integration cycle.

## Ownership

- GAUSS defines scientific metric meaning, units, scope, completeness and
  certification. Missing or partial quantities remain explicit, never zero.
- ATLAS owns typed `Results`, provenance, source/solve/frame/selection identity,
  dependency validation and archive round trips.
- FORGE owns network geometry, physical interface extraction, FULL_SOLVENT
  competitor membership and completeness evidence.
- VECTOR owns export serialization, scoped meshes/tables, canonical logs and
  archive I/O. Export consumes cached geometry/Results and never calls a solve.
- PRISM owns CLI routing, progress, failure boundaries and user-visible timing.

## Acceptance gates

1. `-i` builds only the requested interface; `-i -c system` additionally builds
   one full-system source and derives the interface without a second geometry
   solve.
2. Ordered atom selections and source identity are preserved. The wet 1BNA
   input hash and group membership are recorded before any pilot.
3. Direct contact area and total network area are separate named measurements;
   curvature/topology and alpha-derived values carry honest status and units.
4. FULL_SOLVENT preserves competing solvent generators and does not weaken
   completeness rules. Partial/unresolved states remain exportable but are not
   promoted to certified Results.
5. The default `interface` preset exports physical interface OFF geometry,
   contact atom PDBs, scientific summary/contact tables, completeness fields and
   canonical typed-Results logs only. Archives, legacy logs, dual/visualization
   products and full-system exports remain explicit.
6. Geometry, analysis, export and optional archive-save timings are separately
   observable. Archive-reloaded export is cache-only and must not rebuild.
7. No solver or scientific formula changes are accepted as export optimizations.

The deferred DNA-side/ligand-side reconstructed dual-contact patches remain
documentation-only and are not an acceptance gate for this cycle.
