# NEXUSA Core Stabilization Acceptance Contract

Status: shared integration contract, 2026-10-19.  This document freezes the
scope for EDTA baseline restoration.  It does not authorize 1BNA, p53tet, or
large molecular benchmarks.

## Ownership boundaries

| Owner | Owns | Must not silently change |
|---|---|---|
| FORGE | solver geometry, primary Voronoi cells, topology, completeness production, solvent extraction | scientific formulas or primary geometry to satisfy an exporter |
| PRISM | CLI parsing, calculation planning, selection and routing, progress/timing display | solver or metric definitions |
| GAUSS | definitions and acceptance of cells, surfaces, edges, vertices, area, curvature and statuses | CLI/export serialization |
| ATLAS | dual-complex identity, typed Results, provenance and dependency validation | primary geometry or classic log format |
| VECTOR | classic logs, output presets, meshes, serialization and archive I/O | solver geometry, completeness decisions or scientific formulas |
| NEXUSA | integration, bounded tests, acceptance evidence and root `results.txt` | specialist-owned implementation unless coordinated |

Shared interfaces require an additive, backward-compatible field or an
explicit handoff recorded here before implementation.  No specialist changes
`vorpy/src/interface/export.py` except VECTOR; FORGE owns geometry and
completeness production.

## Authoritative CLI behavior for this stabilization cycle

| Invocation | Networks calculated | Solvent participation |
|---|---|---|
| `python vorpy edta` | one ordinary full-system/default group solve | wet inputs use implicit `SOLVENT_CONTEXT`; the molecular group remains the selected cell owner and the loaded structure supplies context |
| `... -c system` | the same ordinary full-system group route when no interface is requested | same as the default; `-i -c system` is the separate source-network case |
| `... -g <selector>` | one independent selected-group solve | selected group semantics only |
| `... -i` with two `-g` selectors | interface network only | interface policy; no implicit independent group export |
| `... -i -c system` with two `-g` selectors | one full-system source solve plus an interface view/reuse | source solve and interface policy are explicit and separately reported |
| `... -s nt aw|pow|prm` | selected scheme only | no cross-scheme substitution |
| `... --solvent-mode MOLECULE_ONLY|SOLVENT_CONTEXT|FULL_SOLVENT` | same requested networks | explicit mode is preserved in settings/provenance; FULL_SOLVENT may fail closed when completeness is unproven |
| `... --profile geometry|interface|analysis|full` | same requested networks | profile controls optional analysis only; it cannot erase primary geometry |

The classic network/group logs remain available through explicit classic log
exports.  Canonical typed Results logs remain a separate representation.  The
default full-system path must not silently route through interface-only code.
Interface exports use the compact physical/scientific preset by default;
`small`, `medium`, `large`, `all`, classic `logs`, and archive saving remain
explicit expansions.  Archive reload must preserve source identity, native
measurements, typed Results and partial/unresolved statuses without rebuilding
primary geometry.

## EDTA baseline gates

The historical reference is the exact documented 32-cell family, not a test
expectation: EDTA PDB, 32 selected molecular atoms, 612 source/context balls,
no boundary, `surf_res=0.2`, `box_size=1.25`, and the recorded network/mv
setting.  The Power timing record is explicitly `-s nt pow -s mv 5` and
reports 504 vertices, 896 edges, 425 surfaces and 32 complete selected cells.
The older AW record reports 500 vertices, 892 edges, 424 surfaces and 32
complete cells with its recorded AW settings.  These are method-specific
references and must not be conflated.

Acceptance requires independent counts for cells, surfaces, edges and
vertices, plus membership and solver settings.  A complete-cell count is not
the count of all constructed cells and cannot by itself certify a solve.
Default and explicit-system EDTA must be equivalent.  Primary geometry must
remain present when dual/typed Results work is unavailable.  Curvature is
either scientifically computed with evidence or explicitly
`NOT_CALCULATED`; no zero is inferred from missing analysis.

## Validation and change discipline

Before/after evidence must include the command, input hash, solvent mode and
active/context membership, per-stage timings, output file counts/bytes, and
equivalence of native geometry/metrics.  Use bounded EDTA and synthetic tests
only.  Do not modify the dirty tracked PDB, overwrite existing outputs, stage,
commit or push.  Deferred reconstructed dual-contact patches are documented
only and are outside this cycle.
