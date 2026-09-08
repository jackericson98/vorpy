# VorPy Network-Build Performance and Parallelization Audit

**Status:** investigation and planning only (2026-09-08). No production network-building code was modified during this audit. Existing user changes and untracked files were preserved.

## Executive summary

The current implementation already has meaningful staged timers and several recent optimizations: spatial-query caching, squared-distance ordering, cached verification arrays, Numba Power/Primitive and AW geometry kernels, and vectorized surrounding-array construction. The evidence supports further work, but not a blind multiprocessing or native-extension rewrite.

The five most important measured bottlenecks are:

1. Power/Primitive candidate traversal and Python object/list filtering. In the Cambrin selected-group/full-context run this consumed 0.680 s (Power) and 0.595 s (Primitive), across approximately 436k/431k candidate visits.
2. Power/Primitive surrounding verification preparation. It was 1.590 s and 1.560 s respectively after the latest vectorization, about 25–27% of vertex-search wall time.
3. AW candidate filtering plus geometry. DB1976 recorded 2.598 s in the combined AW filtering bucket and 1.790 s in AW geometry; the buckets overlap by design, so they must not be added.
4. Topology construction (`get_build_edges`, `get_build_surfs`) is staged but not yet represented in the same end-to-end benchmark table as vertex search. It contains Python dictionaries, lists, set operations, deterministic sorting, and duplicate-preserving AW/doublet rules.
5. Surface discretization and analysis are serially orchestrated per surface/cell and contain Python/Shapely/object-array work. Their current self-time share is not adequately measured for the complete CLI run.

The best initial parallelization boundary is **across independent trajectory frames**, using process-level workers, isolated output directories, bounded worker counts, and deterministic per-frame manifests. Within one network, the first useful boundary is independent surface construction or analysis tasks only after immutable array inputs and deterministic indexed merging are introduced.

No C/C++/Cython/Rust extension is justified yet. Existing Numba kernels are already effective for stable numeric geometry, and the remaining measured hotspots are primarily object traversal, data layout, topology construction, and unmeasured surface/analysis stages. A native extension should be reconsidered only after a complete profile shows a stable contiguous numeric kernel dominating runtime after algorithmic and data-layout work.

A realistic near-term expectation is **10–25% end-to-end improvement for a single network**, not an order-of-magnitude gain. Trajectory-level process parallelism can approach the number of physical cores when frames are independent and memory/I/O do not dominate.

## Preservation inventory and worktree

Current uncommitted files at audit start:

- `VorPy_Workbench_Functionality_Blueprint.md` — preserved, not modified.
- `pow_verts.txt` and `prm_verts.txt` — preserved, not modified.

Recent performance commits inspected include:

- `3076d542` edge spatial-query cache;
- `8dea48c5` squared-distance candidate ordering;
- `2fc2df35` surrounding verification-array cache;
- `9815e125` repeated four-ball geometry cache;
- `4b9f6348` Numba flat-vertex solver;
- `d88e5e82` cached AW verification arrays;
- `008f47c2` compiled AW vertex solver;
- `4afe9119` indexed AW solver inputs;
- `3c7411c5` vectorized surrounding-array construction;
- `ed282b9a` AW timing instrumentation.

The audit followed the repository instructions in `AGENTS.md`. No profiling package was installed and no medium/large benchmark was started during this audit.

## Current pipeline and call flow

```text
CLI argv
  -> command.interpret / legacy argv resolution
  -> Command._run_pipeline
      -> System(file=...) / input parser and molecular classifications
      -> parse_commands, load_files, apply_settings
      -> command.group / ggroup -> Group selections
      -> Group.build -> Network.build
          -> System spatial-index reuse or Network.sort_balls
          -> find_net_verts
              -> cached vertex state load
              -> find_verts / find_v0 / fast or slow site search
              -> edge candidate traversal, verification, reseeding
              -> encapsulation and doublet expansion
              -> vertex DataFrame and vertex export
          -> Network.connect -> build_net.build
              -> ball->vertex index
              -> doublet edges
              -> regular edges / triple pairing
              -> edge filtering and adjacency
              -> surface candidate maps / closure validation
              -> surface adjacency and packaging
          -> Network.build_edges -> curved/straight edge geometry
          -> Network.build_surfaces -> per-surface mesh/triangulation
              -> perimeter, COM, projection, triangulation
              -> hyperboloid projection or flat unprojection
              -> curvature and surface descriptors
          -> Network.analyze
              -> surface/cell gathering and completeness
              -> geometry, neighbors, contacts, COM/MOI, curvature summaries
          -> completion summary and metrics
      -> Command.run_exports / output presets
          -> info/logs/vertices/edges/surfaces/atoms/mesh/etc.
```

## Stage map

| Stage | Current owner | Inputs/outputs and mutation | Dominant work / complexity | Current timing | Parallel assessment |
|---|---|---|---|---|---|
| CLI/settings | `command/vpy_cmnd.py`, `command/interpret.py`, `command/set.py` | argv -> settings; mutates `Command` and `System` settings | parsing and validation, O(argv) | not separately timed | inherently serial, negligible |
| Input parsing | `system/system.py`, `inputs/*` | files -> `System`, balls, residues/chains; mutates system | parsing, classification, radii; O(atoms) | not consistently timed | serial per file; frames independent |
| Ball construction | `System`, `Group` | atom records -> DataFrames/lists | object allocation and pandas setup | partially in setup | serial setup; vectorizable data conversion |
| Spatial indexing | `Network.sort_balls`, `calculations/sorting.py` | locations/radii -> box dict and globals; mutates network/system cache | O(N) assignment, Python dict/list writes | progress only; not a reconciled stage timer | parallel writes require partition/merge; low priority |
| Seed/vertex search | `find_net_verts`, `find_verts`, `find_v0`, `fast.py`, `slow.py` | arrays/lists/topology state -> vertex lists and ball adjacency; heavily mutates shared state | candidate traversal, spatial queries, four-ball geometry, validation; approximately edge × candidates | detailed `net.vert_timing` and fast counters | edge traversal is conditionally parallel only after isolated state/merge refactor; current DFS/order is serial |
| Encapsulation/reseed | `find_net_verts` | mutable sphere-check list and vertex state | repeated spatial checks and reseeding | timed inside vertex stage | conditionally parallel by disconnected seed, but shared topology and deterministic ordering risk |
| Doublets | `find_net_verts`, `mark_doublets`, `build_net.doublify` | vertex lists -> duplicated vertices/edges | list insertion and candidate scans | timed in vertex/topology stages | small and ordering-sensitive; leave serial |
| Edge construction | `build_net.get_build_edges`, `Network.build_edges`, `build_edge1.py` | vertex topology -> edge records/geometry; mutates DataFrames | triple indexing, pair matching, curved edge geometry | topology timing exists; geometry timing is incomplete | regular triples can be parallelized after immutable index; duplicate AW/doublet merge must be deterministic |
| Edge adjacency | `build_net.add_build_edges` | edges -> ball/vertex adjacency lists | O(edges) list appends | timed as edge adjacency | parallel local accumulation then deterministic concatenate is possible |
| Surface topology | `build_net.get_build_surfs`, `add_build_surfs` | edge/vertex maps -> surfaces and reverse adjacency | dict maps, combinations, closure checks; O(E+V) plus candidate combinations | timed as surfaces/surface adjacency | surface keys independent after maps; deterministic sorted merge required |
| Surface geometry | `build_surf.py`, `perimeter.py`, `fill.py`, `triangulate.py` | loc/rad/edge data -> points, triangles, curvature | per-surface Python geometry, Shapely, triangulation, projection | detailed local timing exists; not reconciled into full CLI table | strongest within-network parallel candidate after isolated surface records |
| Analysis | `analyze.py`, calculations caches | network DataFrames -> cell/surface metrics and summaries; mutates DataFrames | per-ball loops, adjacency, contacts, COM/MOI, curvature | detailed `ANALYSIS TIMING` exists | per-ball tasks can be isolated after immutable array preparation; reductions need deterministic merge |
| Export | `output/*`, `group/export.py`, `output/output.py` | network/group -> files; filesystem side effects | serialization, mesh preparation, text and binary I/O | `export_timing` exists; often mixed with CLI progress | independent artifacts can be process tasks, but I/O bandwidth limits benefit |

## Timing boundaries and gaps

The intended categories are:

```text
Input/setup -> spatial indexing -> network construction -> analysis -> geometry preparation -> export
```

Current implementation has useful but separate timing trees:

- vertex search: `net.vert_timing`, including seed discovery, site-container search, reseeding, candidate validation, and vertex export;
- topology: `build_net` timing buckets for index creation, doublets, regular edges, filtering, adjacency, surfaces, validation, packaging;
- surface geometry: `build_surf` and `triangulate` timing dictionaries;
- analysis: `analyze` timing buckets plus explicit “other/loop overhead” calculation;
- export: `ExportProgress.timings` and `sys.export_timing`.

Gaps that prevent complete reconciliation:

1. input parsing and molecular classification are not included in the same hierarchy;
2. `Network.sort_balls`/spatial indexing has progress but no parent-child timing bucket;
3. `Network.build_edges` geometry is not charged into the topology timing tree consistently;
4. `Network.build_surfaces` orchestration and all surface timing are not always attached to `net.metrics`;
5. `Group.get_info`, layers, and summary/export work can occur outside `Network.metrics['tot']` depending on caller;
6. CLI total, analysis/build completion, and export progress are not currently reconciled into one trace.

The first implementation task should be **timing-only instrumentation at these boundaries**, not an algorithm rewrite. Every stage should have `start/end`, child buckets, call counts, and an explicit unattributed remainder.

## Evidence collected so far

Environment observed in the audit shell:

- Linux x86_64, kernel 6.14.0-36-generic;
- Python 3.14.6;
- NumPy 2.5.2;
- SciPy 1.18.1;
- Numba 0.67.0;
- pandas 3.0.5.

Existing measured data (all 5 Å maximum-vertex runs unless noted) is in `network_build_baselines.csv` and `network_build_hotspots.csv`.

Key results:

- DB1976 vertex search: AW 722 vertices / 5.85 s; Power 731 / 3.01 s; Primitive 703 / 2.72 s.
- Cambrin first ten residue records with all 9,190 atoms retained as context: Power 1,730 vertices / 6.35 s; Primitive 1,708 / 5.87 s.
- User-provided full CLI Cambrin run: 674 complete cells, 7,676 vertices, 7,496 surfaces; network completion around 90.25 s and export progress around 146.55 s. This is not stage-reconciled and should not be mixed with vertex-only timings.
- A selected-group AW benchmark failed to find a seed after a long search. This is a valid diagnostic of the current seed path, not an AW throughput measurement.

The Cambrin Power/Primitive metrics show that vectorized surrounding-array construction reduced surrounding preparation from approximately 2.7/2.6 s to 1.59/1.56 s in comparable selected-group runs. This is evidence for data-layout work, not evidence that all network stages scale similarly.

## Proposed benchmark plan (not yet run)

The user requested staged fixtures and review before expensive runs. Proposed commands:

### Small: EDTA

```bash
PYTHONPATH=. .venv/bin/python scripts/benchmark_vertex_search.py vorpy/data/EDTA.pdb --residues 0 --network pow
PYTHONPATH=. .venv/bin/python scripts/benchmark_vertex_search.py vorpy/data/EDTA.pdb --residues 0 --network prm
```

Use three repetitions after one warm-up for each network type. Add a separate complete `Network.build()` harness after timing boundaries are unified.

### Medium: Cambrin

```bash
PYTHONPATH=. .venv/bin/python scripts/benchmark_vertex_search.py vorpy/data/cambrin.pdb --residues 10 --full-context --network pow
PYTHONPATH=. .venv/bin/python scripts/benchmark_vertex_search.py vorpy/data/cambrin.pdb --residues 10 --full-context --network prm
```

The AW whole-system command is currently not a bounded useful fixture because it remains in seed discovery. Do not rerun it until the seed-path mismatch with the successful CLI workflow is explained.

### Protein: p53tet

Proposed only after small/medium timing review:

```bash
vorpy p53tet -s mv 5
```

Capture build, analysis, and export separately; do not infer vertex-search speed from the CLI total.

### Large: NCP

Defer until the timing-only harness and p53tet profile identify a specific question:

```bash
vorpy NCP -s mv 5
```

For every fixture, record atom/ball count, max vertex, scheme, vertices, edges, surfaces, complete cells, wall/CPU time, peak RSS, stage timings, call counts, output preset, and software/CPU metadata. Use `time.perf_counter()` for wall time and `time.process_time()` for CPU time. Peak memory can use `/usr/bin/time -v` on Linux without adding a Python dependency.

## Profiling plan

1. Use existing `net.vert_timing`, `build_net` timings, `build_surf` timings, `analyze` timings, and `export_timing` first.
2. Add temporary boundary instrumentation only in a separate audit branch or local patch; do not leave noisy prints in production.
3. Run `python -m cProfile -o /tmp/vorpy-small.prof ...` only on EDTA after the benchmark command is reviewed.
4. Inspect cumulative and self time with `python -m pstats` or a small standard-library reader.
5. Use `tracemalloc` only around setup/analysis/export if allocation evidence is needed; it does not represent native allocations perfectly.
6. Do not install `py-spy`, Scalene, line-profiler, or other profiling packages without approval.

## Algorithmic findings

### Vertex search

The current spatial grid reduces expected candidate scope, but the search remains a stateful depth-first traversal. Candidate filtering repeatedly traverses Python lists, sets, DataFrame-derived object arrays, and adjacency lists. The highest-value algorithmic opportunity is a structure-of-arrays candidate representation and batched candidate validation, not simply more workers.

The current cache work removes repeated spatial and geometry work but has limited reuse on unique-edge fixtures. Cache hit rates should be measured before expanding cache scope.

### Topology

`get_build_edges` builds a triple index and then performs deterministic pairing, including exact matching for AW duplicate multiplicities. `get_build_surfs` builds edge and vertex maps and performs closure checks. These are good candidates for immutable numeric arrays and deterministic map/reduce, but duplicate-preserving doublet semantics and stable ordering are correctness constraints.

### Surface construction

Each surface is largely independent after topology is fixed. The work includes perimeter generation, projection, triangulation, Shapely containment, curved projection, and curvature. The algorithmic question is whether perimeter/triangulation work can reuse surface-level geometry and avoid repeated Python conversions. Parallelism should come only after each surface returns an isolated record keyed by original surface ID.

### Analysis

`analyze` has already removed much pandas overhead by caching arrays and preallocating outputs. Remaining work is mixed: Python adjacency traversal and completeness checks, numeric kernels, contact accumulation, and COM/MOI. Per-ball parallelism is plausible after read-only array preparation, but surface contact reductions require deterministic indexed accumulation.

### Export

Export is serial filesystem work and can dominate CLI wall time. It must be excluded from network-build comparisons. Independent artifacts can be scheduled separately, but parallel export may saturate disk bandwidth and increase memory pressure.

## Parallelization dependency graph

| Region | Unit | Shared state | Merge | Recommendation |
|---|---|---|---|---|
| Trajectory frames | frame | none if input/cache immutable | manifest/status records | Best initial process boundary; isolated output dirs |
| Groups/interfaces | group/interface | system spatial index read-only; each Network mutable | group result records | Process-level only after memory measurement and output isolation |
| Vertex edges | edge/candidate | current DFS vertex/adjacency state | deterministic vertex/adjacency merge | Not safe in current architecture; refactor first |
| Edge triples | triple key | read-only vertex arrays | stable sorted edge list | Good conditional parallel target after topology SoA refactor |
| Surfaces | surface ID | read-only topology and geometry arrays | indexed surface records | Best within-network target after stage timing |
| Analysis | ball/cell | read-only surfaces; per-ball outputs | indexed arrays/reductions | Conditional parallel target; preserve completeness/contact semantics |
| Export artifacts | file/artifact | filesystem and shared progress | manifest | Coarse process tasks, only when I/O permits |
| Numeric kernels | contiguous arrays | none | array outputs/reductions | Numba/vectorization first; threads only if kernel releases GIL |

Python threads should not be assumed to help CPU-bound loops. Current Python orchestration and object traversal retain the GIL. Numba kernels can release the GIL during compiled execution, but the current decorators do not establish `parallel=True`/`prange` semantics. Multiprocessing introduces serialization and duplicated memory unless arrays are moved to shared memory.

## Acceleration decision table

See `network_build_opportunities.csv` for machine-readable entries. The short decision is:

- **Algorithm/data layout:** recommended first for candidate traversal and topology.
- **NumPy/Numba:** recommended for contiguous geometry, verification, curvature, and indexed reductions after profiling.
- **Multiprocessing:** recommended first across trajectory frames; conditional for surfaces/analysis.
- **C/C++/Cython/Rust:** not currently justified. The remaining hotspots are not yet proven to be stable numeric kernels that dominate end-to-end time after data-layout work.
- **Leave unchanged:** doublet matching and correctness-sensitive ordering until profiles show a meaningful share.

## Compiled-language criteria

A future extension candidate must have:

- measured dominant self time;
- contiguous numeric inputs and outputs;
- stable mathematical semantics;
- many operations per Python/native boundary crossing;
- no Python object graph traversal;
- deterministic output;
- independent unit tests and a pure-Python/Numba fallback; and
- reliable Windows/Linux/macOS packaging.

Likely future candidates, in order, are batched triangle curvature or surface projection only if Numba cannot handle the measured self time. Current object-heavy vertex/topology stages are poor extension boundaries until converted to structure-of-arrays inputs.

## Correctness and regression gates

Every optimization must compare, at minimum:

- vertex counts, coordinates, radii, and doublet solutions;
- edge endpoint and defining-ball sets;
- surface ball/edge/vertex sets and closure;
- complete/incomplete cells;
- volumes, surface areas, contacts, and overlap values;
- mean/Gaussian/integrated curvature;
- representative surface energy;
- group/interface membership;
- deterministic ordering where externally visible; and
- normalized output-file contents.

The existing ALA78 CA/C multi-edge regression is in `vorpy/tests/network/test_build_net.py` and must remain in every topology optimization gate. Existing network tests currently cover 15+ cases after the recent performance work; the project-wide suite still has environment collection blockers (`hypothesis` and `mesh`) unrelated to this audit.

## Amdahl-law estimates

These are planning ranges, not measured promises:

| Proposal | Assumed fraction of single-network build | Kernel speedup | Estimated end-to-end effect |
|---|---:|---:|---:|
| Candidate SoA/batched filtering | 20–35% | 2–4x | 10–20% |
| Surrounding-array reuse/indexing | 20–27% of selected vertex stage | 1.5–2x | 5–12% |
| Surface per-ID parallel map | 20–40% of full build after measurement | 2–6x | 10–25% |
| Analysis per-ball parallel map | 10–30% of full run after measurement | 2–4x | 5–15% |
| Independent trajectory frames | 80–95% of multi-frame workload | near core count before I/O/memory limits | 2x: ~1.7–1.9x; 4x: ~3.0–3.6x; 8x: ~5–6.5x |

For a single network, memory bandwidth, topology merges, serial setup, and export will cap speedup. For frames, process startup and output contention determine break-even; a minimum of several medium frames is likely needed before parallel scheduling pays off.

## Ranked implementation roadmap

### Tier 1 — low-risk, high-confidence

1. **Unify timing boundaries and counts.** Add temporary or agreed production timing for parsing, sorting, edge geometry, surface orchestration, analysis, and export. Reconcile parent/child totals and expose “other.” Tests: timing-schema tests and no topology change. Fixture: EDTA then Cambrin. Rollback: remove instrumentation only.
2. **Complete Power/Primitive candidate data-layout pass.** Measure cache hit rates, convert candidate/topology state to indexed numeric arrays, and batch safe filtering. Tests: ALA78 topology plus coordinate/metric equivalence. Expected: 8–15% end-to-end. Prerequisite: three warm repetitions.
3. **Reuse immutable surface/verification arrays.** Continue vectorized indexing and distinguish cache hits from misses. Expected: 5–10% in vertex stage where reuse exists.

### Tier 2 — parallel architecture

4. **Trajectory-frame process runner.** Define immutable input/result manifests, isolated output directories, bounded workers, cancellation, and deterministic aggregation. Tests: frame ordering, failure recovery, output collision checks. Expected: near-linear coarse-grained scaling.
5. **Per-surface task interface.** Return `{surface_id, points, tris, descriptors, timing, status}` and merge by ID. Preserve surface order and closure. Only implement after full surface timing.

### Tier 3 — compiled numerical kernels

6. **Batch curvature/projection kernel only if profiling justifies it.** Define contiguous array API, status/error codes, fallback, and tolerance tests. Do not start a C/C++/Rust extension before Numba self-time evidence.

### Tier 4 — structural redesign

7. **Network structure-of-arrays model.** Replace DataFrame/object-list hot-path traversal with immutable numeric arrays plus explicit adjacency CSR-like offsets. This is high effort and should happen only if Tier 1/2 profiles show object traversal dominates.

## Recommended first implementation task after approval

Implement the timing-only audit harness and run EDTA with three warm repetitions for one network scheme. It should produce one JSON/CSV record with:

- input/setup;
- spatial indexing;
- vertex search;
- topology/connect;
- edge geometry;
- surface geometry;
- analysis;
- export;
- total and unattributed time;
- wall and CPU time;
- peak RSS;
- counts and environment metadata.

Do not change numerical code in that task. Once reviewed, run Cambrin Power/Primitive, then p53tet. Defer NCP until the profile identifies a concrete question.

## Unresolved questions

1. Which CLI settings correspond exactly to the successful `vorpy cambrin -s mv 5` run, including network scheme and build/export preset?
2. Why the direct AW benchmark enters the slow seed path while the CLI completes; this must be resolved before using AW benchmark numbers.
3. Which `Network.build_surfaces` and export paths are included in the user’s reported 90/146-second CLI timestamps?
4. What physical core count, RAM, and storage characteristics apply on the Ubuntu benchmark machine?
5. Are trajectory frames independent in the user’s production workflow, and can output directories be isolated per frame?
6. Which exact EDTA/p53tet/NCP commands and output presets are authoritative for comparison?

## Links

- [Baseline measurements](network_build_baselines.csv)
- [Measured hotspots](network_build_hotspots.csv)
- [Opportunities and decisions](network_build_opportunities.csv)
