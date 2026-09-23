# VorPy performance handoff

## Repository checkpoint

- Repository: `/home/jack/PycharmProjects/vorpy`
- Branch/commit: `master` at `a2421e1c50882c40d24801d84d226ee26171bc8a`
  (`Add network performance audit roadmap`)
- Repository instructions: no `AGENTS.md` exists inside this repository. The
  sibling Workbench has its own `AGENTS.md`; it does not govern VorPy code.
- Uncommitted files (user-owned; preserved):
  `VorPy_Workbench_Functionality_Blueprint.md`, `pow_verts.txt`, and
  `prm_verts.txt`.

This document is a performance-only handoff. It does not authorize production
algorithm changes. The committed report and CSVs below are the detailed audit:

- [audit report](../performance/NETWORK_BUILD_PERFORMANCE_AUDIT.md)
- [baselines](../performance/network_build_baselines.csv)
- [hotspots](../performance/network_build_hotspots.csv)
- [opportunities](../performance/network_build_opportunities.csv)

## Current architecture and verified execution path

The normal CLI enters `vorpy/__main__.py`, then
`vorpy/src/command/vpy_cmnd.py:Command._run_pipeline`. Input parsing creates a
`System`; command settings select groups and schemes; `Group.build()` creates a
mutable `Network`. The current build path is:

```text
CLI/settings -> input parsing/classification -> ball arrays/DataFrames
-> spatial indexing/sorting -> find_net_verts/find_verts/find_v0
-> doublets/encapsulation -> connect/build_net topology
-> build_edges -> build_surfaces/triangulation
-> analyze (cells, contacts, geometry, curvature, energy)
-> summaries/logs -> export/output artifacts
```

Important modules and entry points:

- `vorpy/__main__.py`, `vorpy/src/command/vpy_cmnd.py`: CLI and pipeline.
- `vorpy/src/system/system.py`, `vorpy/src/inputs/`: structure and ball setup.
- `vorpy/src/network/network.py`: `Network.build`, stage orchestration.
- `vorpy/src/network/find_net_verts.py`, `find_verts.py`, `find_v0.py`,
  `fast.py`, `slow.py`: vertex search and validation.
- `vorpy/src/network/build_net.py`: edges, adjacency, surfaces, doublets.
- `vorpy/src/network/build_surf.py`, `triangulate.py`, `perimeter.py`, `fill.py`:
  surface geometry.
- `vorpy/src/network/analyze.py` and `vorpy/src/calculations/`: analysis and
  curvature/energy calculations.
- `vorpy/src/output/`, `vorpy/src/group/export.py`: file export.
- `scripts/benchmark_vertex_search.py`: bounded vertex-search benchmark.

The implementation already contains spatial-query/cache work, squared-distance
ordering, cached verification arrays, Numba geometry kernels, indexed AW input,
and vectorized surrounding-array construction. Preserve those changes.

## Completed performance work and measured evidence

Recent performance commits include `3076d542`, `8dea48c5`, `2fc2df35`,
`9815e125`, `4b9f6348`, `d88e5e82`, `008f47c2`, `4afe9119`, `60189924`,
`b9c61c45`, `80ad6c50`, `b61d0f51`, `35a25b11`, and `3c7411c5`.

Measured 5-A maximum-vertex data currently available:

| Fixture/scope | Scheme | Result | Wall time | Main measured cost |
|---|---:|---:|---:|---|
| DB1976 vertex search, 41 balls | AW | 722 vertices | 5.85 s | site-container 5.57 s |
| DB1976 vertex search, 41 balls | Power | 731 | 3.01 s | site-container 2.81 s |
| DB1976 vertex search, 41 balls | Primitive | 703 | 2.72 s | site-container 2.54 s |
| Cambrin, first 10 residues, all 9,190 atoms as context | Power | 1,730 | 6.35 s | surrounding 1.59 s; candidate filtering 0.68 s |
| Same Cambrin scope | Primitive | 1,708 | 5.87 s | surrounding 1.56 s; candidate filtering 0.59 s |

The user’s full Cambrin CLI run reported 674 complete cells, 7,676 vertices,
7,496 surfaces, network completion around 90.25 s, and export progress around
146.55 s. It is not stage-reconciled and must not be used as a vertex-kernel
benchmark.

AW selected-group/full-context benchmarking currently stalls in seed discovery,
while the normal CLI succeeds. Treat that as an unresolved benchmark-path
mismatch, not an AW throughput result. Do not launch whole-system AW profiling
without a bounded timeout and an explanation of this difference.

## Timing/profiling status

Existing timers are distributed across `net.vert_timing`, `build_net` timing
buckets, surface/triangulation timing, `analyze` timing, and export timing.
They are useful but do not yet reconcile one end-to-end hierarchy. Missing or
ambiguous boundaries include input/classification, spatial indexing, edge
geometry, surface orchestration, group/info generation, and the relationship
between network completion and export.

The first approved implementation should be timing-only instrumentation:

1. `Input/setup`
2. `Spatial indexing`
3. `Network construction`
4. `Network analysis`
5. `Geometry preparation`
6. `File export`
7. end-to-end total

Every parent must reconcile child timings and report unattributed remainder.
Use `perf_counter` and `process_time`; use `/usr/bin/time -v` for Linux peak
RSS. Use `cProfile`/`pstats` and `tracemalloc` only after small fixtures are
approved. Do not install `py-spy`, Scalene, line-profiler, or other packages
without user approval.

## Bottlenecks and candidate improvements

The five best-supported hotspots are candidate traversal/object filtering,
surrounding verification preparation, AW filtering/geometry, unmeasured but
likely Python-heavy topology construction, and surface/analysis geometry that
is not yet reconciled in the complete CLI profile.

Priority roadmap:

### Tier 1 — low-risk, high-confidence

- Complete hierarchical timers and three-repetition small-fixture baselines.
- Measure cache hit/miss rates before extending caches.
- Convert candidate state toward immutable structure-of-arrays views and reuse
  contiguous index arrays; retain current Python fallback.
- Batch/vectorize candidate filtering and verification only where numerical
  equivalence is demonstrated.
- Add complete build/analysis timing to the existing benchmark harness.

### Tier 2 — parallel architecture

- Best first boundary: independent trajectory frames, each in a bounded process
  with isolated output directory and deterministic per-frame manifest.
- Later conditional boundaries: independent surfaces, then isolated analysis
  records. Both require immutable inputs, indexed results, deterministic merge,
  cancellation, and memory measurements.
- Groups/interfaces may use processes only after proving system caches are
  read-only and output paths cannot collide.
- Do not parallelize current vertex DFS, doublet handling, or AW ordering yet:
  they mutate shared adjacency/topology state and ordering affects results.

### Tier 3 — compiled numerical kernels

Existing Numba kernels should be measured and extended first. A C/C++/Cython/
Rust extension is **not currently justified**: the remaining measured costs are
mostly Python object traversal, data layout, topology, and unmeasured surface
orchestration. Reconsider only for a stable contiguous numeric kernel that is a
large self-time fraction after Tier 1. Candidate API shape would be contiguous
float64 coordinate/radius arrays plus int32 connectivity/index arrays returning
fixed-shape numeric arrays and status flags, with a pure-Python/Numba fallback.

### Tier 4 — structural redesign

Only if profiling proves necessary: replace object-heavy vertex/topology state
with structure-of-arrays records, then add deterministic map/reduce boundaries.
Do not begin with a wholesale rewrite.

## Correctness invariants and tests

Every optimization must preserve vertex counts/coordinates, edge endpoint and
defining-ball sets, surface closure and memberships, complete/incomplete cells,
volumes, surface areas, contacts, curvature, representative surface energy,
group/interface membership, deterministic externally visible ordering, and
normalized output-file semantics. Keep the ALA78 CA/C multi-edge regression in
`vorpy/tests/network/test_build_net.py` in every topology gate.

Current validation commands:

```bash
PYTHONPATH=. .venv/bin/python -m pytest -q
PYTHONPATH=. .venv/bin/python -m compileall -q vorpy scripts
```

Known environment caveat from the audit: project-wide collection can be blocked
by missing `hypothesis` and `mesh` modules; report that separately from test
failures. Existing network tests and the performance benchmark are the focused
validation surface.

Benchmark fixtures and bounded commands:

```bash
PYTHONPATH=. .venv/bin/python scripts/benchmark_vertex_search.py vorpy/data/EDTA.pdb --residues 0 --network pow
PYTHONPATH=. .venv/bin/python scripts/benchmark_vertex_search.py vorpy/data/EDTA.pdb --residues 0 --network prm
PYTHONPATH=. .venv/bin/python scripts/benchmark_vertex_search.py vorpy/data/cambrin.pdb --residues 10 --full-context --network pow
PYTHONPATH=. .venv/bin/python scripts/benchmark_vertex_search.py vorpy/data/cambrin.pdb --residues 10 --full-context --network prm
vorpy p53tet -s mv 5
```

Run one warm-up then at least three EDTA repetitions; use fewer repetitions for
medium/protein fixtures based on cost. Defer NCP until the timing harness answers
a specific question.

## Incomplete work and uncertainties

- No complete reconciled timing trace exists.
- AW benchmark seed behavior differs between the benchmark path and CLI.
- Surface construction and analysis self-time need a complete profile.
- Cache reuse and memory scaling are not quantified.
- Threading benefits are unproven; Python orchestration generally holds the GIL.
- Export can dominate CLI time and must remain separate from build claims.

## Explicit do-not-do constraints

- Do not change numerical tolerances, geometric behavior, validation checks, or
  output ordering for speed without a reviewed invariant plan.
- Do not introduce multiprocessing, compiled extensions, or production
  algorithm rewrites before the audit plan is approved.
- Do not profile expensive whole-system AW/NCP runs without a staged review and
  timeout.
- Do not treat total CLI time as network-build time when export/analysis are
  included.
- Do not modify the sibling Workbench repository from this workstream.

## Ready-to-paste opening prompt

```text
Read docs/handoffs/VORPY_PERFORMANCE_HANDOFF.md and
docs/performance/NETWORK_BUILD_PERFORMANCE_AUDIT.md completely. Inspect the
current VorPy worktree and preserve all user changes. Do not modify production
network algorithms yet. First implement only the reviewed timing-boundary
harness for EDTA, with child/parent reconciliation, wall and CPU time, call
counts, and explicit unattributed time. Run one warm-up and three repetitions,
compare against the existing CSV baselines, and report any mismatch. Do not
install profiling packages, run NCP, or parallelize anything without approval.
```
