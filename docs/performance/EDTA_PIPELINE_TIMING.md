# EDTA CLI pipeline timing baseline

Measured 2026-09-08 against VorPy `a2421e1c50882c40d24801d84d226ee26171bc8a`.
Production algorithms and Workbench code are unchanged. Instrumentation lives in
`scripts/benchmark_pipeline.py` and `scripts/performance_timing.py`; wrappers are
installed temporarily and restored even when the CLI raises an exception.

The median complete pipeline time was **4.732 s in a fresh process** and
**3.477 s after an additional same-process warm-up**. These are two execution
conditions for the same code, not an optimization speedup. All measured builds
matched their uninstrumented controls exactly for network data and normalized
output files.

## Reproduction and scope

```bash
# Each output path must be new. Each worker has a 120-second outer timeout.
PYTHONPATH=. .venv/bin/python scripts/benchmark_pipeline.py vorpy/data/EDTA.pdb --output /tmp/vorpy-edta-pow-fresh --timeout 120
PYTHONPATH=. .venv/bin/python scripts/benchmark_pipeline.py vorpy/data/EDTA.pdb --output /tmp/vorpy-edta-pow-warm --timeout 120 --warm-kernels
```

Both commands run one uninstrumented control followed by three measured trials.
Each trial uses a separate process and isolated directory. `--warm-kernels`
also runs an uninstrumented CLI build in that worker before its measured build,
checks equivalence to that warm-up, and collects garbage before measurement.
This warms signatures that Numba does not cache on disk. Each measured build
still creates a new System and Network; system/network caches are not reused.
No runs execute concurrently.

The harness calls the ordinary `Command.run()` pipeline with:

```text
vorpy <absolute EDTA.pdb> -s nt pow -s mv 5 -s sr 0.2 -s bs 1.25 -e small
```

The CLI `small` export selects the **tiny** preset: system info, PDB, PyMOL atoms,
group info, shell surfaces and logs. Vertex search also writes its vertex file.
The input is loaded/classified through `System`, not the Workbench loader.
Output-directory metadata is redirected before input loading; the remaining
CLI behavior is preserved, including progress output, which is captured in a
per-worker log.

Scope is **32 selected atoms with 612 context balls**, producing 504 vertices,
896 edges, 425 surfaces and 32 complete cells. The remaining 580 context ball
records are marked incomplete. This is not the old vertex harness's
`--residues 0` scope, which selects all input atoms.

## Measurement semantics

- Every recorded call has its parent ID, wall/CPU start and end, inclusive
  duration, call count, immediate-child total, and unattributed remainder.
- For each clock, `parent = immediate children + remainder`. Category totals
  sum only these remainders, so nested calls are not double-counted. A leaf's
  remainder is the work measured by that boundary, not an additional mystery
  cost. Legacy nested timing dictionaries are retained separately and are not
  added into this hierarchy.
- Setup includes input/classification, settings and group/network creation.
  Spatial indexing is separated even if invoked within another timed stage.
- Network construction includes vertex search, topology, edge and surface
  geometry, plus the enclosing build remainder. Analysis and vertex output
  are subtracted into their own categories.
- Geometry preparation measures `Group.get_info()` and `Group.get_layers()`;
  nested layer calls are counted once. `file_export` contains exporter work
  after those measured preparation calls are removed, plus vertex-file output.
  It includes serialization and any preparation still inside exporter bodies;
  it is **not a pure disk-I/O measurement**.
- Pipeline total starts immediately before `Command()` and ends when
  `Command.run()` returns. It excludes imports, wrapper installation,
  correctness hashing and report generation. Those costs, and optional
  warm-up, are included in the separately retained GNU time worker totals.
- `/usr/bin/time -v` peak RSS covers the entire worker, including imports,
  snapshots and optional warm-up. It is not a per-stage memory allocation
  measurement. Linux process-group termination bounds the entire worker,
  including input, sorting, validation and warm-up.

## Results

All values below are seconds unless stated otherwise. Category values are
arithmetic means of three trials, so they reconcile to the mean pipeline total;
rounding may affect the last displayed digit.

| Exclusive category | Fresh wall | Warm wall | Warm CPU |
|---|---:|---:|---:|
| Input/setup | 0.011283 | 0.009018 | 0.009005 |
| Spatial indexing | 0.003362 | 0.003269 | 0.003264 |
| Network construction | 4.167910 | 2.893708 | 2.890523 |
| Network analysis | 0.064522 | 0.035764 | 0.035699 |
| Geometry preparation | 0.200604 | 0.202185 | 0.201884 |
| File export | 0.280959 | 0.342156 | 0.341918 |
| Pipeline orchestration remainder | 0.000061 | 0.000068 | 0.000082 |
| **Mean pipeline total** | **4.728701** | **3.486168** | **3.482375** |

Fresh pipeline range: **4.694–4.760 s**, median **4.732 s**.
Warm pipeline range: **3.457–3.525 s**, median **3.477 s**.
Maximum worker RSS: **287.3 MiB fresh**, **304.9 MiB warm** (warm includes two
builds and validation in one worker). Fresh full-worker elapsed time was
6.17–6.28 s; use `resource.txt` for full-process CPU and elapsed measurements.

Selected inclusive stage medians (these are nested, not another additive table):

| Stage | Fresh wall | Warm wall |
|---|---:|---:|
| Vertex search, including vertex-file output | 1.576 | 1.396 |
| Topology/connect | 0.061 | 0.012 |
| Edge geometry | 0.775 | 0.386 |
| Surface geometry | 1.700 | 1.117 |
| Entire Network.build, including analysis/output/sorting | 4.259 | 2.936 |

In warm runs the vertex stage visits 192,954 candidates and performs 3,968,382
verification-ball checks, with 722 site-container calls. Surrounding preparation
is approximately 0.249 s and candidate filtering approximately 0.180 s.
These are existing detail counters, not non-overlapping additions to the table.

## Important finding: first-call numerical compilation

`vorpy/src/calculations/calcs.py:calc_tetra_vol` uses `@jit(nopython=True)` without
`cache=True`. Surface construction calls it for the first ball and then the
second ball for each surface. The first volume bucket averaged **0.614 s fresh**
and **0.035 s warm**; the second stayed around **0.033 s**. This supports
first-call compilation as a major contributor to the apparent fresh-process
surface-volume hotspot. It does not establish that every cold/warm difference
comes from Numba; imports performed lazily, allocator state, garbage collection
and other initialization may also contribute.

Warm surface work still averages about 1.117 s, including about 0.824 s in
`build_surf` and 0.126 s in surface-area calculation. Vertex search is the largest
warm stage. The export category was higher in the warm series; this has not been
attributed. Do not assume all stages benefit from warm-up.

## Correctness and validation

Every measured run matched its control for:

- every column and row index of balls, vertices, edges and surfaces, preserving
  exact numerical values and sequence order;
- group membership, complete/incomplete counts, topology, surface triangles,
  geometry, contacts, curvature and energy fields present in those tables;
- all seven output files, with only trial-directory paths, known log timing
  fields/completion timestamp, and group-info timing lines normalized.

No numerical tolerances were loosened. Sets are canonicalized as unordered sets;
list ordering remains significant. This is a Power/EDTA preservation check, not
validation of AW throughput or nonzero curved-surface behavior.

Validation completed:

```text
Focused network + harness suite: 20 passed (includes ALA78 multi-edge regression)
compileall vorpy scripts: passed
Project-wide pytest: collection blocked by missing hypothesis and mesh imports
```

The harness tests cover nested wall/CPU reconciliation, repeated call counts,
exception cleanup, and preservation of numerical/output-order differences in
correctness checks. No packages were installed.

## Relationship to the previous audit

There is no EDTA entry in `network_build_baselines.csv`, so no historical
speedup ratio can be computed. Its Cambrin CLI row also has 17 fields for 16
headers; the harness reports that malformed row rather than interpreting it.
Existing historical CSVs are preserved. Several rows of the opportunities CSV
also need schema repair before machine consumption; they were not used here.

This environment is Python **3.12.3**, NumPy **2.0.2**, Numba **0.60.0**, SciPy
**1.18.1**, pandas **3.0.5**, and Shapely **2.0.7**, on an AMD Ryzen 9 3900X
(24 logical CPUs). The old audit lists Python 3.14.6, NumPy 2.5.2 and Numba 0.67.0.
That is another reason not to compare timings directly. The input hash, source
revision, harness hashes, actual settings and environment accompany the records.

## Recommended next task

Run a bounded EDTA Power `cProfile` comparison of fresh and same-process-warmed
pipeline calls, with particular attention to Numba compilation/cache loading in
surface volume and edge geometry. Use cumulative and self time to separate
startup from repeated geometry work; retain the exact existing correctness gate.
If compilation dominates the fresh CLI cost, test disk caching of eligible
existing Numba kernels as a small follow-up change before rewriting geometry.
For steady-state optimization, use the warmed vertex/surface profiles to choose
between candidate traversal and triangulation work. Do not select a rewrite or
parallel implementation from the fresh surface-volume bucket alone.

No AW, Cambrin, p53tet or NCP runs, package installation, multiprocessing of
network computation, or compiled-extension work were performed in this task.

## Retained evidence

- [Fresh-process summary](edta_power_pipeline/summary.json) and
  [CSV](edta_power_pipeline/summary.csv).
- [Warm-kernel summary](edta_power_pipeline/warm-kernels/summary.json) and
  [CSV](edta_power_pipeline/warm-kernels/summary.csv).
- Each dataset includes `control/` and `run-1/` through `run-3/` records with
  timing trees, legacy details, correctness hashes and GNU time reports.
  Warm records also retain their same-process uninstrumented warm-up result.
- Original generated files and console logs remain under
  `/tmp/vorpy-edta-power-audit-final` and `/tmp/vorpy-edta-power-warm-final`.
