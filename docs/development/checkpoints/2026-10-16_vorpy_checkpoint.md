# VorPy development checkpoint — 2026-10-16

## Checkpoint status

This document records the frozen shared-tree state for the targeted checkpoint
commit. The repository is on `wip/experimental-development` at
`361cc1dd5431f573d72eda02b4fb6ff39f785bb8`. No reset, clean, stash, rebase,
push, or molecular solve was performed.

PRISM's final report records the standalone progress fixture repair. The
targeted progress test now passes; production progress code was not changed.
The remaining dirty worktree is handled by an explicit staging allowlist below.
Excluded research data and generated artifacts remain untouched. The six-agent
freeze is in force for this checkpoint.

## Functionality completed in this cycle

- Full-source wet interface semantics and materialized source views are covered
  by bounded tests; selected interface exports are scoped to physical facets.
- Partial physical geometry is exported with typed, honest incomplete status;
  unsupported geometry remains fail-closed.
- `recover_network_source()` reloads and validates a checkpoint without calling
  `Network.build`.
- Benchmark preflight validates wet/dry hashes, group membership, parser plans,
  distinct destinations, resources, and competing VorPy processes.
- The preflight script now uses the analysis profile, timing report, `-i` plus
  `-e dir`, and explicit archive-save settings for command verification.
- PRISM's progress fixture now uses an explicitly interactive stream, preserving
  redraw behavior while keeping redirected output free of carriage returns.

## Inputs and baseline evidence

| Input | SHA-256 | Expected atoms |
|---|---|---:|
| `vorpy/data/1bna_DB1976.pdb` | `B24AAF71BBEC96A32A0E6BDBBE89175AD9C52A99459AB1A696B7B7DD9D66FF4D` | 10,494 |
| `vorpy/data/1bna_DB1976_no_sol.pdb` | `0243D7E3CC55E420C34CFFF31792A1398BEEC24336E48BB7954C1595AE8E5863` | 805 |
| `vorpy/data/1AAY.pdb` | `061D1861CFBE6D65DAAAB71B0B192CF8C81EC3B484A0BB6B1B5B295CBAE80EBD` | local input |
| `vorpy/data/2KAI.pdb` | `2E869BEE8330CFFD4DDFFE62F1E6063726FD1F8EBEE59F5834BF04411D275E6A` | local input |

The 1BNA wet/dry pair has 805 matching molecular coordinates (maximum delta
0.0 in the prepared preflight), DNA A+B has 758 atoms, and C/MOL has 47.
Wet-only solvent is 9,669 water atoms plus 20 sodium ions. Selenium inference
and radii are covered by chemistry tests.

The user's POW-dry retest solved 805 generators in 47.08 s and then failed
during Results/export after observing 191 physical surfaces, 512 edges, 322
vertices, and provisional area 179.944168 Å². Geometry was incomplete. The
preserved checkpoint is
`output/benchmarks/1bna_DB1976/retest_20261016_155326/pow_dry/1bna_DB1976_no_sol/1bna_DB1976_no_sol_full_system.geometry.vpy`
(15,194,871 bytes; SHA-256
`838A047F6CADD472AA821E98CA4E5A879D371C509FD39C21C1C7F5308E518F27`). An
earlier checkpoint remains at
`output/benchmarks/1bna_DB1976/pow_dry/1bna_DB1976_no_sol/`.

Historical wet timing observations were 635.4 s for POW and 1666.2 s for AW;
they are not a complete four-run phase-resolved baseline. No current AW/POW
wet/dry molecular run was launched for this checkpoint.

## Integration tests

Passed bounded tests:

- Focused archive/interface/recovery/export suite: **33 passed**.
- Current-tree CLI, EDTA, reuse, archive, Results, lineage, ordered-interface,
  output-preset, and source-scope matrix: **198 passed**.
- Watchdog process/failure handling: **10 passed**.
- Broader integration matrix before the final progress fixture repair: **276
  passed, 1 failed, 1 warning**.
- Exact PRISM standalone progress test after the repair: **1 passed**.
- Benchmark preflight against a fresh destination: **passed**; wet 10,494,
  dry 805, matching molecular coordinates, parser-valid plans, no competing
  VorPy process.
- Checkpoint recovery smoke test: **passed**; 805 generators recovered,
  525 complete cells, 191 surfaces, 512 edges, 322 vertices, no rebuild.

FINAL TARGETED CHECKPOINT RUN:

- The repaired standalone progress test passed: **1 passed**.
- The bounded staged-source interface/CLI/archive/output matrix completed with
  exit code 0; it exercised chemistry inference, 1BNA command parsing,
  materialized AW/POW views, partial exports, Results/lineage, archive
  recovery, output reservation, timing, watchdogs, and output routing.
- A deliberately broader sweep was not used as the commit gate: it reached the
  repository timeout in the pre-existing 1AAY molecular-group case, and a GUI
  scientific-session import reached the optional PyVista/VTK stack. Those
  cases were not rerun and are not represented as passing checkpoint evidence.
- Staged Python syntax validation passed for 89 files, and
  `git diff --cached --check` passed after removing three new blank-at-EOF
  errors. The staged exclusion scan found no research input, output, archive,
  report, virtual-environment, or generated-mesh path.

The previous failure in
`vorpy/tests/network/test_views.py::test_standalone_progress_pads_replaced_line`
is resolved in the test fixture by using an interactive stream. `git diff
--check` on the full worktree still reports trailing whitespace in modified
PDB inputs (expected fixed-width PDB padding); those inputs are excluded from
the checkpoint commit and are not normalized.

Skipped/not attempted: all molecular solves, all 1BNA full solves, trajectory
frames, 1AAY/2KAI solves, large archive reloads, full solvent visualization,
performance optimization, and origin push. These are deliberate checkpoint
restrictions, not passing results.

## Scientific limitations and performance

Partial exports report observed supported geometry only; missing geometry,
global topology, total interface area, and unavailable curvature remain
unresolved rather than zero. A source-derived view is not claimed equivalent
to a direct solve under a different generator universe. The failed wet AW run
reached 75,887 vertices before termination and did not complete network
construction; no mathematically equivalent localized algorithm has been
demonstrated, so no optimization claim is made.

## Archive and excluded local artifacts

The source checkpoint above is recoverable read-only through
`recover_network_source()` and was not rewritten. Generated `output/`, `.vpy`
archives, logs, OFF/CSV/PDB exports, `.agent-reports/`, local benchmark
manifests/transcripts, virtual environments, scratch data, and untracked
molecular inputs remain local and are intentionally excluded from this
checkpoint commit. Root `results.txt` is also intentionally excluded.

Important local paths include `output/benchmarks/`, `output/1bna_DB1976_*`,
`output/1AAY_*`, `output/2KAI_*`, `vorpy/data/1BNA_DB1976_frames/`,
`vorpy/data/1AAY_frames/`, `vorpy/data/2KAI_frames/`, and the preserved
checkpoint paths above. They must be copied separately and hash-verified.

## Reproduction commands

Use PowerShell from the repository root with `.venv\\Scripts\\python.exe`.
The read-only preflight is:

```powershell
& .\\.venv\\Scripts\\python.exe scripts\\preflight_1bna_db1976_benchmark.py --manifest docs\\benchmarks\\1bna_db1976\\benchmark_manifest.json --output-root output\\benchmarks\\1bna_DB1976\\retest_20261016_export_repair_v1 --repeat-id retest_20261016_export_repair_v1
```

The bounded focused matrix is:

```powershell
& .\\.venv\\Scripts\\python.exe -m pytest -q vorpy\\tests\\test_cli_edta_exports.py vorpy\\tests\\test_cli_1bna_build_mode.py vorpy\\tests\\test_full_system_reuse.py vorpy\\tests\\test_nexus_end_to_end_checkpoint.py vorpy\\tests\\test_output_presets.py vorpy\\tests\\output\\test_source_view_export_scope.py vorpy\\tests\\io\\test_interface_result_archive.py vorpy\\tests\\io\\test_network_recovery.py vorpy\\tests\\network\\test_pow_source_view_geometry.py vorpy\\tests\\test_nexus_aw_interface_integration.py vorpy\\tests\\test_nexus_interface_output_contract.py vorpy\\tests\\test_nexus_pow_watchdog.py
```

The progress gate is:

```powershell
& .\\.venv\\Scripts\\python.exe -m pytest -q vorpy\\tests\\network\\test_views.py -k standalone_progress
```

The four benchmark commands remain manual-only and must use the versioned
manifest/output roots; they were not executed here.

## Environment and recovery

The project requires Python `>=3.9`; the prepared environment is the repository
`.venv`. Dependencies are declared in `pyproject.toml`: NumPy, SciPy,
Matplotlib, pandas, Shapely, Pillow, Numba, with pytest/pytest-cov and related
development tools. On another machine, install the project and development
extras in a fresh environment, then copy the excluded input/checkpoint files
from the listed local paths and verify the hashes in this document before any
solve. Do not assume Git contains those files.

## Next steps and ownership

1. NEXUSA: stage and review the explicit source/test/documentation allowlist,
   rerun its focused matrix, and run `git diff --check` on the staged set.
2. ATLAS: revalidate checkpoint/archive round trips from the committed source.
3. FORGE/GAUSS/VECTOR: retain source-view, completeness, and physical export
   contracts; no solver optimization during baseline retest.
4. After approval: manually run one isolated benchmark at a time, record
   per-phase timings and memory, then compare against the historical wet
   timings.

The intended commit message is
`checkpoint: preserve VorPy interface and benchmark development`.
