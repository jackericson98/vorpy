# EDTA edge curvature optimization

Measured 2026-09-12 against the fused AW kernel at commit
17580535cb222c37f2b0fc62ff0fff51840a8a55.

The ordinary `.venv/bin/python vorpy edta` path builds 32 complete cells,
500 vertices, 892 edges and 424 surfaces. This is AW; the older Power timing
report measures a different calculation.

## Change

Batch the directed face cross products at each quadrature point and clip the
three mean-curvature cosines together. Preserve individual dot products,
quadrature order, accumulation order, orientation logic and validation checks.
No new dependencies or numerical tolerance changes.

Also fix verbose kernel overhead accounting: subtract the timing increments
for the current edge, rather than cumulative measurements from previous edges.

## Evidence

Captured all 892 calls to `aw_edge_curvature_measures` during the normal EDTA
CLI pipeline. Replayed the same geometry and quadrature inputs through the
original and revised functions, with timing instrumentation disabled. Each
version received one warm-up and three timed passes, run sequentially.

| Kernel | Three passes (s) | Median (s) |
| --- | --- | --- |
| Original | 4.7972, 4.8301, 4.8368 | 4.8301 |
| Batched | 3.3759, 3.4277, 3.3778 | 3.3778 |

The kernel took **30.1% less time (1.43x speedup)**. All returned mean and
Gaussian dictionaries matched exactly for all 892 edges. This measures kernel
replay, excluding network topology lookup, shared geometry resolution and
exports; it is not an end-to-end speedup claim. The original cProfile replay
recorded 88,824 `np.cross` calls with 3.467 s cumulative time out of 7.791 s
profiled execution. Profiling times are separate from the table above.

Both the original requested CLI command and the revised full CLI run completed.
The latter used `-e dir /tmp/vorpy-edta-curvature-after` to isolate exports.

Validation: 71 focused tests passed:

```bash
.venv/bin/python -m pytest -q vorpy/tests/calculations/test_edge_geometry.py vorpy/tests/network/test_build_net.py
```

Added comparisons against independent mean and oriented Gaussian integrals for
straight/curved edges in both directions, plus deterministic repeated-call timer
accounting coverage. The existing ALA78 topology regression is included.

Session artifacts (temporary): `/tmp/profile_edta_edges.py`,
`/tmp/compare_edta_edges.py`, `/tmp/vorpy-edge-geometry-before.py`,
`/tmp/edta-edge-inputs.pkl`, `/tmp/edta-edges.prof`,
`/tmp/edta-edge-comparison.log`, and `/tmp/vorpy-edta-curvature-after.log`.
The capture script wraps the network module's kernel entry point and calls
`Command().run()` with `edta`; the comparison script reuses the current geometry
classes for both kernel implementations. These temporary files are local
benchmark evidence, not portable fixtures.
