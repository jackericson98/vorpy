"""Bounded benchmark for the pre-fusion and fused surface geometry paths."""

from __future__ import annotations

import subprocess
import sys
import time
from statistics import median

import numpy as np
from numba import jit

from vorpy.src.calculations import calc_surf_sa, calc_tetra_vol
from vorpy.src.network.build_surfs import _surface_tetra_volumes


WARM_REPETITIONS = 7
STRICT_TOLERANCE = 1.0e-14


@jit(nopython=True)
def _legacy_volume_kernel(references, points, triangles):
    """Volume-only kernel used before area was fused into the traversal."""
    totals = np.zeros(len(references), dtype=np.float64)
    for triangle in triangles:
        p1 = points[triangle[0]]
        p2 = points[triangle[1]]
        p3 = points[triangle[2]]
        for reference_index in range(len(references)):
            p0 = references[reference_index]
            r01x, r01y, r01z = p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]
            r02x, r02y, r02z = p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]
            r03x, r03y, r03z = p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]
            cross_x = r01y * r02z - r01z * r02y
            cross_y = r01z * r02x - r01x * r02z
            cross_z = r01x * r02y - r01y * r02x
            triple_product = r03x * cross_x + r03y * cross_y + r03z * cross_z
            totals[reference_index] += abs(triple_product) / 6.0
    return totals


def _legacy_separate(references, points, triangles):
    references = np.asarray(references, dtype=float)
    points = np.asarray(points, dtype=float)
    triangles = np.asarray(triangles, dtype=np.int64)
    volumes = _legacy_volume_kernel(references, points, triangles)
    area = calc_surf_sa(tris=triangles, points=points)
    return volumes, area


def _flat_surface(nx=96, ny=96):
    x, y = np.meshgrid(
        np.linspace(-2.0, 2.0, nx + 1),
        np.linspace(-2.0, 2.0, ny + 1),
        indexing="ij",
    )
    points = np.column_stack((x.ravel(), y.ravel(), np.zeros(x.size)))
    triangles = []
    stride = ny + 1
    for i in range(nx):
        for j in range(ny):
            a = i * stride + j
            b, c, d = a + 1, a + stride, a + stride + 1
            triangles.extend(((a, b, c), (b, d, c)))
    return points, np.asarray(triangles, dtype=np.int64)


def _curved_surface(n_lat=72, n_lon=144, radius=2.0):
    points = [np.array([0.0, 0.0, radius])]
    for i in range(1, n_lat):
        theta = np.pi * i / n_lat
        for j in range(n_lon):
            phi = 2.0 * np.pi * j / n_lon
            points.append(
                [
                    radius * np.sin(theta) * np.cos(phi),
                    radius * np.sin(theta) * np.sin(phi),
                    radius * np.cos(theta),
                ]
            )
    south = len(points)
    points.append(np.array([0.0, 0.0, -radius]))

    triangles = []
    for j in range(n_lon):
        triangles.append((0, 1 + j, 1 + (j + 1) % n_lon))
    for i in range(n_lat - 2):
        upper = 1 + i * n_lon
        lower = upper + n_lon
        for j in range(n_lon):
            a, b = upper + j, upper + (j + 1) % n_lon
            c, d = lower + j, lower + (j + 1) % n_lon
            triangles.extend(((a, c, b), (b, c, d)))
    last = 1 + (n_lat - 2) * n_lon
    for j in range(n_lon):
        triangles.append((last + j, south, last + (j + 1) % n_lon))
    return np.asarray(points, dtype=float), np.asarray(triangles, dtype=np.int64)


def _expected_volumes(references, points, triangles):
    return np.asarray(
        [
            sum(
                calc_tetra_vol(
                    reference,
                    points[triangle[0]],
                    points[triangle[1]],
                    points[triangle[2]],
                )
                for triangle in triangles
            )
            for reference in references
        ]
    )


def _timed(call):
    start = time.perf_counter()
    result = call()
    return time.perf_counter() - start, result


def _benchmark(name, points, triangles):
    references = np.asarray(
        [[0.25, -0.5, 1.75], [2.0, 1.5, -0.75]], dtype=float
    )
    legacy_cold, legacy_result = _timed(
        lambda: _legacy_separate(references, points, triangles)
    )
    fused_cold, fused_result = _timed(
        lambda: _surface_tetra_volumes(references, points, triangles)
    )

    expected_area = calc_surf_sa(tris=triangles, points=points)
    expected_volumes = _expected_volumes(references, points, triangles)
    np.testing.assert_allclose(
        fused_result[1], expected_area, rtol=STRICT_TOLERANCE, atol=STRICT_TOLERANCE
    )
    np.testing.assert_allclose(
        fused_result[0], expected_volumes,
        rtol=STRICT_TOLERANCE,
        atol=STRICT_TOLERANCE,
    )
    np.testing.assert_allclose(
        fused_result[1], legacy_result[1],
        rtol=STRICT_TOLERANCE,
        atol=STRICT_TOLERANCE,
    )

    def warm_median(call):
        call()
        samples = [_timed(call)[0] for _ in range(WARM_REPETITIONS)]
        return median(samples)

    legacy_warm = warm_median(lambda: _legacy_separate(references, points, triangles))
    fused_warm = warm_median(
        lambda: _surface_tetra_volumes(references, points, triangles)
    )
    improvement = 100.0 * (legacy_warm - fused_warm) / legacy_warm
    print(
        f"{name}: points={len(points):,}, triangles={len(triangles):,}; "
        f"cold legacy={legacy_cold:.6f}s, fused={fused_cold:.6f}s; "
        f"warm median legacy={legacy_warm:.6f}s, fused={fused_warm:.6f}s; "
        f"improvement={improvement:+.1f}%"
    )
    print(
        f"  max area diff={abs(fused_result[1] - expected_area):.3e}; "
        f"max volume diff={np.max(abs(fused_result[0] - expected_volumes)):.3e}"
    )


def _run_case(name):
    factory = _flat_surface if name == "flat" else _curved_surface
    _benchmark(name, *factory())


if __name__ == "__main__":
    if len(sys.argv) == 2:
        _run_case(sys.argv[1])
    else:
        for case in ("flat", "curved"):
            subprocess.run([sys.executable, __file__, case], check=True)
