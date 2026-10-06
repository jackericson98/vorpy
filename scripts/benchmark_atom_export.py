"""Small export-preparation benchmark; does not solve a molecular network.

Run from the repository root: python scripts/benchmark_atom_export.py
Reference mode restores repeated row lookups and unused OFF face fields.
"""
import json
import sys
from pathlib import Path
from time import perf_counter
from types import SimpleNamespace

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from vorpy.src.output import curvature_colors, edges, verts


def main():
    locations = np.random.default_rng(7).normal(size=(256, 3))
    net = SimpleNamespace(
        settings={'net_type': 'aw', 'surf_scheme': 'int_mean_curv',
                  'surf_col': 'coolwarm', 'scheme_factor': 'log'},
        surfs=pd.DataFrame(),
        edges=pd.DataFrame([
            {'balls': [0, 1, 2], 'points': [locations[i % 256], locations[(i + 1) % 256]],
             'int_mean_curv_by_ball': {0: 1., 1: -1., 2: 2.}}
            for i in range(512)]),
        verts=pd.DataFrame([
            {'loc': location, 'edges': [(i + j) % 512 for j in range(4)]}
            for i, location in enumerate(locations)]),
    )
    # Warm compiled drawing code before measuring either mode.
    edges.prepare_edges(net, list(range(512)), color_limit=3.)
    results = {}
    reference_samples = {}
    for reference in (True, False):
        label = 'reference' if reference else 'optimized'
        results[label] = {}
        with curvature_colors.export_color_cache(net):
            if reference:
                del net._export_component_rows
            for name, prepare, size in [('edges', edges.prepare_edges, 512),
                                        ('vertices', verts.prepare_verts, 256)]:
                start = perf_counter()
                for cell in range(120):
                    indices = [(cell + offset) % size for offset in range(32)]
                    mesh = prepare(net, indices, color_limit=3., target_cells=[cell % 2],
                                   color_mode='cell', include_face_data=reference)
                    if cell == 0:
                        if reference:
                            reference_samples[name] = mesh
                        else:
                            expected = reference_samples[name]
                            for field in ('points', 'triangles', 'face_colors'):
                                np.testing.assert_array_equal(getattr(mesh, field), getattr(expected, field))
                results[label][name] = perf_counter() - start
            if reference:
                net._export_component_rows = {}
    results['speedup'] = {name: results['reference'][name] / results['optimized'][name]
                          for name in ('edges', 'vertices')}
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
