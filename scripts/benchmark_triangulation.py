"""Compare small surface meshes against the committed performance baseline."""
import importlib
import json
from pathlib import Path
import subprocess
import sys
from time import perf_counter
from types import ModuleType

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main():
    baseline_code = subprocess.check_output([
        'git', '-c', f'safe.directory={ROOT.as_posix()}', 'show',
        'efcb7d4f:vorpy/src/network/triangulate.py'], cwd=ROOT, text=True)
    baseline = ModuleType('triangulation_baseline')
    exec(compile(baseline_code, '<committed triangulation baseline>', 'exec'), baseline.__dict__)
    optimized = importlib.import_module('vorpy.src.network.triangulate')
    surfaces = []
    for index in range(600):
        angles = np.arange(32 + index % 24) * (2 * np.pi / (32 + index % 24))
        radius = .8 + .1 * (index % 10) + .12 * np.cos(3 * angles)
        surfaces.append([np.asarray(point) for point in np.column_stack(
            (radius * np.cos(angles), radius * .7 * np.sin(angles)))])
    results, reference = {}, []
    for label, module in [('baseline', baseline), ('optimized', optimized)]:
        timing = {}
        start = perf_counter()
        for index, perimeter in enumerate(surfaces):
            points, triangles = module.triangulate_2D_Surface(perimeter, .2, timing=timing)
            if label == 'baseline':
                reference.append((points, triangles))
            else:
                # Verify outside the measured interval below.
                expected_points, expected_triangles = reference[index]
                verify_start = perf_counter()
                np.testing.assert_array_equal(points, expected_points)
                np.testing.assert_array_equal(triangles, expected_triangles)
                start += perf_counter() - verify_start
        results[label] = {'seconds': perf_counter() - start,
                          'grid': timing['tri_spiderweb'], 'classification': timing['tri_sort'],
                          'filtering': timing['tri_point_filter'], 'qhull': timing['tri_delaunay']}
    results['speedup'] = results['baseline']['seconds'] / results['optimized']['seconds']
    results['identical_meshes'] = len(surfaces)
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
