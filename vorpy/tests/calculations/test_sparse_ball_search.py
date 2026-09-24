import itertools

import pytest

from vorpy.src.calculations.sorting import get_balls


@pytest.mark.parametrize("cells, dist, extra", [
    ([[0, 0, 0]], 0, 0),
    ([[7, 8, 9]], 2.5, 0),
    ([[2, 4, 6], [12, 10, 9]], 0.5, 2),
    ([[19, 19, 19]], 100, 0),
    ([[40, 40, 40]], 0, 0),
])
def test_sparse_search_matches_dense_iteration_order(cells, dist, extra):
    n = 20
    # Deliberately reverse insertion order and preserve order inside each cell.
    keys = list(itertools.product(range(0, n, 3), repeat=3))[::-1]
    matrix = {key: [2 * i + 1, 2 * i] for i, key in enumerate(keys)}
    matrix[-1, -1, -1] = [n]
    size = [0.5, 1.0, 2.0]
    reach = int(dist / min(size)) + 3
    bounds = [range(max(0, min(c[a] for c in cells) - reach - extra),
                    min(n, max(c[a] for c in cells) + reach + extra)) for a in range(3)]
    expected = [ball for key in itertools.product(*bounds) for ball in matrix.get(key, [])]
    assert get_balls(cells, dist, extra, matrix, size, 1) == expected


def test_large_sparse_search_never_probes_empty_cells():
    class OccupiedOnly(dict):
        def __getitem__(self, key):
            assert key in self, "Large search probed an empty grid cell"
            return super().__getitem__(key)

    matrix = OccupiedOnly({(-1, -1, -1): [1000], (999, 999, 999): [9], (0, 0, 0): [0]})
    assert get_balls([500, 500, 500], 1000, 0, matrix, [1, 1, 1], 1) == [0, 9]
