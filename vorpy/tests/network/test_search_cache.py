import numpy as np

from vorpy.src.network import fast
from vorpy.src.network.search_cache import BoundedCache, retained_size, cache_bucket


def test_byte_budget_evicts_oldest_and_hits_refresh_recency():
    value = np.zeros((100, 3))
    entry_size = retained_size(('a', value)) + 128
    cache = BoundedCache(max_bytes=entry_size * 2, max_entries=10)
    cache['a'] = value
    cache['b'] = value.copy()
    assert cache['a'] is value
    cache['c'] = value.copy()
    assert 'a' in cache and 'c' in cache and 'b' not in cache
    assert cache.evictions == 1
    assert cache.bytes <= cache.max_bytes


def test_oversized_geometry_is_returned_without_retaining_it():
    cache = {'geometry': BoundedCache(max_bytes=1)}
    calls = []
    def calculate():
        calls.append(1)
        return np.ones((100, 3))
    for _ in range(2):
        np.testing.assert_array_equal(fast._cached_geometry(cache, 'key', calculate), 1)
    assert len(calls) == 2
    assert len(cache['geometry']) == 0
    assert cache['geometry'].skipped == 2


def test_replacements_and_entry_limit_keep_accounting_bounded():
    cache = BoundedCache(max_bytes=10000, max_entries=2)
    cache['a'] = np.zeros(10)
    cache['a'] = np.zeros(20)
    assert cache.bytes == retained_size(('a', cache['a'])) + 128
    for i in range(20):
        cache[i] = i
        assert len(cache) <= 2
        assert cache.bytes <= cache.max_bytes
    assert cache.evictions > 0


def test_numpy_views_count_retained_base_allocation():
    base = np.zeros((1000, 3))
    view = base[:1]
    assert retained_size(view) >= base.nbytes


def test_aw_results_match_uncached_after_neighborhood_eviction(monkeypatch):
    rng = np.random.default_rng(18)
    locs = rng.uniform(-20, 20, (4097, 3))
    rads = np.ones(4097)
    monkeypatch.setattr(fast, 'ball_search_reach', lambda distance: distance)
    monkeypatch.setattr(fast, 'box_search', lambda point: [0, 0, 0])
    monkeypatch.setattr(fast, 'get_balls', lambda *args, **kwargs: list(range(4097)))
    cache = {'aw_verification': BoundedCache(max_bytes=1024 * 1024, max_entries=1)}
    for radius in (0., 1., 2., 0.):
        for point in (np.zeros(3), np.full(3, 30.)):
            expected = fast.verify_aw_local(point, radius, [0, 1, 2, 3], locs, rads, 1, None)
            actual = fast.verify_aw_local(point, radius, [0, 1, 2, 3], locs, rads, 1, cache)
            assert actual == expected
            assert cache['aw_verification'].bytes <= cache['aw_verification'].max_bytes
    assert cache['aw_verification'].evictions >= 3


def test_grid_equivalent_radii_reuse_arrays_but_keep_exact_clearance(monkeypatch):
    from vorpy.src.calculations import sorting
    monkeypatch.setattr(sorting, 'sub_box_size', [10., 10., 10.], raising=False)
    monkeypatch.setattr(fast, 'box_search', lambda point: [0, 0, 0])
    calls = []
    def get_balls(*args, **kwargs):
        calls.append(kwargs['dist'])
        return list(range(4097))
    monkeypatch.setattr(fast, 'get_balls', get_balls)
    locs = np.full((4097, 3), 100.)
    locs[4] = [1.5, 0., 0.]
    rads = np.ones(4097)
    cache = {}
    assert fast.verify_aw_local(np.zeros(3), .4, [0, 1, 2, 3], locs, rads, 1., cache)
    assert not fast.verify_aw_local(np.zeros(3), .6, [0, 1, 2, 3], locs, rads, 1., cache)
    assert len(calls) == 1
    assert not fast.verify_aw_local(np.zeros(3), 9.1, [0, 1, 2, 3], locs, rads, 1., cache)
    assert len(calls) == 2


def test_fast_sizes_conservatively_cover_search_payloads():
    balls = list(range(1000))
    arrays = (np.zeros((1000, 3)), np.ones(1000), {ball: i for i, ball in enumerate(balls)})
    for name, value in [('candidates', balls), ('aw_verification', arrays),
                        ('surrounding', (balls,) + arrays)]:
        bucket = cache_bucket({}, name)
        key = ((0, 1, 2), 3.)
        bucket[key] = value
        assert bucket.bytes >= retained_size((key, value)) + 128
        assert bucket.bytes <= bucket.max_bytes
