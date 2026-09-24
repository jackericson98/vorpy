import numpy as np
import pytest

from vorpy.src.network import fast


@pytest.mark.parametrize("count", [32, 2500, 4096, 4097])
@pytest.mark.parametrize("radius", [-1.5, 0.0, 2.0])
def test_aw_verification_matches_full_clearance(monkeypatch, count, radius):
    rng = np.random.default_rng(16)
    locs = rng.uniform(-20, 20, (count, 3))
    rads = rng.uniform(0.2, 2, count)
    defining = [0, 2, 5, 8]
    monkeypatch.setattr(fast, "box_search", lambda point: [0, 0, 0])

    def neighborhood(*args, **kwargs):
        assert count > 4096, "Small systems must not rebuild neighborhoods"
        return list(range(count))

    monkeypatch.setattr(fast, "get_balls", neighborhood)
    cache = {}
    for point in (np.zeros(3), locs[7], np.full(3, 30.0)):
        other = np.ones(count, dtype=bool)
        other[defining] = False
        expected = np.all(np.linalg.norm(locs[other] - point, axis=1) - rads[other] >= radius)
        for search_cache in (None, cache):
            assert fast.verify_aw_local(point, radius, defining, locs, rads, max(rads), search_cache) == expected
    if count <= 4096:
        assert "aw_verification" not in cache


def test_aw_verification_skips_definers_and_checks_box(monkeypatch):
    locs = np.zeros((4, 3))
    rads = np.ones(4)
    monkeypatch.setattr(fast, "box_search", lambda point: [0, 0, 0])
    assert fast.verify_aw_local(np.zeros(3), 1, [0, 1, 2, 3], locs, rads, 1)
    monkeypatch.setattr(fast, "box_search", lambda point: None)
    assert not fast.verify_aw_local(np.zeros(3), 1, [0, 1, 2, 3], locs, rads, 1)


def test_choose_vertex_keeps_valid_secondary_when_primary_is_blocked(monkeypatch):
    monkeypatch.setattr(fast, "box_search", lambda point: [0, 0, 0])
    locs = np.zeros((5, 3))
    rads = np.ones(5)
    candidate = {"balls": [0, 1, 2, 3], "loc": np.zeros(3), "rad": 1,
                 "loc2": np.array([10., 0., 0.]), "rad2": 1}
    result, rejected = fast.choose_vert(candidate, [0, 1, 2], None, locs, rads, None, max_ball_rad=1)
    assert rejected is None
    assert np.array_equal(result[0]["loc"], [10, 0, 0])
    assert result[0]["loc2"] is None
