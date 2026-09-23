import numpy as np
import importlib

from vorpy.src.network import slow

find_v0_module = importlib.import_module("vorpy.src.network.find_v0")


def _run_container(monkeypatch, limits, success_at=None, full=40.0):
    calls = []

    monkeypatch.setattr(slow, "box_search", lambda loc: 0)
    monkeypatch.setattr(slow, "get_balls", lambda cells, dist: [0, 1, 2, 3])

    def fake_find_site(**kwargs):
        calls.append(kwargs["max_vert"])
        if success_at is not None and kwargs["max_vert"] >= success_at:
            return ([{"balls": [0, 1, 2, 3], "loc": np.zeros(3), "rad": 1.0}], kwargs["invalid_ndxs"])
        return None, kwargs["invalid_ndxs"]

    monkeypatch.setattr(slow, "find_site", fake_find_site)
    metrics = {}
    result = slow.find_site_container_slow(
        [0, 1, 2], np.zeros((4, 3)), np.ones(4), [[] for _ in range(4)], [],
        limits, "aw", metrics=metrics, full_max_vert=full,
    )
    return result, calls, metrics


def test_aw_seed_retries_progressively_and_clamps(monkeypatch):
    result, calls, metrics = _run_container(monkeypatch, 4.0, success_at=30.0, full=30.0)

    assert result is not None
    assert metrics["seed_search_limits"] == [4.0, 8.0, 16.0, 30.0]
    assert metrics["seed_success_max_vert"] == 30.0
    assert calls[0] == 4.0
    assert 8.0 in calls and 16.0 in calls and 30.0 in calls


def test_aw_seed_stops_on_first_attempt(monkeypatch):
    result, calls, metrics = _run_container(monkeypatch, 4.0, success_at=4.0, full=40.0)

    assert result is not None
    assert metrics["seed_search_limits"] == [4.0]
    assert metrics["seed_success_max_vert"] == 4.0
    assert calls == [4.0]


def test_aw_seed_reaches_exact_final_limit(monkeypatch):
    result, calls, metrics = _run_container(monkeypatch, 4.0, success_at=None, full=40.0)

    assert result is None
    assert metrics["seed_search_limits"] == [4.0, 8.0, 16.0, 32.0, 40.0]
    assert max(calls) == 40.0


def test_find_v0_keeps_fast_initial_aw_limit(monkeypatch):
    observed = []
    candidate = {"balls": [0, 1, 2, 3], "loc": np.zeros(3), "rad": 1.0, "loc2": None}

    monkeypatch.setattr(find_v0_module, "get_balls", lambda cells, dist: [0, 1, 2, 3, 4])
    monkeypatch.setattr(find_v0_module, "box_search", lambda loc: 0)
    monkeypatch.setattr(find_v0_module, "calc_com", lambda locations: np.zeros(3))
    monkeypatch.setattr(find_v0_module, "calc_circ", lambda *args: (None, 1.0))

    def fake_container(*args, **kwargs):
        observed.append((kwargs["max_vert"], kwargs["full_max_vert"]))
        return [candidate, kwargs.get("metrics")]

    monkeypatch.setattr(find_v0_module, "find_site_container_slow", fake_container)
    result = find_v0_module.find_v0(
        locs=np.zeros((5, 3)), rads=np.ones(5), b_verts=[[] for _ in range(5)],
        max_vert=40.0, net_type="aw", b0=0, group_ndxs=list(range(5)),
    )

    assert result == candidate
    assert observed == [(4.0, 40.0)]
