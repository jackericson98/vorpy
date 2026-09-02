"""Regression tests for interface buried-water group classification."""

from vorpy.src.interface import water as water_module


def test_closed_cycle_set_is_buried(monkeypatch):
    graph = {
        "edge_indices": [1, 2, 3, 4],
        "open_vertices": [],
    }
    monkeypatch.setattr(
        water_module,
        "_analyze_component_graph",
        lambda iface, keys, waters: graph,
    )

    cycle_set = {
        "water_keys": {("water", 1), ("water", 2)},
        "edge_indices": [1, 2, 3, 4],
    }

    result = water_module._classify_cycle_set_group(
        iface=None,
        cycle_set=cycle_set,
        waters={},
    )

    assert result["burial_class"] == "buried"
    assert result["extra_noncycle_edge_indices"] == []
    assert result["open_vertex_indices"] == []


def test_cycle_with_open_branch_is_semi_buried(monkeypatch):
    graph = {
        "edge_indices": [1, 2, 3, 4, 5],
        "open_vertices": [99],
    }
    monkeypatch.setattr(
        water_module,
        "_analyze_component_graph",
        lambda iface, keys, waters: graph,
    )

    cycle_set = {
        "water_keys": {("water", 1), ("water", 2)},
        "edge_indices": [1, 2, 3, 4],
    }

    result = water_module._classify_cycle_set_group(
        iface=None,
        cycle_set=cycle_set,
        waters={},
    )

    assert result["burial_class"] == "semi_buried"
    assert result["extra_noncycle_edge_indices"] == [5]
    assert result["open_vertex_indices"] == [99]
