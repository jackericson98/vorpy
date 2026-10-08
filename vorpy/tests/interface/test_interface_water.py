"""Regression tests for interface buried-water group classification."""

from types import SimpleNamespace

from vorpy.src.interface import water as water_module
from vorpy.src.interface.export import export_buried_water_groups


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


def test_empty_water_analysis_and_export_do_not_create_placeholder_directory(tmp_path):
    iface = SimpleNamespace(
        dir=str(tmp_path / "interface"),
        name="A_B_interface",
        water_topology={"buried_cycle_sets": []},
        buried_water_groups=[],
        _buried_water_exports_complete=False,
        sys=None,
    )

    assert water_module.build_buried_water_groups(iface) == []
    export_buried_water_groups(iface)

    assert not (tmp_path / "interface" / "waters").exists()


def test_buried_water_analysis_defers_directory_creation_until_export(tmp_path, monkeypatch):
    residue = SimpleNamespace(name="HOH", seq=1, atoms=[0])

    class FakeGroup:
        def __init__(self, **kwargs):
            self.dir = kwargs["output_directory"]
            self.settings = {}
            self.net = None

        def build(self):
            return None

        def get_info(self):
            return None

    monkeypatch.setattr(water_module, "Group", FakeGroup)
    monkeypatch.setattr(water_module, "_cycle_set_residues", lambda cycle, topology: [residue])
    monkeypatch.setattr(water_module, "_buried_group_name", lambda index, residues: "buried_1")
    monkeypatch.setattr(water_module, "_residue_label", lambda value: "HOH:1")
    monkeypatch.setattr(water_module, "_buried_group_color", lambda index: "viridis")
    monkeypatch.setattr(water_module, "summarize_water_group", lambda group: {})

    iface = SimpleNamespace(
        dir=str(tmp_path / "interface"),
        name="A_B_interface",
        interface_id="A_B",
        sys=SimpleNamespace(update_progress=None),
        group1=SimpleNamespace(settings={}),
        water_topology={
            "buried_cycle_sets": [{
                "water_keys": {"water-1"},
                "edge_indices": [1],
                "cycle_set_id": 1,
            }],
            "waters": {},
        },
    )

    groups = water_module.build_buried_water_groups(iface)

    assert len(groups) == 1
    assert not (tmp_path / "interface" / "waters").exists()
