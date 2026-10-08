from types import SimpleNamespace

import pytest


def _fake_system(tmp_path):
    calls = {"system": [], "group": [], "interface": [], "canonical": [], "archives": []}
    group_net = SimpleNamespace(settings={"net_type": "pow"})
    interface_net = SimpleNamespace(settings={"net_type": "pow"})
    group = SimpleNamespace(
        name="group_1",
        net=group_net,
        settings={"net_type": "pow"},
        dir=None,
        exports=lambda **kwargs: calls["group"].append(kwargs),
    )
    interface = SimpleNamespace(
        name="interface_group_1_group_2",
        net=interface_net,
        dir=str(tmp_path / "interface_group_1_group_2"),
        export=lambda **kwargs: calls["interface"].append(kwargs),
    )
    system = SimpleNamespace(
        name="fixture",
        groups=[group],
        ifaces=[interface],
        files={"dir": str(tmp_path)},
        exports=lambda **kwargs: calls["system"].append(kwargs),
        update_progress=lambda **kwargs: None,
    )
    return system, calls


@pytest.mark.parametrize("preset", ["micro", "tiny"])
def test_micro_and_tiny_are_metadata_only_for_all_network_scopes(tmp_path, monkeypatch, preset):
    output = __import__("vorpy.src.output.output", fromlist=["output"])
    io = __import__("vorpy.src.io", fromlist=["io"])
    system, calls = _fake_system(tmp_path)
    monkeypatch.setattr(io, "save_network", lambda network, path: calls["archives"].append(path))
    monkeypatch.setattr(output, "_move_vert_file", lambda *args: None)
    monkeypatch.setattr(output, "_export_canonical_interface_log", lambda **kwargs: None)

    output.export_preset(system, preset)

    assert calls["archives"] == [str(tmp_path / "fixture.vpy")]
    assert calls["system"] == [{"info": True}]
    assert calls["group"] == [{"info": True}]
    assert calls["interface"] == [{"info": True, "buried_water": False}]


def test_system_only_preset_still_saves_archive(tmp_path, monkeypatch):
    output = __import__("vorpy.src.output.output", fromlist=["output"])
    io = __import__("vorpy.src.io", fromlist=["io"])
    system = SimpleNamespace(
        name="system_only",
        groups=[],
        ifaces=[],
        files={"dir": str(tmp_path)},
        exports=lambda **kwargs: None,
        update_progress=lambda **kwargs: None,
    )
    archives = []
    monkeypatch.setattr(io, "save_network", lambda network, path: archives.append(path))

    output.export_preset(system, "micro")

    assert archives == [str(tmp_path / "system_only.vpy")]


def test_small_exports_full_geometry_logs_and_no_atom_visualization(tmp_path, monkeypatch):
    output = __import__("vorpy.src.output.output", fromlist=["output"])
    io = __import__("vorpy.src.io", fromlist=["io"])
    system, calls = _fake_system(tmp_path)
    monkeypatch.setattr(io, "save_network", lambda network, path: calls["archives"].append(path))
    monkeypatch.setattr(output, "_move_vert_file", lambda *args: None)
    monkeypatch.setattr(output, "_export_canonical_interface_log",
                        lambda **kwargs: calls["canonical"].append(kwargs))

    output.export_preset(system, "small")

    assert calls["archives"] == [str(tmp_path / "fixture.vpy")]
    assert calls["system"] == [{"info": True}]
    assert {"info", "surfs", "logs"} == {
        next(iter(kwargs)) for kwargs in calls["group"]
    }
    assert {"info", "surfs"} == {
        next(iter(kwargs)) for kwargs in calls["interface"]
    }
    interface_info = next(kwargs for kwargs in calls["interface"] if "info" in kwargs)
    assert interface_info["buried_water"] is False
    assert len(calls["canonical"]) == 1
    assert not any("atoms" in kwargs or "set_atoms" in kwargs or "pdb" in kwargs
                   for kwargs in calls["group"] + calls["interface"])


def test_richer_presets_keep_physical_surface_and_atom_script_pairing():
    output = __import__("vorpy.src.output.output", fromlist=["output"])

    for preset in ("medium", "large", "all"):
        assert any(kwargs.get("pdb") for _, kwargs in output.SYSTEM_PRESETS[preset])
        assert any(kwargs.get("set_atoms") for _, kwargs in output.SYSTEM_PRESETS[preset])
        assert any(kwargs.get("surfs") for _, kwargs in output.GROUP_PRESETS[preset])
        assert any(kwargs.get("surfs") for _, kwargs in output.INTERFACE_PRESETS[preset])


def test_rich_dual_derivatives_drop_preset_pymol_launchers(tmp_path):
    import json

    output = __import__("vorpy.src.output.output", fromlist=["output"])
    dual = tmp_path / "dual"

    def export(**kwargs):
        dual.mkdir()
        (dual / "apollonius_interface.pml").write_text("run app\n", encoding="utf-8")
        (dual / "apollonius_interface.py").write_text("from pymol import cmd\n", encoding="utf-8")
        (dual / "interface_dual_summary.json").write_text(
            json.dumps({"pymol_script": str(dual / "apollonius_interface.pml")}),
            encoding="utf-8",
        )

    iface = SimpleNamespace(dir=str(tmp_path), export=export)
    output._export_dual_geometry_without_launcher(iface)

    assert not list(dual.glob("*.pml"))
    assert not list(dual.glob("*.py"))
    assert json.loads((dual / "interface_dual_summary.json").read_text())["pymol_script"] is None
