from pathlib import Path
import csv
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


def test_interface_default_preset_uses_first_class_layout_without_archive(tmp_path):
    from vorpy.src.output.output import export_interface_default
    from vorpy.tests.test_nexus_interface_output_contract import _contract_case

    _, interface = _contract_case(tmp_path)
    interface.sys.balls["rad"] = 1.7
    interface.sys.exports = lambda **kwargs: None
    interface.sys.update_progress = lambda **kwargs: None

    export_interface_default(interface.sys)

    directory = tmp_path / "interface"
    assert {
        "interface_surfs.off", "interface_edges.off", "interface_verts.off",
        "interface_info.txt", "group_a_atoms.pdb", "group_b_atoms.pdb",
        "surr_atoms.pdb", "surr_residues.pdb",
    } <= {path.name for path in directory.iterdir()}
    assert not (directory / "visualization_manifest.json").exists()
    assert not list(tmp_path.glob("*.geometry.vpy"))
    assert not (tmp_path / "synthetic-interface.vpy").exists()


def test_interface_default_reuses_cache_and_emits_classic_logs(tmp_path, monkeypatch):
    from vorpy.src.interface import export as interface_export
    from vorpy.src.network import Network
    from vorpy.src.output.export_scope import ExportScope
    from vorpy.src.output.output import export_interface_default
    from vorpy.tests.test_nexus_interface_output_contract import _contract_case

    optimized_root = tmp_path / "optimized"
    optimized_root.mkdir()
    _, interface = _contract_case(optimized_root)
    interface.sys.exports = lambda **kwargs: None
    interface.sys.update_progress = lambda **kwargs: None
    calls = {"scope": 0, "mesh": 0}
    original_scope = interface_export._interface_physical_scope
    original_mesh = ExportScope.mesh_network

    def counted_scope(*args, **kwargs):
        calls["scope"] += 1
        return original_scope(*args, **kwargs)

    def counted_mesh(scope, *args, **kwargs):
        calls["mesh"] += 1
        return original_mesh(scope, *args, **kwargs)

    monkeypatch.setattr(interface_export, "_interface_physical_scope", counted_scope)
    monkeypatch.setattr(ExportScope, "mesh_network", counted_mesh)
    monkeypatch.setattr(Network, "build", lambda *args, **kwargs: pytest.fail(
        "export rebuilt an existing network"
    ))

    export_interface_default(interface.sys)
    optimized = Path(interface.dir)
    assert calls == {"scope": 1, "mesh": 1}
    assert (optimized / "logs.csv").is_file()
    classic = optimized / "aw_logs.csv"
    assert classic.is_file()
    rows = list(csv.reader(classic.open(newline="", encoding="utf-8")))
    atoms_header = rows[rows.index(["Atoms"]) + 1]
    assert atoms_header[:5] == ["Index", "Name", "Chain", "Residue", "Residue Sequence"]
    assert "Complete Cell?" in atoms_header
    assert "Maximum Mean Curvature" in atoms_header

    reference_root = tmp_path / "reference"
    reference_root.mkdir()
    _, reference = _contract_case(reference_root)
    from vorpy.src.interface.export import interface_exports
    interface_exports(
        reference,
        atoms=True,
        surfs=True,
        edges=True,
        verts=True,
        logs=True,
        info=True,
    )
    reference_dir = Path(reference.dir)
    equivalent = (
        "interface_surfs.off", "interface_edges.off", "interface_verts.off",
        "interface_info.txt", "network_info.txt",
        "interface_surfaces.csv", "interface_edges.csv", "interface_vertices.csv",
        "group_a_atoms.pdb", "group_b_atoms.pdb", "surr_atoms.pdb", "surr_residues.pdb",
    )
    for name in equivalent:
        if name.endswith(".txt"):
            actual = (optimized / name).read_text(encoding="utf-8")
            expected = (reference_dir / name).read_text(encoding="utf-8")
            actual = actual.replace(str(optimized_root), "<ROOT>")
            expected = expected.replace(str(reference_root), "<ROOT>")
            assert actual == expected
        else:
            assert (optimized / name).read_bytes() == (reference_dir / name).read_bytes()


def test_interface_cli_default_and_explicit_exports_stay_separate(monkeypatch):
    from vorpy.src.command.command_export import argv_export
    from vorpy.src.output import output

    default_calls = []
    monkeypatch.setattr(
        output, "export_interface_default",
        lambda system: default_calls.append(system),
    )
    default_system = SimpleNamespace(_compact_interface_workflow=True)
    argv_export(default_system, [])
    assert default_calls == [default_system]

    explicit_calls = []
    interface = SimpleNamespace(
        name="A_B_interface", net=SimpleNamespace(), dir=None,
        export=lambda **kwargs: explicit_calls.append(kwargs),
    )
    explicit_system = SimpleNamespace(
        groups=[], ifaces=[interface], update_progress=lambda **kwargs: None,
    )
    output.other_exports(explicit_system, "surfaces")
    assert explicit_calls == [{"surfs": True}]


def test_explicit_interface_preset_parses_with_directory_and_exports_files(tmp_path):
    from vorpy.src.command.vpy_cmnd import Command
    from vorpy.tests.test_nexus_interface_output_contract import _contract_case

    _, interface = _contract_case(tmp_path)
    system = interface.sys
    system.balls["rad"] = 1.7
    system.exports = lambda **kwargs: None
    system.update_progress = lambda **kwargs: None
    system.verbose = False
    system.file_type = "off"

    destination = tmp_path / "output" / "benchmarks" / "1bna_DB1976" / "pow_dry"
    command = Command(sys=system)
    command.base_file = str(tmp_path / "synthetic_interface.pdb")
    command._arguments = [
        command.base_file, "-i", "-e", "interface", "and", "dir", str(destination)
    ]
    command.parse_commands()
    command._prepare_interface_output_directories()
    command.run_exports()

    assert command.exports == [["interface"]]
    assert Path(system.files["dir"]) == destination / system.name
    assert Path(interface.dir) == destination / system.name / "dna_ligand_interface"
    assert {
        "interface_surfs.off", "interface_edges.off", "interface_verts.off",
        "interface_info.txt", "group_a_atoms.pdb", "group_b_atoms.pdb",
        "surr_atoms.pdb", "surr_residues.pdb",
    } <= {path.name for path in Path(interface.dir).iterdir()}
    assert not (destination / "visualization_manifest.json").exists()
    assert not list(destination.glob("*.geometry.vpy"))


def test_interface_directory_and_internal_checkpoint_are_segregated(tmp_path, monkeypatch):
    from vorpy.src.command.vpy_cmnd import Command

    interface = SimpleNamespace(name="DNA_c_2_interface")
    system = SimpleNamespace(
        name="1bna_DB1976", files={"dir": str(tmp_path)}, ifaces=[interface],
    )
    command = Command(sys=system)
    command._prepare_interface_output_directories()
    assert Path(interface.dir) == tmp_path / "DNA_c_2_interface"
    refusing = Command(sys=system)
    refusing.fail_if_exists = True
    with pytest.raises(SystemExit, match="existing interface output directory"):
        refusing._prepare_interface_output_directories()

    captured = {}
    monkeypatch.setattr(
        "vorpy.src.io.save_network",
        lambda network, path, **kwargs: captured.setdefault("path", Path(path)) or path,
    )
    command.interface_mode = True
    command.execution_profile = None
    network = SimpleNamespace(group_name="1bna_DB1976_full_system", archive_provenance={})
    checkpoint = command._persist_geometry_checkpoint(network, scope="full_system_source")

    assert checkpoint == captured["path"]
    assert checkpoint.parent == tmp_path / ".vorpy" / "checkpoints"
    assert not (tmp_path / "1bna_DB1976_full_system.geometry.vpy").exists()
