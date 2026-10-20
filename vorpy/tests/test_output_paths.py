from pathlib import Path
from types import SimpleNamespace

import pytest

from vorpy.src.command.command_export import argv_export
from vorpy.src.command.vpy_cmnd import Command
from vorpy.src.output.output import _set_group_directory
from vorpy.src.output.set_system_directory import default_output_root, set_sys_dir


def _system(name="fixture", root=None):
    return SimpleNamespace(
        name=name,
        files={"dir": None, "vpy_dir": None if root is None else str(root)},
    )


def test_default_output_is_project_root_output_from_another_cwd(tmp_path, monkeypatch):
    project = tmp_path / "project"
    other_cwd = tmp_path / "other-cwd"
    project.mkdir()
    other_cwd.mkdir()
    system = _system(root=project)

    monkeypatch.chdir(other_cwd)
    set_sys_dir(system)

    destination = project / "output" / "fixture"
    assert Path(system.files["dir"]) == destination
    assert destination.is_dir()
    assert not (other_cwd / "fixture").exists()
    assert not (other_cwd / "output").exists()


def test_missing_vpy_dir_falls_back_to_repository_output_root():
    root = default_output_root(_system())

    assert root == Path(__file__).resolve().parents[2] / "output"


@pytest.mark.parametrize("destination_kind", ["relative", "absolute"])
def test_explicit_cli_directory_is_created_at_requested_location(
    tmp_path, monkeypatch, destination_kind
):
    working = tmp_path / "working"
    working.mkdir()
    monkeypatch.chdir(working)
    requested = (
        Path("relative-results")
        if destination_kind == "relative"
        else tmp_path / "absolute-results"
    )
    system = _system()
    command = Command(sys=system)
    command._arguments = [
        "fixture.pdb", "-e", "dir", str(requested),
    ]
    command.parse_commands()

    expected_root = (
        working / requested if destination_kind == "relative" else requested
    )
    expected = expected_root / "fixture"
    assert Path(system.files["dir"]) == expected
    assert expected.is_dir()


def test_interface_and_group_exports_stay_beneath_resolved_root(tmp_path):
    root = tmp_path / "output" / "fixture"
    root.mkdir(parents=True)
    interface = SimpleNamespace(name="A_B_interface", dir=None)
    system = SimpleNamespace(files={"dir": str(root)}, ifaces=[interface])
    command = Command(sys=system)
    command._prepare_interface_output_directories()

    group = SimpleNamespace(name="A", ball_ndxs=[0], dir=None)
    system.balls = [0, 1]
    system.groups = [group]
    _set_group_directory(system, group)
    Path(interface.dir, "interface_surfs.off").write_text("OFF\n", encoding="ascii")
    Path(group.dir, "group_atoms.pdb").write_text("ATOM\n", encoding="ascii")

    assert Path(interface.dir) == root / "A_B_interface"
    assert Path(group.dir) == root / "A"
    assert (Path(interface.dir) / "interface_surfs.off").is_file()
    assert (Path(group.dir) / "group_atoms.pdb").is_file()
    assert not (tmp_path / "A_B_interface").exists()
    assert not (tmp_path / "A").exists()


def test_interface_collision_and_failed_export_do_not_escape_root(tmp_path):
    root = tmp_path / "output" / "fixture"
    existing = root / "A_B_interface"
    existing.mkdir(parents=True)
    interface = SimpleNamespace(name="A_B_interface", dir=None)
    system = SimpleNamespace(files={"dir": str(root)}, ifaces=[interface])
    command = Command(sys=system)
    command._prepare_interface_output_directories()

    assert Path(interface.dir) == root / "A_B_interface_2"

    def failed_export():
        Path(interface.dir, "interface_info.txt").write_text(
            "partial\n", encoding="ascii"
        )
        raise RuntimeError("synthetic export failure")

    with pytest.raises(RuntimeError, match="synthetic export failure"):
        failed_export()

    assert (Path(interface.dir) / "interface_info.txt").is_file()
    assert not (tmp_path / "interface_info.txt").exists()


def test_implicit_interface_export_does_not_use_visualization_or_cwd(monkeypatch):
    from vorpy.src.output import output

    calls = []
    monkeypatch.setattr(
        output, "export_interface_default", lambda system: calls.append(system)
    )
    system = SimpleNamespace(_compact_interface_workflow=True)

    argv_export(system, [])

    assert calls == [system]
