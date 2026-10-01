import numpy as np
import pandas as pd
import pytest
from types import SimpleNamespace

from vorpy.src.boundary import (
    BoundaryConfig,
    BoundaryMode,
    has_explicit_solvent,
    resolve_boundary,
    virtual_shell,
    network_geometry,
)


def atoms(residues=("ALA", "GLY")):
    return pd.DataFrame({
        "loc": [np.array([0.0, 0.0, 0.0]), np.array([2.0, 1.0, -1.0])],
        "rad": [1.8, 1.6],
        "res_name": list(residues),
    })


def test_central_solvent_detection_and_default_mode():
    assert not has_explicit_solvent(atoms())
    assert has_explicit_solvent(atoms(("ALA", "hOh")))
    assert has_explicit_solvent(atoms(("NA", "ALA")))
    assert resolve_boundary(BoundaryConfig(), False).mode is BoundaryMode.NONE
    assert resolve_boundary(BoundaryConfig(), True).mode is BoundaryMode.NONE


def test_explicit_boundary_without_solvent_and_future_modes_fail_cleanly():
    with pytest.raises(ValueError, match="no explicit solvent was detected"):
        resolve_boundary(BoundaryConfig("explicit"), False)
    with pytest.raises(NotImplementedError, match="planned but not implemented"):
        resolve_boundary(BoundaryConfig("ses"), False)


def test_virtual_shell_is_deterministic_encloses_atomic_union_and_uses_water_oxygen_radius():
    structure = atoms()
    config = BoundaryConfig("shell", shell_spacing=2.5)
    first = virtual_shell(structure, config)
    second = virtual_shell(structure, config)
    assert first == second
    assert len(first) > 0
    assert all(site.is_boundary_generator and site.radius == 1.5 for site in first)
    coords = np.asarray([site.position for site in first])
    molecular_coords = np.asarray(list(structure["loc"]))
    molecular_radii = np.asarray(structure["rad"])
    lower = np.min(molecular_coords - molecular_radii[:, None], axis=0)
    upper = np.max(molecular_coords + molecular_radii[:, None], axis=0)
    assert np.all(coords.min(axis=0) < lower)
    assert np.all(coords.max(axis=0) > upper)


def test_network_geometry_keeps_virtual_sites_separate_from_physical_system():
    physical = atoms().assign(mass=[12.0, 14.0])
    system = SimpleNamespace(
        balls=physical,
        boundary_config=BoundaryConfig("shell"),
        has_explicit_solvent=False,
    )
    locs, radii, masses, boundary_indices, resolved = network_geometry(system)
    assert len(system.balls) == 2
    assert len(locs) == len(radii) == len(masses) == 2 + len(boundary_indices)
    assert boundary_indices == tuple(range(2, len(locs)))
    assert masses[:2] == [12.0, 14.0]
    assert masses[2:] == [0.0] * len(boundary_indices)
    assert resolved.mode is BoundaryMode.SHELL
    assert len(system.boundary_generators) == len(boundary_indices)


def test_explicit_solvent_keeps_network_geometry_unchanged():
    physical = atoms(("ALA", "HOH")).assign(mass=[12.0, 16.0])
    system = SimpleNamespace(
        balls=physical,
        boundary_config=BoundaryConfig(),
        has_explicit_solvent=True,
    )
    locs, radii, masses, boundary_indices, resolved = network_geometry(system)
    assert len(locs) == len(radii) == len(masses) == len(physical)
    assert boundary_indices == ()
    assert system.boundary_generators == ()
    assert resolved.mode is BoundaryMode.NONE


def test_cli_boundary_options_are_noninteractive_and_removed_from_legacy_args():
    from vorpy.src.command.vpy_cmnd import boundary_cli_options

    args, config, requested = boundary_cli_options([
        "molecule.pdb", "-b", "all", "--boundary", "shell",
        "--probe-radius", "1.6", "--shell-spacing", "2.0", "--shell-offset", "3.2",
    ])
    assert args == ["molecule.pdb", "-b", "all"]
    assert requested and config.mode is BoundaryMode.SHELL
    assert (config.probe_radius, config.shell_spacing, config.shell_offset) == (1.6, 2.0, 3.2)


def test_tiny_unsolvated_structure_solves_with_finite_physical_cells(tmp_path, monkeypatch):
    from vorpy.workbench.services.vorpy_backend import VorPyBackend, VorPySolveSettings

    path = tmp_path / "tetra.pdb"
    atoms_to_write = (
        ("C", (0., 0., 0.)), ("N", (1.4, 0.3, 0.2)), ("O", (-0.5, 1.5, 0.4)),
        ("S", (0.4, -1.2, 1.7)), ("C", (-1.5, -0.7, -0.3)), ("P", (0.8, 0.7, -1.6)),
    )
    lines = []
    for i, (element, (x, y, z)) in enumerate(atoms_to_write, start=1):
        lines.append(
            f"ATOM  {i:5d} {element:<4s} ALA A   1    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"  1.00  0.00          {element:>2s}\n"
        )
    path.write_text("".join(lines) + "END\n")
    backend = VorPyBackend(VorPySolveSettings(max_vertices=20, build_surfaces=False, boundary_mode="shell"))
    result = backend.solve(path, lambda *_: None, lambda: False)
    physical_cells = result.export_group.net.balls.iloc[:len(atoms_to_write)]
    assert physical_cells["sa"].gt(0).all()
    assert np.isfinite(physical_cells["vol"].to_numpy(dtype=float)).all()
    assert result.complete_cells == len(atoms_to_write)
    assert len(result.atoms) == len(atoms_to_write)
    assert sum(atom.mass for atom in result.atoms) > 0
    assert all(not atom.is_boundary_generator for atom in result.atoms)
    boundary_info = dict(result.info_sections["Boundary"])
    assert boundary_info["Mode"] == "Virtual solvent shell"
    assert int(boundary_info["Boundary generators"]) > 0

    # The network has virtual rows beyond the physical System table. Logs must
    # retain physical atom metadata without indexing those rows into System.
    from vorpy.src.output.logs import write_logs
    monkeypatch.chdir(tmp_path)
    write_logs(result.export_group)
    assert (tmp_path / f"{result.export_group.settings['net_type']}_logs.csv").exists()


@pytest.mark.parametrize("answer, expected", [("", "none"), ("0", "none"), ("1", "shell")])
def test_boundary_prompt_waits_for_terminal_choice(monkeypatch, answer, expected):
    from vorpy.src.command.vpy_cmnd import prompt_boundary
    monkeypatch.setattr("sys.stdin.isatty", lambda: True)
    calls = []
    monkeypatch.setattr("builtins.input", lambda prompt: calls.append(prompt) or answer)
    assert prompt_boundary(BoundaryConfig()).mode.value == expected
    assert len(calls) == 1


def test_boundary_batch_default_and_explicit_none(monkeypatch):
    from vorpy.src.command.vpy_cmnd import prompt_boundary, boundary_cli_options
    monkeypatch.setattr("sys.stdin.isatty", lambda: False)
    def unexpected_prompt(*args):
        pytest.fail("Batch runs must not prompt")
    monkeypatch.setattr("builtins.input", unexpected_prompt)
    assert resolve_boundary(prompt_boundary(BoundaryConfig()), False).mode is BoundaryMode.NONE
    _, config, requested = boundary_cli_options(["--boundary", "none"])
    assert requested and config.mode is BoundaryMode.NONE
    system = SimpleNamespace(balls=atoms().assign(mass=[12., 14.]),
                             boundary_config=config, has_explicit_solvent=False)
    locs, _, _, indices, _ = network_geometry(system)
    assert len(locs) == len(system.balls) and indices == ()
