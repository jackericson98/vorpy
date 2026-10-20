import pandas as pd
import pytest

from vorpy.src.boundary import (
    BoundaryGenerator,
    BoundaryConfig,
    SolventMode,
    prepare_network_geometry,
)
from vorpy.src.group import Group
from vorpy.src.interface import Interface
from vorpy.src.network import Network
from vorpy.src.network.views import NetworkView
from vorpy.src.group.shared_view import analyze_group_view
from vorpy.src.system.system import System
from vorpy.src.command.vpy_cmnd import solvent_cli_options
from vorpy.workbench.services.vorpy_backend import VorPySolveSettings


def _balls():
    return pd.DataFrame({
        "loc": [(0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (0.0, 1.5, 0.0)],
        "rad": [1.7, 1.5, 1.5],
        "mass": [12.0, 16.0, 18.0],
        "num": [100, 101, 102],
        "res_name": ["ALA", "HOH", "GLY"],
        "name": ["CA", "O", "CA"],
    })


def test_molecule_only_preserves_system_rows_but_excludes_solvent_generators():
    system = System(spheres=_balls(), name="solvent-mode", make_dir=False)

    geometry = prepare_network_geometry(system)

    assert system.solvent_mode is SolventMode.MOLECULE_ONLY
    assert len(system.balls) == 3
    assert geometry.system_indices == (0, 2)
    assert geometry.solvent_local_indices == ()
    assert geometry.solvent_generator_ids == ()
    assert geometry.solvent_cells_complete is False


def test_explicit_modes_record_membership_and_completeness():
    system = System(spheres=_balls(), name="solvent-mode", make_dir=False)

    context = prepare_network_geometry(system, SolventMode.SOLVENT_CONTEXT)
    assert context.system_indices == (0, 1, 2)
    assert context.solvent_local_indices == (1,)
    assert context.participating_solvent_generator_digest
    assert context.solvent_cells_complete is False

    full = prepare_network_geometry(system, SolventMode.FULL_SOLVENT)
    network = Network(
        locs=full.locations,
        rads=full.radii,
        masses=full.masses,
        group=[0, 1, 2],
        settings={"net_type": "aw", "build_type": "all"},
        generator_system_indices=full.system_indices,
        generator_stable_ids=full.stable_ids,
        solvent_local_indices=full.solvent_local_indices,
        solvent_generator_ids=full.solvent_generator_ids,
        solvent_mode=full.solvent_mode,
    )
    assert network.solvent_provenance["solvent_mode"] == "FULL_SOLVENT"
    assert network.solvent_provenance["solvent_cells_complete"] is True
    assert network.settings["solvent_generator_digest"]
    assert network.generator_membership_evidence[1][1] == 1
    assert network.generator_membership_evidence[1][3] is True


def test_group_default_network_maps_local_generators_to_system_ids():
    system = System(spheres=_balls(), name="solvent-mode", make_dir=False)
    group = Group(
        system,
        name="solute",
        atoms=[0, 2],
        make_net=True,
        build_net=False,
    )

    assert len(group.net.balls) == 2
    assert group.net.group == [0, 1]
    assert group.net.generator_system_indices == (0, 2)
    assert group.net.solvent_provenance["solvent_mode"] == "MOLECULE_ONLY"
    assert group.net.solvent_provenance["solvent_generator_digest"] is None

    with pytest.raises(ValueError, match="excluded by solvent_mode"):
        Group(
            system,
            name="water",
            atoms=[1],
            make_net=True,
            build_net=False,
        )


def test_boundary_identity_is_local_and_deterministic():
    system = System(
        spheres=_balls(),
        name="solvent-mode",
        make_dir=False,
        boundary_config=BoundaryConfig("shell", shell_spacing=4.0),
    )
    first = prepare_network_geometry(system)
    second = prepare_network_geometry(system)

    assert first.system_indices[:2] == (0, 2)
    assert first.system_indices[2:] == (None,) * len(first.boundary_indices)
    assert first.stable_ids[2:] == tuple(
        f"boundary:{index}" for index in range(len(first.boundary_indices))
    )
    assert first.boundary_indices == second.boundary_indices
    assert first.stable_ids == second.stable_ids


def test_cli_requires_explicit_opt_in_for_non_default_solvent_modes():
    remaining, mode = solvent_cli_options(["input.pdb", "--solvent-context"])
    assert remaining == ["input.pdb"]
    assert mode is SolventMode.SOLVENT_CONTEXT

    remaining, mode = solvent_cli_options(
        ["input.pdb", "--solvent-mode", "FULL_SOLVENT"]
    )
    assert remaining == ["input.pdb"]
    assert mode is SolventMode.FULL_SOLVENT


def test_workbench_settings_share_cli_solvent_and_retention_contract():
    defaults = VorPySolveSettings()
    assert defaults.solvent_mode == SolventMode.MOLECULE_ONLY.value
    assert defaults.retained_system_indices == ()

    settings = VorPySolveSettings(
        solvent_mode="solvent-context",
        retained_system_indices=(2, 0, 2),
    )
    assert settings.solvent_mode == SolventMode.SOLVENT_CONTEXT.value
    assert settings.retained_system_indices == (0, 2)


def test_retained_solvent_generator_is_included_in_group_scope():
    system = System(
        spheres=_balls(),
        name="retained-solvent",
        make_dir=False,
        retained_system_indices=(1,),
    )
    group = Group(
        system,
        name="solute-with-retained-water",
        atoms=[0, 2],
        make_net=True,
        build_net=False,
    )

    assert group.net.generator_system_indices == (0, 1, 2)
    assert group.net.retained_system_indices == (1,)
    assert group.net.retained_generator_ids == ("system:1",)
    assert group.net.get_solve_context_key().digest


def test_interface_scope_rejects_overlapping_group_membership():
    with pytest.raises(ValueError, match="overlap"):
        Network(
            locs=[(0.0, 0.0, 0.0)] * 3,
            rads=[1.0] * 3,
            group=[0, 1, 2],
            iface_grps=({0, 1}, {1, 2}),
            settings={"net_type": "aw"},
        )


def test_dedicated_interface_network_preserves_distinct_group_sides(tmp_path):
    balls = _balls().copy()
    balls["res_name"] = ["GLY", "ALA", "SER"]
    system = System(
        spheres=balls,
        name="synthetic-interface",
        output_directory=str(tmp_path),
        make_dir=False,
    )
    group_a = Group(
        system,
        name="chains-A-B",
        atoms=[0, 1],
        make_net=False,
        build_net=False,
    )
    group_b = Group(
        system,
        name="chain-I",
        atoms=[2],
        make_net=False,
        build_net=False,
    )
    system._compact_interface_workflow = True
    interface = Interface(system, group_a, group_b)
    network = interface.make_net()

    assert network.iface_grps == ({0, 1}, {2})
    assert set(network.group) == {0, 1, 2}
    assert network.network_mode == "interface"


def test_full_solvent_interface_keeps_water_as_global_competitor(tmp_path):
    system = System(
        spheres=_balls(),
        name="full-solvent-interface",
        output_directory=str(tmp_path),
        make_dir=False,
        solvent_mode=SolventMode.FULL_SOLVENT,
    )
    group_a = Group(
        system, name="solute-A", atoms=[0], make_net=False, build_net=False,
    )
    group_b = Group(
        system, name="solute-B", atoms=[2], make_net=False, build_net=False,
    )
    system._compact_interface_workflow = True
    interface = Interface(system, group_a, group_b)
    network = interface.make_net()

    assert network.generator_system_indices == (0, 1, 2)
    assert network.solvent_local_indices == (1,)
    assert network.group == [0, 2]
    assert network.iface_grps == ({0}, {2})
    assert network.solvent_mode == SolventMode.FULL_SOLVENT.value


def test_full_solvent_interface_keeps_water_and_ion_competitors_outside_owners(tmp_path):
    balls = _balls().copy()
    balls.loc[2, "res_name"] = "NA"
    system = System(
        spheres=balls,
        name="full-solvent-water-ion-interface",
        output_directory=str(tmp_path),
        make_dir=False,
        solvent_mode=SolventMode.FULL_SOLVENT,
    )
    group_a = Group(system, name="solute-A", atoms=[0], make_net=False, build_net=False)
    group_b = Group(system, name="solute-B", atoms=[1], make_net=False, build_net=False)
    system._compact_interface_workflow = True
    interface = Interface(system, group_a, group_b)
    network = interface.make_net()

    assert network.generator_system_indices == (0, 1, 2)
    assert network.solvent_local_indices == (1, 2)
    assert network.solvent_generator_ids == ("system:1", "system:2")
    assert network.group == [0, 1]
    assert network.iface_grps == ({0}, {1})
    assert set(interface.group1_indices) | set(interface.group2_indices) == {0, 1}
    assert set(network.generator_system_indices) - {0, 1} == {2}


def test_molecule_only_view_accepts_incomplete_numerical_boundary_support():
    locations = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (1.0, 1.0, 2.0)]
    balls = pd.DataFrame({
        "loc": locations,
        "rad": [1.5, 1.5, 1.5],
        "num": [0, 1, 2],
        "system_num": [10, 11, -1],
        "stable_id": ["atom-10", "atom-11", "boundary:2"],
        "is_boundary_generator": [False, False, True],
        "complete": [True, True, False],
        "vol": [3.0, 4.0, 0.0],
    })
    vertices = pd.DataFrame({
        "loc": [(0.0, 0.0, 0.5)] * 3,
        "rad": [0.2] * 3,
        "balls": [[0, 1, 2]] * 3,
        "complete": [True] * 3,
    })
    edges = pd.DataFrame({
        "balls": [[0, 1, 2]] * 3,
        "verts": [[0, 1], [1, 2], [2, 0]],
        "complete": [True] * 3,
    })
    surfaces = pd.DataFrame({
        "balls": [[0, 1, 2]],
        "verts": [[0, 1, 2]],
        "edges": [[0, 1, 2]],
        "points": [[(0.0, 0.0, 0.0)] * 3],
        "tris": [[(0, 1, 2)]],
        "sa": [2.0],
        "complete": [True],
    })
    network = Network(
        locs=locations,
        rads=[1.5] * 3,
        group=[0, 1],
        settings={"net_type": "prm", "surf_res": 0.2, "max_vert": 5.0},
        balls=balls,
        verts=vertices,
        edges=edges,
        surfs=surfaces,
        system=type("StaticSystem", (), {"frame_index": 0})(),
        solvent_mode="MOLECULE_ONLY",
    )
    network.boundary_indices = (2,)
    network.boundary_generators = (
        BoundaryGenerator(position=locations[2]),
    )

    view = NetworkView.for_group(network, [10])
    result = analyze_group_view(view)

    assert result.complete is True
    assert result.system_ids == (10,)
    assert result.area == pytest.approx(2.0)
    assert result.volume == pytest.approx(3.0)
    assert view.boundary_generator_ids == ("boundary:2",)


def test_solved_molecule_only_network_supports_two_views_without_rebuild(monkeypatch):
    rows = [
        ((0.0, 0.0, 0.0), 1.7, 12.0, "C"),
        ((2.0, 0.0, 0.0), 1.5, 14.0, "N"),
        ((0.0, 2.0, 0.0), 1.5, 16.0, "O"),
        ((0.0, 0.0, 2.0), 1.5, 32.0, "S"),
        ((-1.8, -1.0, -1.0), 1.5, 12.0, "C"),
        ((2.0, 2.0, 2.0), 1.5, 14.0, "N"),
        ((1.0, 1.0, 1.0), 1.4, 18.0, "HOH"),
    ]
    balls = pd.DataFrame({
        "loc": [row[0] for row in rows],
        "rad": [row[1] for row in rows],
        "mass": [row[2] for row in rows],
        "res_name": [row[3] for row in rows],
        "num": list(range(len(rows))),
        "name": ["A"] * len(rows),
    })
    system = System(
        spheres=balls,
        name="f10-solved",
        make_dir=False,
        boundary_config=BoundaryConfig("shell", shell_spacing=4.0),
    )
    group = Group(
        system,
        name="all-solute",
        atoms=list(range(6)),
        make_net=True,
        build_net=False,
    )
    group.build(calculate_curvature=False)
    network = group.net
    balls_before = network.balls.copy(deep=True)

    monkeypatch.setattr(
        Network,
        "build",
        lambda *args, **kwargs: pytest.fail("view construction rebuilt Network"),
    )
    first = NetworkView.for_group(network, [0, 1, 2])
    second = NetworkView.for_group(network, [3, 4, 5])

    assert network.solvent_mode == SolventMode.MOLECULE_ONLY.value
    assert network.solvent_local_indices == ()
    assert first.completeness.proven and second.completeness.proven
    assert first.system_ids == (0, 1, 2)
    assert second.system_ids == (3, 4, 5)
    assert first.selection_id != second.selection_id
    assert first.solve_universe_id == second.solve_universe_id
    assert first.local_to_system[:6] == tuple((index, index) for index in range(6))
    pd.testing.assert_frame_equal(network.balls, balls_before)


def test_f10_metadata_survives_valid_molecular_archive_round_trip(tmp_path):
    from vorpy.src.io.network_archive import load_network, save_network

    atoms = (
        ("ATOM", "C", "ALA", "A", 1, (0.0, 0.0, 0.0)),
        ("ATOM", "N", "ALA", "A", 1, (2.0, 0.0, 0.0)),
        ("ATOM", "O", "ALA", "A", 1, (0.0, 2.0, 0.0)),
        ("ATOM", "S", "ALA", "A", 1, (0.0, 0.0, 2.0)),
        ("ATOM", "C", "ALA", "A", 1, (-1.8, -1.0, -1.0)),
        ("ATOM", "N", "ALA", "A", 1, (2.0, 2.0, 2.0)),
        ("HETATM", "O", "HOH", "W", 2, (1.0, 1.0, 1.0)),
    )
    path = tmp_path / "f10.pdb"
    lines = [
        f"{record:<6s}{index:5d} {name:<4s} {residue:>3s} {chain}   {sequence:1d}    "
        f"{loc[0]:8.3f}{loc[1]:8.3f}{loc[2]:8.3f}  1.00  0.00          {name:>2s}\n"
        for index, (record, name, residue, chain, sequence, loc) in enumerate(atoms, start=1)
    ]
    path.write_text("".join(lines) + "END\n")

    system = System(
        file=str(path),
        output_directory=str(tmp_path),
        make_dir=False,
        boundary_config=BoundaryConfig("shell", shell_spacing=4.0),
    )
    group = Group(
        system,
        name="all-solute",
        atoms=list(range(6)),
        make_net=True,
        build_net=False,
    )
    group.build(calculate_curvature=False)
    network = group.net
    context_before_save = network.get_solve_context_key(frame_id="static")
    context_digest = context_before_save.digest
    archive = tmp_path / "f10.vpy"
    save_network(network, archive)
    loaded = load_network(archive)

    assert loaded.settings["solvent_mode"] == "MOLECULE_ONLY"
    assert loaded.settings["solvent_generator_digest"] is None
    assert loaded.settings["solvent_cells_complete"] is False
    loaded_context = loaded.get_solve_context_key(frame_id="static")
    assert loaded_context.digest == context_digest
    assert tuple(loaded.balls["system_num"].iloc[:6]) == tuple(range(6))
    assert len(loaded.balls) == len(network.balls)
