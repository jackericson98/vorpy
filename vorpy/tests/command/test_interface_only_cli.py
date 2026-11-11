import pandas as pd
import pytest
from pathlib import Path

from vorpy.src.command.vpy_cmnd import Command
from vorpy.src.output.layout import canonical_interfaces_directory, slug
from vorpy.src.system.system import System


@pytest.mark.parametrize("scheme", ["aw", "pow", "prm"])
def test_interface_only_cli_builds_only_direct_interface_network(
    tmp_path, capsys, scheme,
):
    system = System(
        spheres=pd.DataFrame({
            "loc": [
                (0.0, 0.0, 0.0), (2.0, 0.0, 0.0),
                (0.0, 2.0, 0.0), (0.0, 0.0, 2.0),
                (0.8, 0.8, 0.8),
            ],
            "rad": [0.75, 0.5, 0.75, 0.75, 0.35],
            "mass": [1.0] * 5,
            "num": list(range(5)),
            "res_name": ["ALA", "ALA", "ALA", "ALA", "HOH"],
            "name": ["CA"] * 5,
        }),
        name="synthetic-interface-cli",
        output_directory=str(tmp_path),
        make_dir=False,
        solvent_mode="MOLECULE_ONLY",
    )
    command = Command(sys=system)
    command.fail_if_exists = True
    command.base_file = str(tmp_path / "synthetic.pdb")
    command._arguments = [
        command.base_file,
        "-g", "a", "0",
        "-g", "a", "4",
        "-i",
        "-s", "nt", scheme, "and", "mv", "40",
        "--profile", "geometry",
        "--boundary", "none",
    ]
    command.parse_commands()
    command.run_exports = lambda: None
    command._save_archive = lambda: None

    command._process_system()

    assert command._full_system_group is None
    assert len(system.groups) == 2
    assert all(group.net is None for group in system.groups)
    assert len(system.ifaces) == 1

    network = system.ifaces[0].net
    assert network.network_mode == "interface"
    assert network.generator_system_indices == (0, 1, 2, 3, 4)
    assert network.solvent_local_indices == (4,)
    assert network.iface_grps == ({0}, {4})
    assert set(network.group) == {0, 4}
    assert Path(system.ifaces[0].dir).is_dir()

    output = capsys.readouterr().out
    assert f"Solving {scheme.upper()} interface network..." in output
    assert f"Solving {scheme.upper()} network..." not in output


def test_interface_only_cli_rejects_preexisting_destination_before_constructor(
    tmp_path,
):
    system = System(
        spheres=pd.DataFrame({
            "loc": [
                (0.0, 0.0, 0.0), (2.0, 0.0, 0.0),
                (0.0, 2.0, 0.0), (0.0, 0.0, 2.0),
            ],
            "rad": [0.75, 0.5, 0.75, 0.75],
            "mass": [1.0] * 4,
            "num": list(range(4)),
            "res_name": ["ALA"] * 4,
            "name": ["CA"] * 4,
        }),
        name="synthetic-interface-collision",
        output_directory=str(tmp_path),
        make_dir=False,
        solvent_mode="MOLECULE_ONLY",
    )
    command = Command(sys=system)
    command.fail_if_exists = True
    command.base_file = str(tmp_path / "synthetic.pdb")
    command._arguments = [
        command.base_file,
        "-g", "a", "0",
        "-g", "a", "1",
        "-i",
        "-s", "nt", "pow", "and", "mv", "40",
        "--profile", "geometry",
        "--boundary", "none",
    ]
    command.parse_commands()
    command.apply_settings()
    command.create_groups()
    group_a, group_b = command._planned_interface_pairs()[0]
    interface_id = "__".join(sorted((group_a.group_id, group_b.group_id)))
    destination = canonical_interfaces_directory(system) / slug(interface_id)
    destination.mkdir(parents=True)
    command.create_groups = lambda: None
    command.run_exports = lambda: None
    command._save_archive = lambda: None

    with pytest.raises(SystemExit, match="existing interface output directory"):
        command._process_system()
    assert system.ifaces == []


def test_surrounding_interface_reservation_matches_materialized_interface_id(
    tmp_path,
):
    system = System(
        spheres=pd.DataFrame({
            "loc": [
                (0.0, 0.0, 0.0), (2.0, 0.0, 0.0),
                (0.0, 2.0, 0.0), (0.0, 0.0, 2.0),
                (0.8, 0.8, 0.8),
            ],
            "rad": [0.75, 0.5, 0.75, 0.75, 0.35],
            "mass": [1.0] * 5,
            "num": list(range(5)),
            "res_name": ["ALA", "ALA", "ALA", "ALA", "HOH"],
            "name": ["CA"] * 5,
        }),
        name="synthetic-surrounding-reservation",
        output_directory=str(tmp_path),
        make_dir=False,
        solvent_mode="MOLECULE_ONLY",
    )
    system.type = "balls"
    command = Command(sys=system)
    command.fail_if_exists = True
    command.base_file = str(tmp_path / "synthetic-surrounding.pdb")
    command._arguments = [
        command.base_file,
        "-g", "a", "0",
        "-i",
        "-s", "nt", "pow", "and", "mv", "40",
        "--profile", "geometry",
        "--boundary", "none",
    ]
    command.parse_commands()
    command.run_exports = lambda: None
    command._save_archive = lambda: None

    command._process_system()

    interface = system.ifaces[0]
    expected = canonical_interfaces_directory(system) / slug(interface.interface_id)
    assert Path(interface.dir).resolve() == expected.resolve()
    assert Path(interface.dir).resolve() in {
        Path(path).resolve()
        for path in system._cli_reserved_interface_directories
    }


def test_full_solvent_direct_cli_keeps_water_and_uses_partial_results_scope(
    tmp_path, capsys, monkeypatch,
):
    system = System(
        spheres=pd.DataFrame({
            "loc": [
                (0.0, 0.0, 0.0), (2.0, 0.0, 0.0),
                (0.0, 2.0, 0.0), (0.0, 0.0, 2.0),
                (0.8, 0.8, 0.8),
            ],
            "rad": [0.75, 0.5, 0.75, 0.75, 0.35],
            "mass": [1.0] * 5,
            "num": list(range(5)),
            "res_name": ["ALA", "ALA", "ALA", "ALA", "HOH"],
            "name": ["CA"] * 5,
        }),
        name="synthetic-wet-interface-cli",
        output_directory=str(tmp_path),
        make_dir=False,
        solvent_mode="FULL_SOLVENT",
    )
    command = Command(sys=system)
    command.base_file = str(tmp_path / "synthetic-wet.pdb")
    command._arguments = [
        command.base_file,
        "-g", "a", "0",
        "-g", "a", "1",
        "-i",
        "-s", "nt", "pow", "--profile", "geometry",
        "--boundary", "none",
    ]
    command.parse_commands()
    command.run_exports = lambda: None
    command._save_archive = lambda: None

    original_make_interfaces = system.make_interfaces

    def make_interfaces(*args, **kwargs):
        original_make_interfaces(*args, **kwargs)
        # Keep this regression independent of the optional GUDHI alpha plugin.
        from vorpy.src.interface.geometry_analysis import InterfaceGeometryAnalysis

        system.ifaces[0].geometry_analysis = InterfaceGeometryAnalysis(
            metadata={"result_status": "PARTIAL"}, voronoi_side_1={}
        )

    monkeypatch.setattr(system, "make_interfaces", make_interfaces)

    command._process_system()

    network = system.ifaces[0].net
    assert network.network_mode == "interface"
    assert network.generator_system_indices == (0, 1, 2, 3, 4)
    assert network.solvent_local_indices == (4,)
    assert network.solvent_mode == "FULL_SOLVENT"
    assert command.interface_results
    result = command.interface_results[0]
    assert result.identity.solvent_mode == "SOLVENT_CONTEXT"
    assert result.provenance.details["source_solvent_provenance"][
        "solvent_mode"
    ] == "FULL_SOLVENT"
    assert result.provenance.details["results_solvent_provenance"][
        "scope"
    ] == "direct_interface"
    assert "Solving POW interface network..." in capsys.readouterr().out
