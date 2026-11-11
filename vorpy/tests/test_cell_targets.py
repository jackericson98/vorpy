from types import SimpleNamespace

import pandas as pd
import pytest

from vorpy.src.boundary import (
    CellTarget,
    SolventMode,
    molecular_target_system_indices,
    normalize_cell_target,
)
from vorpy.src.command.vpy_cmnd import Command, cell_target_cli_options
from vorpy.src.group import Group
from vorpy.src.system import System


def _synthetic_system():
    """Return a molecular system with two molecule and two water atoms."""
    rows = pd.DataFrame(
        {
            "loc": [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0),
                    (0.0, 2.0, 0.0), (2.0, 2.0, 0.0)],
            "rad": [1.5, 1.5, 1.4, 1.4],
            "res_name": ["MOL", "MOL", "SOL", "SOL"],
        }
    )
    molecule = SimpleNamespace(atoms=(0, 1), name="MOL", seq=1)
    return SimpleNamespace(
        type="molecular",
        balls=rows,
        residues=(molecule,),
    )


def test_cell_target_normalization_and_synthetic_competitor_contract():
    system = _synthetic_system()

    assert normalize_cell_target(None) is CellTarget.MOLECULAR
    assert normalize_cell_target("all-cells") is CellTarget.ALL
    assert molecular_target_system_indices(system, range(4)) == (0, 1)


def test_group_keeps_water_competitors_but_targets_molecular_cells(tmp_path):
    system = System.__new__(System)
    system.type = "molecular"
    system.balls = _synthetic_system().balls
    system.residues = list(_synthetic_system().residues)
    system.chains = []
    system.groups = []
    system.files = {"dir": str(tmp_path)}
    system.name = "synthetic"
    system.foam_box = None
    system.round_to = 6
    system.solvent_mode = SolventMode.SOLVENT_CONTEXT
    system.retained_system_indices = ()
    system.boundary_generators = []
    system.boundary_config = None

    molecular = Group(
        system,
        name="molecular",
        atoms=[0, 1],
        solvent_mode=SolventMode.SOLVENT_CONTEXT,
        make_net=True,
        build_net=False,
        print_metrics=False,
    )
    assert len(molecular.net.balls) == 4
    assert tuple(molecular.net.group) == (0, 1)
    assert molecular.net.solvent_local_indices == (2, 3)

    all_cells = Group(
        system,
        name="all",
        atoms=list(range(4)),
        solvent_mode=SolventMode.SOLVENT_CONTEXT,
        cell_target=CellTarget.ALL,
        make_net=True,
        build_net=False,
        print_metrics=False,
    )
    assert tuple(all_cells.net.group) == (0, 1, 2, 3)

    water = Group(
        system,
        name="water",
        atoms=[2, 3],
        solvent_mode=SolventMode.SOLVENT_CONTEXT,
        make_net=True,
        build_net=False,
        print_metrics=False,
    )
    assert tuple(water.net.group) == (2, 3)


@pytest.mark.parametrize("net_type", ["aw", "pow", "del"])
def test_target_contract_is_method_independent(tmp_path, net_type):
    system = System.__new__(System)
    source = _synthetic_system()
    system.type = source.type
    system.balls = source.balls
    system.residues = list(source.residues)
    system.chains = []
    system.groups = []
    system.files = {"dir": str(tmp_path / net_type)}
    system.name = f"synthetic_{net_type}"
    system.foam_box = None
    system.round_to = 6
    system.solvent_mode = SolventMode.SOLVENT_CONTEXT
    system.retained_system_indices = ()
    system.boundary_generators = []
    system.boundary_config = None
    group = Group(
        system,
        name="molecular",
        atoms=[0, 1],
        settings={"net_type": net_type},
        solvent_mode=SolventMode.SOLVENT_CONTEXT,
        make_net=True,
        build_net=False,
        print_metrics=False,
    )
    assert tuple(group.net.group) == (0, 1)
    assert len(group.net.balls) == 4


def test_cli_cell_target_is_independent_of_solvent_mode():
    remaining, target = cell_target_cli_options(
        ["input.pdb", "--solvent-context", "--cell-target", "all"]
    )
    assert remaining == ["input.pdb", "--solvent-context"]
    assert target is CellTarget.ALL


def test_edta_default_policy_has_32_targets_and_612_competitors():
    system = System(
        file="vorpy/data/EDTA.pdb",
        make_dir=False,
        print_actions=False,
        solvent_mode=SolventMode.SOLVENT_CONTEXT,
    )
    command = Command(sys=system)
    command.base_file = "vorpy/data/EDTA.pdb"
    command.solvent_mode = SolventMode.SOLVENT_CONTEXT
    command._resolve_generator_policy()

    policy = command._source_generator_policy
    assert policy["active_generator_count"] == 612
    assert policy["target_cell_count"] == 32
    assert len(command._cell_target_system_indices(command._source_system_indices)) == 32
    assert policy["target_system_indices"] == tuple(range(32))
