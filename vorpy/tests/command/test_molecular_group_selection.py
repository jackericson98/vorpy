from pathlib import Path

import pytest

from vorpy.src.boundary import SolventMode, prepare_network_geometry
from vorpy.src.command.group import validate_group_selection
from vorpy.src.command.vpy_cmnd import Command
from vorpy.src.system.system import System


DATA = Path(__file__).resolve().parents[2] / "data"


@pytest.mark.parametrize(
    ("filename", "argv_groups", "expected_groups", "mode", "active", "required"),
    [
        (
            "2KAI_no_sol.pdb",
            (("c A and c B",), ("c I",)),
            (3482, 858),
            SolventMode.MOLECULE_ONLY.value,
            4340,
            (),
        ),
        (
            "1AAY_no_sol.pdb",
            (("c A",), ("c B and c C",)),
            (1422, 694),
            SolventMode.SOLVENT_CONTEXT.value,
            2119,
            ("ZN",),
        ),
    ],
)
def test_quoted_canonical_chain_groups_materialize_from_pdb(
    filename, argv_groups, expected_groups, mode, active, required,
):
    """Quoted Windows selectors must reach the real VorPy group API correctly."""
    source = DATA / filename
    system = System(
        file=str(source),
        make_dir=False,
        print_actions=False,
        solvent_mode=mode,
    )
    command = Command(sys=system)
    command.base_file = str(source)
    command._arguments = [
        str(source),
        "-g", argv_groups[0][0],
        "-g", argv_groups[1][0],
        "-i",
    ]
    command.parse_commands()
    command.create_groups()

    assert tuple(len(group.ball_ndxs) for group in system.groups) == expected_groups
    assert all(group.ball_ndxs for group in system.groups)
    assert command.groups == {
        0: [["c", "A"], ["c", "B"]]
        if filename.startswith("2KAI") else [["c", "A"]],
        1: [["c", "I"]]
        if filename.startswith("2KAI") else [["c", "B"], ["c", "C"]],
    }

    geometry = prepare_network_geometry(system, solvent_mode=mode)
    assert len(geometry.system_indices) == active
    if required:
        zn_rows = {
            index for index in range(len(system.balls))
            if str(system.balls.iloc[index]["res_name"]).strip().upper() in required
        }
        assert len(zn_rows) == 3
        assert zn_rows <= set(geometry.system_indices)


def test_group_preflight_rejects_empty_or_mismatched_groups():
    source = DATA / "2KAI_no_sol.pdb"
    system = System(file=str(source), make_dir=False, print_actions=False)
    with pytest.raises(ValueError, match="membership mismatch"):
        validate_group_selection(
            system,
            ((["c", "A"],), (["c", "I"],)),
            (3482, 858),
            solvent_mode=SolventMode.MOLECULE_ONLY.value,
            expected_active_generators=4340,
        )


def test_canonical_1aay_default_materializes_molecular_groups_and_zinc():
    source = DATA / "1AAY.pdb"
    system = System(file=str(source), make_dir=False, print_actions=False)
    command = Command(sys=system)
    command.base_file = str(source)
    command._arguments = [
        str(source), "-g", "protein", "-g", "dna", "-i", "-c", "system",
    ]
    command.parse_commands()
    command.apply_settings()
    command._resolve_default_solvent_mode()
    command.create_groups()

    assert tuple(len(group.ball_ndxs) for group in system.groups) == (1422, 694)
    geometry = prepare_network_geometry(system)
    assert len(geometry.system_indices) == 2119
    zinc_rows = {
        index for index in range(len(system.balls))
        if str(system.balls.iloc[index]["res_name"]).strip().upper() == "ZN"
    }
    assert zinc_rows == set(system.retained_system_indices)
    assert zinc_rows <= set(geometry.system_indices)
