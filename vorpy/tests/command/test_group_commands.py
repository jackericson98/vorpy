import pandas as pd

from vorpy.src.command.group import get_group_spheres


def _atoms():
    # Deliberately make 'num' 1-based like PDB serials so an off-by-one
    # regression cannot accidentally pass.
    return pd.DataFrame(
        {
            "num": [735, 736, 737, 738],
            "name": ["C1'", "H1'", "N9", "C8"],
            "res_name": ["DG", "DG", "DG", "DG"],
            "res_seq": [24, 24, 24, 24],
            "element": ["C", "H", "N", "C"],
        }
    )


def test_numeric_atom_selection_is_system_position():
    atoms = _atoms()

    assert get_group_spheres(atoms, ["0"]) == [0]
    assert get_group_spheres(atoms, ["1"]) == [1]
    assert get_group_spheres(atoms, ["2"]) == [2]

    # Critical regression: selector 2 chooses row 2 (N9), not atoms['num'] == 2.
    assert atoms.iloc[get_group_spheres(atoms, ["2"])[0]]["name"] == "N9"


def test_numeric_atom_selection_does_not_match_num_column():
    atoms = _atoms()

    # 736 exists in atoms['num'], but it is not a valid positional index in
    # this four-row test system.  It must NOT select H1'.
    assert get_group_spheres(atoms, ["736"]) == []


def test_atom_range_is_positional_and_inclusive():
    atoms = _atoms()
    assert get_group_spheres(atoms, ["1-3"]) == [1, 2, 3]


def test_atom_name_returns_system_positions():
    atoms = _atoms()
    assert get_group_spheres(atoms, ["N9"]) == [2]


def test_element_returns_system_positions():
    atoms = _atoms()
    assert get_group_spheres(atoms, ["carbon"]) == [0, 3]
