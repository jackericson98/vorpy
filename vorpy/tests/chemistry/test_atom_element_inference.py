from pathlib import Path

from vorpy.src.chemistry import element_radii
from vorpy.src.objects.atom import get_element, get_radius
from vorpy.src.system import System


DATA = Path(__file__).resolve().parents[2] / "data"


def test_missing_pdb_element_column_preserves_explicit_selenium_atom_name():
    atom = {
        'name': 'SE21',
        'res_name': 'MOL',
        'res': None,
        'element': None,
        'rad': None,
    }

    atom['element'] = get_element(atom)
    assert atom['element'] == 'SE'
    assert get_radius(atom) == element_radii['SE'] == 1.90


def test_db1976_blank_element_column_uses_se_radius():
    system = System(
        file=str(DATA / "1bna_DB1976.pdb"),
        make_dir=False,
        print_actions=False,
    )

    selenium = system.balls[
        (system.balls["name"].astype(str).str.strip().str.upper() == "SE21")
    ]
    assert len(selenium) == 1
    assert selenium.iloc[0]["element"] == "SE"
    assert selenium.iloc[0]["rad"] == element_radii["SE"] == 1.90
