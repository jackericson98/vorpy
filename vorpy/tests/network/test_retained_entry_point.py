from pathlib import Path

import pytest

from vorpy.src.boundary import SolventMode, prepare_network_geometry
from vorpy.src.system import System


EDTA = Path(__file__).resolve().parents[2] / "data" / "EDTA.pdb"


def test_system_api_retains_edta_mg_without_readding_water():
    system = System(
        file=str(EDTA),
        make_dir=False,
        retained_system_indices=(32,),
    )

    geometry = prepare_network_geometry(system, SolventMode.MOLECULE_ONLY)

    assert system.retained_system_indices == (32,)
    assert geometry.system_indices == tuple(range(33))
    assert geometry.retained_system_indices == (32,)
    assert geometry.retained_generator_ids == ("system:32",)
    assert geometry.solvent_local_indices == ()


@pytest.mark.parametrize(
    "value, error",
    [
        ("32", ValueError),
        ((-1,), IndexError),
    ],
)
def test_system_api_rejects_invalid_retained_generator_arguments(value, error):
    with pytest.raises(error):
        System(spheres=[], make_dir=False, retained_system_indices=value)


def test_prepare_geometry_rejects_retained_id_outside_loaded_system():
    system = System(file=str(EDTA), make_dir=False, retained_system_indices=(999,))

    with pytest.raises(IndexError, match="retained System index 999"):
        prepare_network_geometry(system, SolventMode.MOLECULE_ONLY)
