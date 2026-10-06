from types import SimpleNamespace

import pandas as pd
import pytest

from vorpy.src.geometry.filtrations import build_alpha_filtration


def _network(scheme, locations, radii, surfaces=()):
    return SimpleNamespace(
        settings={"net_type": scheme},
        group_name="alpha-adapter-test",
        balls=pd.DataFrame(
            {"loc": locations, "rad": radii, "complete": [True] * len(radii)},
            index=[10 + 3 * index for index in range(len(radii))],
        ),
        surfs=pd.DataFrame(surfaces, columns=["balls", "edges", "verts"]),
        edges=pd.DataFrame(columns=["balls", "verts", "points"]),
        verts=pd.DataFrame(columns=["balls", "loc"]),
    )


def test_primitive_adapter_preserves_squared_radius_native_values_and_ids():
    network = _network("prm", [(0, 0, 0), (2, 0, 0)], [1, 1])
    filtration = build_alpha_filtration(network)
    assert filtration.metadata["native_alpha_units"] == "A^2"
    assert filtration.simplex_birth((10,)).filtration_birth == pytest.approx(0.0)
    assert filtration.simplex_birth((13,)).filtration_birth == pytest.approx(0.0)
    edge = filtration.simplex_birth((13, 10))
    assert edge.generator_tuple == (10, 13)
    assert edge.filtration_birth == pytest.approx(1.0)
    assert filtration.simplices_at(1.0).counts == {0: 2, 1: 1, 2: 0, 3: 0}


def test_power_adapter_reuses_weighted_backend_and_records_exact_weights():
    network = _network("pow", [(0, 0, 0), (4, 0, 0)], [1, 2])
    filtration = build_alpha_filtration(network)
    assert filtration.metadata["weights_A2"] == [1.0, 4.0]
    assert filtration.metadata["native_alpha_units"] == "A^2"
    assert filtration.metadata["physical_offset_A"] is None
    assert filtration.simplex_birth((10,)).birth_source == "gudhi_weighted_alpha"
    edge = filtration.simplex_birth((10, 13))
    assert edge is not None
    assert edge.filtration_birth is not None
    assert filtration.interface_at(edge.filtration_birth, {10}, {13}) == [edge]


def test_aw_adapter_delegates_births_and_retains_unresolved_records():
    network = _network(
        "aw", [(0, 0, 0), (8, 0, 0)], [2, 2],
        surfaces=[([10, 13], [], [])],
    )
    filtration = build_alpha_filtration(network)
    assert filtration.metadata["experimental"] is True
    assert filtration.metadata["native_alpha_units"] == "A"
    assert filtration.metadata["physical_offset_A"] == "native alpha"
    assert filtration.simplex_birth((10,)).filtration_birth == pytest.approx(-2.0)
    pair = filtration.simplex_birth((10, 13))
    assert pair is not None and pair.filtration_birth is None
    assert pair.supported is False
    assert filtration.interface_at(100.0, {10}, {13}) == []
    assert filtration.interface_at(100.0, {10}, {13}) == []


def test_aw_unresolved_pair_is_not_fabricated_into_alpha_subcomplex():
    network = _network(
        "aw", [(0, 0, 0), (8, 0, 0)], [2, 2],
        surfaces=[([10, 13], [], [])],
    )
    filtration = build_alpha_filtration(network)
    subcomplex = filtration.simplices_at(100.0)
    assert subcomplex.counts == {0: 2, 1: 0, 2: 0, 3: 0}
    assert len(subcomplex.blocked) == 1
