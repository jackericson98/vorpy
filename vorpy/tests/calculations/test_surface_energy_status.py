import numpy as np
import pandas as pd
from types import SimpleNamespace

from vorpy.src.calculations.surface_energy import (
    calc_helfrich_energy,
    calc_morphometric_energy,
    calc_surface_energy_geometry_from_curvatures,
    sum_surface_energy_geometry,
)
from vorpy.src.group.sort import get_info
from vorpy.tests.calculations.test_polyhedral_gaussian_curvature import _cube_network


def test_missing_curvature_does_not_become_area_only_helfrich_energy():
    geometry = sum_surface_energy_geometry([{"Area": 2.0}])

    assert geometry["Curvature Status"] == "NOT_CALCULATED"
    assert geometry["Integrated Gaussian Curvature"] is None
    energy = calc_helfrich_energy(geometry, 1.0, 2.0, 3.0)
    assert energy["Area Energy"] == 2.0
    assert energy["Bending Energy"] is None
    assert energy["Gaussian Curvature Energy"] is None
    assert energy["Total Energy"] is None
    assert energy["Status"] == "NOT_CALCULATED"

    morphometric = calc_morphometric_energy(geometry, 1.0, 2.0, 3.0)
    assert morphometric["Total Energy"] is None
    assert morphometric["Status"] == "NOT_CALCULATED"


def test_calculated_zero_curvature_remains_numeric_zero():
    geometry = {
        "Area": 2.0,
        "Integrated Mean Curvature": 0.0,
        "Integrated Mean Curvature Squared": 0.0,
        "Integrated Gaussian Curvature": 0.0,
        "Curvature Status": "CERTIFIED",
    }

    energy = calc_helfrich_energy(geometry, 1.0, 2.0, 3.0)
    assert energy["Status"] == "CERTIFIED"
    assert energy["Bending Energy"] == 0.0
    assert energy["Gaussian Curvature Energy"] == 0.0
    assert energy["Total Energy"] == 2.0


def test_invalid_triangle_curvature_is_unresolved_but_zero_is_certified():
    points = np.asarray(((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)))
    tris = ((0, 1, 2),)

    unresolved = calc_surface_energy_geometry_from_curvatures(
        points, tris, (np.nan,), (0.0,)
    )
    assert unresolved["Curvature Status"] == "UNRESOLVED"
    assert unresolved["Integrated Mean Curvature"] is None
    assert unresolved["Integrated Gaussian Curvature"] is None

    certified = calc_surface_energy_geometry_from_curvatures(
        points, tris, (0.0,), (0.0,)
    )
    assert certified["Curvature Status"] == "CERTIFIED"
    assert certified["Integrated Mean Curvature"] == 0.0
    assert certified["Integrated Gaussian Curvature"] == 0.0


def test_group_aggregation_marks_missing_gaussian_curvature_unresolved():
    network = _cube_network("aw")
    locations = [
        (-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
        (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1), (0, 0, 0),
    ]
    network.balls = pd.DataFrame({
        "loc": [np.asarray(point, dtype=float) for point in locations],
        "num": list(range(9)), "rad": [1.0] * 9,
        "complete": [True] * 9, "vol": [1.0] * 9,
        "vdw_vol": [1.0] * 9, "mass": [1.0] * 9,
        "com": [np.asarray(point, dtype=float) for point in locations],
    })
    network.surfs["balls"] = [[0, 8]] * 6
    network.surfs["sa"] = [1.0] * 6
    network.surfs["int_mean_curv"] = [0.0] * 6
    network.surfs["int_mean_curv_sq"] = [0.0] * 6
    network.surfs["int_gauss_curv"] = [None] + [0.0] * 5
    network.surfs["flat"] = [True] * 6
    network.edges["length"] = [1.0] * len(network.edges)
    group = SimpleNamespace(
        net=network,
        ball_ndxs=list(range(8)),
        layer_surfs=[list(range(6))],
        layer_net_atoms=[list(range(8))],
        get_layers=lambda max_layers=1: None,
    )

    get_info(group)

    assert group.gaussian_curvature_status == "UNRESOLVED"
    assert group.int_gauss_curv is None
