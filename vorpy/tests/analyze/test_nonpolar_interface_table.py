from types import SimpleNamespace

import pandas as pd
import pytest

from vorpy.src.analyze.nonpolar_interface import build_geometry_tables


def network_fixture():
    balls = pd.DataFrame([
        {"num": 0, "name": "CA", "res_name": "ALA", "res_seq": 1, "chain_name": "A",
         "element": "C", "rad": 1.7, "vol": 4., "sa": 10., "complete": True},
        {"num": 1, "name": "CB", "res_name": "LEU", "res_seq": 2, "chain_name": "A",
         "element": "C", "rad": 1.7, "vol": 5., "sa": 11., "complete": True},
        {"num": 2, "name": "O", "res_name": "SER", "res_seq": 3, "chain_name": "B",
         "element": "O", "rad": 1.4, "vol": 6., "sa": 12., "complete": True},
    ])
    surfs = pd.DataFrame([
        {"balls": [0, 1], "sa": 2., "mean_curv": .1, "avg_mean_curv": .2,
         "gauss_curv": .03, "avg_gauss_curv": .04, "int_mean_curv": .2,
         "int_mean_curv_sq": .02, "int_gauss_curv": .06},
        {"balls": [0, 2], "sa": 3., "mean_curv": .3, "avg_mean_curv": .4,
         "gauss_curv": .05, "avg_gauss_curv": .06, "int_mean_curv": .9,
         "int_mean_curv_sq": .27, "int_gauss_curv": .15},
    ])
    return SimpleNamespace(
        balls=balls,
        surfs=surfs,
        iface_grps=({0, 1}, {2}),
        group_name="fixture",
        sys=SimpleNamespace(balls=balls),
    )


def test_tables_label_cross_side_surface_and_aggregate_area():
    result = build_geometry_tables(network_fixture(), system_id="complex", frame_id=4)
    surfaces = result["surfaces"].set_index("surface_id")
    atoms = result["atoms"].set_index("atom_index")

    assert surfaces.loc[0, "interface_label"] == 0
    assert surfaces.loc[1, "interface_label"] == 1
    assert surfaces.loc[1, "nonpolar_surface_label"] == "mixed"
    assert atoms.loc[0, "interface_area"] == pytest.approx(3.)
    assert atoms.loc[0, "interface_fraction"] == pytest.approx(.3)
    assert atoms.loc[2, "interface_label"] == 1
    assert atoms.loc[2, "nonpolar_class"] == "other"
    assert result["metadata"]["surface_count"] == 2


def test_tables_do_not_invent_interface_label_without_groups():
    network = network_fixture()
    network.iface_grps = None
    result = build_geometry_tables(network)
    assert result["surfaces"]["interface_label"].isna().all()
    assert result["atoms"]["interface_label"].isna().all()
    assert result["metadata"]["interface_groups_supplied"] is False


def test_strict_rule_excludes_water_and_ions():
    network = network_fixture()
    network.balls.loc[0, ["res_name", "element"]] = ["SOL", "O"]
    result = build_geometry_tables(network)
    atom = result["atoms"].set_index("atom_index").loc[0]
    assert atom["nonpolar_class"] == "excluded"


def test_unknown_rule_fails_explicitly():
    with pytest.raises(ValueError, match="Unknown nonpolar"):
        build_geometry_tables(network_fixture(), nonpolar_rule="guess")
