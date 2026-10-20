from types import SimpleNamespace

import pandas as pd
import pytest

from vorpy.src.geometry.duals import build_dual
from vorpy.src.geometry.filtrations.alpha import AlphaFiltration, AlphaSimplex
from vorpy.src.geometry.interfaces import build_alpha_interface, build_full_interface


def _network(scheme):
    return SimpleNamespace(
        settings={"net_type": scheme},
        balls=pd.DataFrame(
            {
                "loc": [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 2.0, 0.0)],
                "rad": [1.0, 1.0, 1.0],
                "complete": [True, False, True],
            },
            index=[10, 13, 16],
        ),
        surfs=pd.DataFrame(
            {
                "balls": [[10, 13], [10, 16], [13, 16]],
                "edges": [[7], [7], [7]],
                "verts": [[4, 5], [4, 5], [4, 5]],
                "sa": [3.0, 4.0, 5.0],
            },
            index=[20, 21, 22],
        ),
        edges=pd.DataFrame(
            {"balls": [[10, 13, 16]], "verts": [[4, 5]], "surfs": [[20, 21, 22]]},
            index=[7],
        ),
        verts=pd.DataFrame(
            {
                "balls": [[10, 13, 16], [10, 13, 16]],
                "loc": [(1.0, 1.0, 0.0), (1.0, 1.0, 1.0)],
            },
            index=[4, 5],
        ),
    )


def _alpha(network):
    records = {dimension: {} for dimension in range(4)}
    for generator_id in network.balls.index:
        row = AlphaSimplex(0, (int(generator_id),), 0.0, 0.0, "synthetic", True)
        records[0][row.generator_tuple] = row
    row = AlphaSimplex(1, (10, 13), 0.0, 0.0, "synthetic", True)
    records[1][row.generator_tuple] = row
    return AlphaFiltration(
        network,
        network.settings["net_type"],
        records,
        {"native_alpha_units": "A^2", "filtration_kind": "synthetic"},
    )


@pytest.mark.parametrize("scheme", ["prm", "pow", "aw"])
def test_optional_dual_and_alpha_views_do_not_redefine_primary_network(scheme):
    network = _network(scheme)
    before = {
        name: frame.copy(deep=True)
        for name, frame in (
            ("balls", network.balls),
            ("surfs", network.surfs),
            ("edges", network.edges),
            ("verts", network.verts),
        )
    }
    counts = {name: len(frame) for name, frame in before.items()}
    completeness = tuple(bool(value) for value in before["balls"]["complete"])

    dual = build_dual(network)
    full = build_full_interface(network, dual, {10}, {13, 16})
    assert dual.network is network
    assert full.metadata["physical_geometry_source"] == (
        "network.surfs/network.edges/network.verts"
    )
    if scheme == "pow":
        selected = build_alpha_interface(
            network,
            dual,
            _alpha(network),
            0.0,
            {10},
            {13, 16},
            calculate_curvature=False,
        )
        assert selected.selection_mode == "alpha"
        assert len(selected.selected_surfaces) == 1
        assert len(network.surfs) == counts["surfs"]

    assert {name: len(getattr(network, name)) for name in counts} == counts
    assert tuple(bool(value) for value in network.balls["complete"]) == completeness
    for name, expected in before.items():
        pd.testing.assert_frame_equal(getattr(network, name), expected)
