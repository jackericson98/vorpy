from itertools import combinations
from types import SimpleNamespace

import pandas as pd
import pytest

from vorpy.src.geometry.duals import (
    AWDualBuilder,
    DualComplex,
    PowerDualBuilder,
    PrimitiveDualBuilder,
    build_dual,
)


def _closed_network(scheme):
    points = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
    ]
    balls = pd.DataFrame({
        "loc": points,
        "rad": [1.0, 1.1, 1.2, 1.3],
        "num": [10, 11, 12, 13],
        "complete": [True] * 4,
    })
    verts = pd.DataFrame({
        "balls": [list(range(4))],
        "loc": [(0.25, 0.25, 0.25)],
    })
    edges = pd.DataFrame({
        "balls": [list(triple) for triple in combinations(range(4), 3)],
        "verts": [[0]] * 4,
    })
    surfs = pd.DataFrame({
        "balls": [list(pair) for pair in combinations(range(4), 2)],
        "edges": [[0]] * 6,
        "verts": [[0]] * 6,
    })
    return SimpleNamespace(
        settings={"net_type": scheme}, balls=balls, surfs=surfs,
        edges=edges, verts=verts, group_name="synthetic",
    )


@pytest.mark.parametrize(
    ("scheme", "builder"),
    [("prm", PrimitiveDualBuilder), ("pow", PowerDualBuilder), ("aw", AWDualBuilder)],
)
def test_scheme_builder_maps_network_incidence_to_abstract_simplices(scheme, builder):
    network = _closed_network(scheme)
    dual = builder(network).build()

    assert isinstance(dual, DualComplex)
    assert dual.scheme == scheme
    assert dual.simplex_counts == {0: 4, 1: 6, 2: 4, 3: 1}
    assert dual.simplex(1, (1, 0)).primal_feature_ids == ("surf:0",)
    assert dual.simplex(2, (2, 0, 1)).primal_feature_ids == ("edge:0",)
    assert dual.simplex(3, (3, 2, 1, 0)).primal_feature_ids == ("vert:0",)
    audit = dual.audit()
    assert audit.valid
    assert audit.generator_count == 4
    assert audit.feature_counts == {"cell": 4, "surf": 6, "edge": 4, "vert": 1}
    assert audit.simplex_counts == {0: 4, 1: 6, 2: 4, 3: 1}


def test_dispatch_records_exact_power_weights_and_aw_abstract_geometry_label():
    power = build_dual(_closed_network("pow"))
    assert power.simplices[0][(0,)].generator_weights == (1.0,)
    assert power.metadata["input_radii"] == "exact Network.balls.rad values"

    aw = build_dual(_closed_network("aw"))
    assert aw.simplices[3][(0, 1, 2, 3)].metadata["dual_kind"] == "apollonius_abstract_incidence"
    assert "abstract" in aw.simplices[3][(0, 1, 2, 3)].metadata["physical_realization"]


def test_dual_extraction_does_not_mutate_solved_network():
    network = _closed_network("aw")
    before = {name: frame.copy(deep=True) for name, frame in (
        ("balls", network.balls), ("surfs", network.surfs),
        ("edges", network.edges), ("verts", network.verts),
    )}
    build_dual(network)
    for name, expected in before.items():
        pd.testing.assert_frame_equal(getattr(network, name), expected)


def test_scheme_builder_rejects_mismatched_solved_network():
    with pytest.raises(ValueError, match="requires net_type='pow'"):
        PowerDualBuilder(_closed_network("aw"))


def test_unrepresentable_primal_feature_remains_explicit_in_audit():
    network = _closed_network("aw")
    network.verts.loc[1] = {"balls": [0, 1, 2, 3, 99], "loc": (2.0, 2.0, 2.0)}
    dual = build_dual(network)
    audit = dual.audit()
    assert not audit.valid
    assert audit.unmapped_feature_count == 1
    assert any("vert:1 unresolved" in issue for issue in audit.reverse_mapping_errors)
    assert dual.unresolved_features[0].reason.startswith("Expected 4 distinct generators")
