from types import SimpleNamespace

import pandas as pd
import pytest

from vorpy.src.boundary import BoundaryConfig, BoundaryGenerator
from vorpy.src.network.network import Network
from vorpy.src.network.views import NetworkView, NetworkViewError


def _network(
    *,
    scheme="prm",
    incomplete=False,
    broken=False,
    curved=False,
    source_network_id="producer:test-network",
):
    locations = [
        (0.0, 0.0, 0.0),
        (2.0, 0.0, 0.0),
        (0.0, 2.0, 0.0),
        (0.0, 0.0, 2.0),
        (4.0, 4.0, 4.0),
    ]
    balls = pd.DataFrame({
        "loc": locations,
        "rad": [1.0] * 5,
        "num": [0, 1, 2, 3, 4],
        "system_num": [10, 11, 12, 13, None],
        "stable_id": ["atom-A", "atom-B", "atom-C", "atom-D", None],
        "is_boundary_generator": [False, False, False, False, True],
        "complete": [True, True, True, True, not incomplete],
    })
    verts = pd.DataFrame({
        "balls": [[0, 1, 4], [0, 1, 4], [0, 1, 4], [0, 1, 4]],
        "loc": [
            (0.2, 0.2, 0.2), (0.8, 0.2, 0.2),
            (0.8, 0.8, 0.2), (0.2, 0.8, 0.2),
        ],
        "rad": [0.2] * 4,
        "complete": [True] * 4,
    })
    edges = pd.DataFrame({
        "balls": [[0, 1, 4]] * 4,
        "verts": [[0, 1], [1, 2], [2, 3], [3, 0]],
        "surfs": [[0, 1]] * 4,
        "points": [[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]] * 4,
        "complete": [True] * 4,
    })
    surface_points = (
        [(0.2, 0.2, 0.2), (0.8, 0.2, 0.3), (0.8, 0.8, 0.2), (0.2, 0.8, 0.1)]
        if curved
        else [(0.2, 0.2, 0.2), (0.8, 0.2, 0.2), (0.8, 0.8, 0.2), (0.2, 0.8, 0.2)]
    )
    surfaces = pd.DataFrame({
        "balls": [[0, 1, 4], [0, 2]],
        "verts": [[0, 1, 2, 3], [0, 1, 2, 3]],
        "edges": [[0, 1, 2, 3], [0, 1, 2, 3]],
        "points": [surface_points, surface_points],
        "tris": [[(0, 1, 2)], [(0, 1, 2)]],
        "orientation": [(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0)],
        "mean_curvature": [0.5 if curved else 0.0] * 2,
        "gauss_curvature": [0.25 if curved else 0.0] * 2,
        "complete": [True, True],
    })
    if broken:
        surfaces.at[0, "edges"] = [0, 1, 2, 99]

    network = Network(
        locs=locations,
        rads=[1.0] * 5,
        group=[10, 11, 12, 13],
        settings={
            "net_type": scheme,
            "surf_res": 0.25,
            "max_vert": 5.0,
            "box_size": 1.5,
            "num_splits": 2,
            "build_type": "all",
        },
        balls=balls,
        verts=verts,
        edges=edges,
        surfs=surfaces,
        system=SimpleNamespace(frame_index=7),
        source_network_id=(
            "test-source-network"
            if source_network_id is None else source_network_id
        ),
    )
    network.boundary_config = BoundaryConfig(mode="shell")
    network.boundary_mode = "shell"
    network.boundary_indices = (4,)
    network.boundary_generators = (
        BoundaryGenerator(position=locations[4], radius=1.0),
    )
    return network


@pytest.mark.parametrize("scheme,curved", [("prm", False), ("aw", True)])
def test_group_view_preserves_identity_incidence_orientation_and_scheme(scheme, curved):
    network = _network(scheme=scheme, curved=curved)
    view = NetworkView.for_group(network, [10])

    assert view.scheme == scheme
    assert view.frame_id == "7"
    assert view.system_ids == (10,)
    assert view.local_to_system == ((0, 10), (1, 11), (2, 12), (3, 13), (4, None))
    assert view.cells[0].source_id == 0
    assert view.cells[0].generator_system_ids == (10,)
    assert view.generator_identity == ((0, "atom-A", 10),)
    assert view.surfaces[0].source_id == 0
    assert view.surfaces[0].generator_ids == ("atom-A", "atom-B", "boundary:4")
    assert view.surfaces[0].generator_system_ids == (10, 11, None)
    assert view.surfaces[0].data["orientation"] == (1.0, 0.0, 0.0)
    assert view.surfaces[0].incidence["edges"] == (0, 1, 2, 3)
    assert view.surfaces[0].data["mean_curvature"] == (0.5 if curved else 0.0)
    assert view.boundary_generator_ids == ("boundary:4",)
    assert view.completeness.proven


def test_interface_view_preserves_ordered_sides_and_reversed_identity():
    network = _network()
    forward = network.make_view(iface_grps=([10], [11]))
    reverse = network.make_view(iface_grps=([11], [10]))

    assert forward.selection_kind == "interface"
    assert forward.group_a_ids == (10,)
    assert forward.group_b_ids == (11,)
    assert reverse.group_a_ids == (11,)
    assert reverse.group_b_ids == (10,)
    assert forward.selection_id != reverse.selection_id
    assert forward.source_ids("surface") == reverse.source_ids("surface") == (0,)
    assert forward.surfaces[0].data["orientation"] == (1.0, 0.0, 0.0)


def test_view_is_snapshot_and_never_calls_network_build_stages(monkeypatch):
    network = _network()
    balls_before = network.balls.copy(deep=True)
    surfs_before = network.surfs.copy(deep=True)
    group_before = list(network.group)
    iface_before = network.iface_grps

    def fail(*args, **kwargs):
        raise AssertionError("view creation invoked a network build stage")

    for name in ("build", "connect", "build_edges", "build_surfaces", "find_verts"):
        monkeypatch.setattr(Network, name, fail)

    view = NetworkView.for_interface(network, [10], [11])
    pd.testing.assert_frame_equal(network.balls, balls_before)
    pd.testing.assert_frame_equal(network.surfs, surfs_before)
    assert list(network.group) == group_before
    assert network.iface_grps == iface_before

    network.balls.loc[0, "rad"] = 99.0
    network.surfs.at[0, "orientation"] = (0.0, 0.0, 1.0)

    assert view.cells[0].data["rad"] == 1.0
    assert view.surfaces[0].data["orientation"] == (1.0, 0.0, 0.0)
    assert not hasattr(view, "materialize")
    assert not surfs_before.empty


def test_selection_identity_is_deterministic_for_equal_source_content():
    first = NetworkView.for_group(_network(), [11, 10])
    second = NetworkView.for_group(_network(), [10, 11])

    assert first.selection_id == second.selection_id
    assert first.source_network_id == second.source_network_id
    assert not first.source_network_id_is_provisional
    assert first.source_network_id_is_authoritative
    assert first.solve_universe_id == second.solve_universe_id


def test_view_explicit_source_id_without_archive_is_authoritative():
    network = _network(source_network_id="producer-source")

    view = NetworkView.for_group(network, [10])

    assert view.source_network_id == "producer-source"
    assert view.source_network_id_is_authoritative
    assert not view.source_network_id_is_provisional


def test_view_explicit_source_id_remains_authoritative_when_archived():
    network = _network(source_network_id="producer-source")
    network.archive_id = "archive-artifact"

    view = NetworkView.for_group(network, [10])

    assert view.source_network_id == "producer-source"
    assert view.source_network_id_is_authoritative
    assert view.source_network_id != network.archive_id


def test_view_missing_source_id_without_archive_is_provisional():
    network = _network()
    network.source_network_id = None

    view = NetworkView.for_group(network, [10])

    assert view.source_network_id == f"network:{view.solve_context.digest}"
    assert view.source_network_id_is_provisional
    assert not view.source_network_id_is_authoritative


def test_view_missing_source_id_ignores_archive_identity():
    network = _network()
    network.source_network_id = None
    network.archive_id = "archive-only-identity"

    view = NetworkView.for_group(network, [10])

    assert view.source_network_id == f"network:{view.solve_context.digest}"
    assert view.source_network_id != network.archive_id
    assert view.source_network_id_is_provisional
    assert not view.source_network_id_is_authoritative


def test_view_rejects_conflicting_explicit_source_identity():
    network = _network(source_network_id="producer-source-a")

    with pytest.raises(NetworkViewError, match="conflicts with the source Network"):
        NetworkView.for_group(
            network,
            [10],
            source_network_id="producer-source-b",
        )


def test_view_preserves_explicit_source_network_namespace():
    network = _network(source_network_id="source-network-7")

    view = network.make_view(group_ids=[10])
    repeated = network.make_view(group_ids=[10])
    network.source_network_id = None
    overridden = NetworkView.for_group(
        network,
        [10],
        source_network_id="source-network-8",
    )

    assert view.source_network_id == "source-network-7"
    assert view.network_id == view.source_network_id
    assert not view.source_network_id_is_provisional
    assert repeated.selection_id == view.selection_id
    assert overridden.source_network_id == "source-network-8"
    assert overridden.selection_id != view.selection_id


def test_incomplete_or_broken_selection_fails_closed():
    with pytest.raises(NetworkViewError, match="completeness"):
        NetworkView.for_interface(_network(incomplete=True), [10], [11])

    with pytest.raises(NetworkViewError, match="incidence"):
        NetworkView.for_interface(_network(broken=True), [10], [11])

    with pytest.raises(NetworkViewError, match="unknown"):
        NetworkView.for_group(_network(), [999])

    with pytest.raises(NetworkViewError, match="Unsupported source solve context"):
        NetworkView.for_group(_network(scheme="unsupported"), [10])
