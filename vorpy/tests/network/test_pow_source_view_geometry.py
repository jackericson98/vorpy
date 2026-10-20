import contextlib
import io

import pytest

from vorpy.src.network.network import Network
from vorpy.src.network.views import NetworkView
from vorpy.tests.network.test_views import _network


LOCATIONS = [
    (0.0, 0.0, 0.0),
    (2.0, 0.0, 0.0),
    (0.0, 2.0, 0.0),
    (0.0, 0.0, 2.0),
    (0.8, 0.8, 0.8),
]
RADII = [0.75, 0.5, 0.75, 0.75, 0.35]


def _build_reference(scheme, *, interface, competitor_near, timing=False):
    locations = list(LOCATIONS)
    if competitor_near:
        locations[1] = (0.2, 0.8, 0.2)
    network = Network(
        locs=locations,
        rads=RADII,
        group=[0, 4] if interface else list(range(5)),
        iface_grps=({0}, {4}) if interface else None,
        settings=None,
        sort_balls=False,
    )
    network.default_settings(
        surf_res=0.25,
        box_size=1.5,
        max_vert=40.0,
        build_type="all",
        net_type=scheme,
        num_splits=2,
    )
    network.settings["timing"] = timing
    network.sort_balls()
    with contextlib.redirect_stdout(io.StringIO()):
        network.build(calculate_curvature=False)
    return network


def _interface_feature_signature(network, group_a=(0,), group_b=(4,)):
    """Return physical feature identities, incidence, and area for A/B facets."""
    stable_ids = {
        int(row["num"]): f"system:{int(row['num'])}"
        for _, row in network.balls.iterrows()
    }

    def generators(row):
        return tuple(sorted(stable_ids[int(value)] for value in row["balls"]))

    surfaces = {
        int(index): row
        for index, row in network.surfs.iterrows()
        if set(group_a).intersection(row["balls"])
        and set(group_b).intersection(row["balls"])
    }
    edge_ids = sorted({
        int(edge_id)
        for row in surfaces.values()
        for edge_id in row["edges"]
    })
    vertex_ids = sorted({
        int(vertex_id)
        for row in surfaces.values()
        for vertex_id in row["verts"]
    })
    return {
        "pairs": tuple(sorted(
            (
                next(iter(set(group_a).intersection(row["balls"]))),
                next(iter(set(group_b).intersection(row["balls"]))),
            )
            for row in surfaces.values()
        )),
        "surfaces": tuple(sorted(
            (
                generators(row),
                tuple(int(value) for value in row["verts"]),
                tuple(int(value) for value in row["edges"]),
                float(row["sa"]),
            )
            for row in surfaces.values()
        )),
        "edges": tuple(sorted(
            (
                generators(network.edges.iloc[index]),
                tuple(int(value) for value in network.edges.iloc[index]["verts"]),
                tuple(int(value) for value in network.edges.iloc[index]["surfs"]),
            )
            for index in edge_ids
        )),
        "vertices": tuple(sorted(
            (
                generators(network.verts.iloc[index]),
                tuple(int(value) for value in network.verts.iloc[index]["edges"]),
                tuple(int(value) for value in network.verts.iloc[index]["surfs"]),
            )
            for index in vertex_ids
        )),
        "area": sum(float(row["sa"]) for row in surfaces.values()),
        "surface_count": len(surfaces),
        "edge_count": len(edge_ids),
        "vertex_count": len(vertex_ids),
    }


@pytest.mark.parametrize("scheme", ["aw", "pow", "prm"])
def test_interface_geometry_matches_full_reference_and_competitor_suppression(scheme):
    full = _build_reference(scheme, interface=False, competitor_near=False)
    selected = _build_reference(scheme, interface=True, competitor_near=False)

    full_faces = {
        tuple(sorted(int(ball) for ball in row["balls"])): float(row["sa"])
        for _, row in full.surfs.iterrows()
        if 0 in row["balls"] and 4 in row["balls"]
    }
    selected_faces = {
        tuple(sorted(int(ball) for ball in row["balls"])): float(row["sa"])
        for _, row in selected.surfs.iterrows()
        if 0 in row["balls"] and 4 in row["balls"]
    }

    assert selected_faces == full_faces
    assert len(full_faces) == 1
    assert len(selected.verts) == len(full.verts)
    assert len(selected.edges) == len(full.edges)

    full_features = _interface_feature_signature(full)
    selected_features = _interface_feature_signature(selected)
    assert selected_features == full_features
    assert selected_features["pairs"] == ((0, 4),)
    assert selected_features["area"] == pytest.approx(sum(selected_faces.values()))

    # Both builds expose explicit evidence; the bounded no-boundary fixture is
    # intentionally not promoted to complete-cell certification.
    assert full.balls["complete"].tolist() == selected.balls["complete"].tolist()
    assert not bool(full.balls["complete"].all())

    full_blocked = _build_reference(scheme, interface=False, competitor_near=True)
    selected_blocked = _build_reference(scheme, interface=True, competitor_near=True)

    assert not any(
        0 in row["balls"] and 4 in row["balls"]
        for _, row in full_blocked.surfs.iterrows()
    )
    assert not any(
        0 in row["balls"] and 4 in row["balls"]
        for _, row in selected_blocked.surfs.iterrows()
    )

    # Generator 1 is a third-party competitor, not either selected owner.
    assert 1 not in {0, 4}


@pytest.mark.parametrize("scheme", ["aw", "pow", "prm"])
def test_interface_work_evidence_is_explicit_and_not_claimed_selective(scheme):
    full = _build_reference(
        scheme, interface=False, competitor_near=False, timing=True,
    )
    selected = _build_reference(
        scheme, interface=True, competitor_near=False, timing=True,
    )

    full_counts = full.build_timing["counts"]
    selected_counts = selected.build_timing["counts"]
    assert full_counts["balls"] == selected_counts["balls"] == 5
    assert full_counts["vertices"] == len(full.verts)
    assert selected_counts["vertices"] == len(selected.verts)
    assert selected_counts["surfaces"] == 1
    assert selected_counts["edges"] == 3
    assert selected_counts["candidate_facets"] == full_counts["candidate_facets"]
    assert selected.group == [0, 4]
    assert len(selected.iface_grps) == 2


def test_pow_source_view_ignores_unrelated_incomplete_exterior_cell():
    source = _network(scheme="pow")
    source.solvent_mode = "FULL_SOLVENT"
    source.settings["solvent_mode"] = "FULL_SOLVENT"
    source.balls.loc[2, "stable_id"] = "water-1"
    source.balls.loc[3, "stable_id"] = "water-2"
    source.balls.loc[3, "complete"] = False
    source.surfs.at[0, "balls"] = [0, 1, 2]
    source.surfs.at[0, "sa"] = 2.5

    view = NetworkView.for_interface(source, [10], [11])

    assert view.source_ids("surface") == (0,)
    assert view.completeness.proven
    assert view.surfaces[0].generator_system_ids == (10, 11, 12)
    assert view.surfaces[0].data["sa"] == 2.5
    assert view.surfaces[0].evidence.supported is True
    assert view.surfaces[0].evidence.complete is True
    assert view.surfaces[0].evidence.bounded is True
    assert not view.completeness.unresolved_features
    assert 3 not in view.support_generator_local_indices
