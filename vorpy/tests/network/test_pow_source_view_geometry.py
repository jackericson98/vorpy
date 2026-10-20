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


def _build_reference(scheme, *, interface, competitor_near):
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
    network.sort_balls()
    with contextlib.redirect_stdout(io.StringIO()):
        network.build(calculate_curvature=False)
    return network


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
