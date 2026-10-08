import time
import json
from types import SimpleNamespace

import numpy as np

from vorpy.src.network.build_net import (
    add_build_edges,
    build,
    get_build_edges,
    get_build_surfs,
)


V_BALLS = [
    [0, 1, 2, 3],
    [0, 1, 2, 4],
    [0, 1, 3, 4],
]
V_LOCS = [
    np.array([0.0, 0.0, 0.0]),
    np.array([1.0, 0.0, 0.0]),
    np.array([0.0, 1.0, 0.0]),
]
V_DUBS = [0, 0, 0]


def _net(timing):
    return SimpleNamespace(
        _timing_enabled=timing,
        settings={'verbose': False},
        update_progress=lambda *args, **kwargs: None,
    )


def _build(*, timing, interface):
    return build(
        V_BALLS,
        V_LOCS,
        V_DUBS,
        num_balls=5,
        my_time=time.time(),
        interface=interface,
        iface_grps=({0}, {1}) if interface else None,
        net=_net(timing),
    )


def test_legacy_build_timing_disabled_and_enabled_match():
    disabled = _build(timing=False, interface=False)
    enabled = _build(timing=True, interface=False)
    assert disabled == enabled


def test_interface_build_timing_disabled_and_enabled_match():
    disabled = _build(timing=False, interface=True)
    enabled = _build(timing=True, interface=True)
    assert disabled == enabled


def test_explicit_none_and_empty_helper_timings_match():
    num_balls = 5
    b_verts = [[] for _ in range(num_balls)]
    for vertex_id, balls in enumerate(V_BALLS):
        for ball_id in balls:
            b_verts[ball_id].append(vertex_id)

    edges_without_timing = get_build_edges(
        b_verts, V_BALLS, V_LOCS, V_DUBS, time.time(), timings=None,
    )
    timings = {}
    edges_with_empty_timing = get_build_edges(
        b_verts, V_BALLS, V_LOCS, V_DUBS, time.time(), timings=timings,
    )
    assert edges_without_timing == edges_with_empty_timing
    assert timings

    e_balls, e_verts = edges_without_timing
    b_edges, v_edges = add_build_edges(
        num_balls, e_balls, len(V_BALLS), e_verts,
    )
    surfaces_without_timing = get_build_surfs(
        b_verts, b_edges, V_BALLS, v_edges, e_balls, time.time(),
        interface=True, iface_grps=({0}, {1}), timings=None,
    )
    surface_timings = {}
    surfaces_with_empty_timing = get_build_surfs(
        b_verts, b_edges, V_BALLS, v_edges, e_balls, time.time(),
        interface=True, iface_grps=({0}, {1}), timings=surface_timings,
    )
    assert surfaces_without_timing == surfaces_with_empty_timing
    assert surface_timings


def test_repeated_builds_are_identical():
    first = _build(timing=False, interface=True)
    second = _build(timing=False, interface=True)
    assert first == second


def test_interface_diagnostic_is_opt_in_and_machine_readable(capsys):
    net = _net(False)
    net.settings['interface_diagnostics'] = True
    result = build(
        V_BALLS,
        V_LOCS,
        V_DUBS,
        num_balls=5,
        my_time=time.time(),
        interface=True,
        iface_grps=({0}, {1}),
        net=net,
    )

    report = net.interface_topology_diagnostics
    assert report['total_interface_facets'] == 1
    assert report['selected_physical_facets'] == 1
    assert report['connected_components'] == 1
    assert report['boundary_edges'] == 3
    assert report['stage_counts']['physical_faces_retained'] == {
        'vertices': 3, 'edges': 3, 'faces': 1,
    }
    assert result[3]['balls'] == [[0, 1]]

    line = next(
        line for line in capsys.readouterr().out.splitlines()
        if line.startswith('VORPY_INTERFACE_DIAGNOSTICS ')
    )
    assert json.loads(line.removeprefix('VORPY_INTERFACE_DIAGNOSTICS ')) == report
