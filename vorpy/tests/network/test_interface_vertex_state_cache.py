import numpy as np
import pandas as pd
import pytest

from vorpy.src.io.network_archive import load_network, save_network
from vorpy.src.network.find_net_verts import (
    _load_cached_vertex_state,
    _store_native_vertex_state,
)
from vorpy.src.network.network import Network
from vorpy.src.system.system import System


def _vertex_state():
    return {
        'balls': [0, 1, 2, 3],
        'loc': [0.25, 0.25, 0.25],
        'rad': 0.2,
        'loc2': [None, None, None],
        'rad2': None,
    }


def _network():
    locations = [
        (0.0, 0.0, 0.0),
        (2.0, 0.0, 0.0),
        (0.0, 2.0, 0.0),
        (0.0, 0.0, 2.0),
    ]
    balls = pd.DataFrame({
        'loc': locations,
        'rad': [1.0] * 4,
        'num': list(range(4)),
        'complete': [True] * 4,
        'stable_id': [f'atom-{index}' for index in range(4)],
    })
    verts = pd.DataFrame({
        'loc': [(0.25, 0.25, 0.25)],
        'rad': [0.2],
        'balls': [[0, 1, 2, 3]],
    })
    edges = pd.DataFrame({
        'balls': [[0, 1, 2]],
        'verts': [[0]],
        'points': [[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]],
    })
    surfs = pd.DataFrame({
        'balls': [[0, 1, 2, 3]],
        'verts': [[0]],
        'edges': [[0]],
        'points': [[
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        ]],
        'tris': [[(0, 1, 2), (0, 2, 3)]],
    })
    system = System(simple=True, name='cache-round-trip', make_dir=False)
    system.balls = balls.copy(deep=True)
    system.balls['res'] = None
    system.balls['chn'] = None
    return Network(
        locs=locations,
        rads=[1.0] * 4,
        group=[0, 1, 2, 3],
        iface_grps=({0, 1}, {2, 3}),
        group_name='synthetic-interface',
        settings={'net_type': 'aw', 'max_vert': 5.0},
        balls=balls,
        verts=verts,
        edges=edges,
        surfs=surfs,
        system=system,
    )


def _store_state(network):
    _store_native_vertex_state(
        network,
        [[0, 1, 2, 3]],
        [[0.25, 0.25, 0.25]],
        [0.2],
        [[None, None, None]],
        [None],
    )


def test_loads_archived_native_interface_vertex_state_and_reports_hit():
    network = _network()
    _store_state(network)

    loaded = _load_cached_vertex_state(network)

    assert loaded is not None
    vert_ndxs, vlocs, vrads, vloc2s, vrad2s, b_verts = loaded
    assert vert_ndxs == [[0, 1, 2, 3]]
    assert np.allclose(vlocs, [[0.25, 0.25, 0.25]])
    assert vrads == [0.2]
    assert vloc2s == [[None, None, None]]
    assert vrad2s == [None]
    assert b_verts == [[0], [0], [0], [0]]
    assert network.interface_vertex_cache_stats == {
        'hits': 1,
        'misses': 0,
        'source': 'interface_vertex_state',
        'records': 1,
        'accepted': 1,
        'reason': None,
    }


def test_loads_legacy_cached_wrapper_shape():
    network = _network()
    _store_state(network)
    network.cached_interface_vertex_state = [
        {'source_index': 4, 'data': _vertex_state()},
    ]
    network.interface_vertex_state = None

    loaded = _load_cached_vertex_state(network)

    assert loaded is not None
    assert loaded[0] == [[0, 1, 2, 3]]
    assert network.interface_vertex_cache_stats['source'] == 'cached_interface_vertex_state'
    assert network.interface_vertex_cache_stats['hits'] == 1


def test_malformed_cached_state_fails_closed_without_native_fallback():
    network = _network()
    _store_state(network)
    network.cached_interface_vertex_state = [
        {'source_index': 0, 'data': {**_vertex_state(), 'loc2': 1.0}},
    ]
    network.interface_vertex_state = [_vertex_state()]

    assert _load_cached_vertex_state(network) is None
    assert network.interface_vertex_cache_stats['hits'] == 0
    assert network.interface_vertex_cache_stats['misses'] == 1
    assert network.interface_vertex_cache_stats['reason'] == 'invalid_secondary_geometry'


def test_interface_cache_is_ignored_outside_interface_scope():
    network = _network()
    network.iface_grps = None
    network.interface_vertex_state = [_vertex_state()]

    assert _load_cached_vertex_state(network) is None
    assert network.interface_vertex_cache_stats == {
        'hits': 0,
        'misses': 1,
        'source': None,
        'records': 0,
        'accepted': 0,
        'reason': 'non_interface_scope',
    }


def test_native_state_without_dependency_fails_closed():
    network = _network()
    network.interface_vertex_state = [_vertex_state()]

    assert _load_cached_vertex_state(network) is None
    assert network.interface_vertex_cache_stats['reason'] == 'dependency_missing'


@pytest.mark.parametrize('change', ('coordinates', 'radii', 'net_type'))
def test_native_state_invalidates_when_solve_dependency_changes(change):
    network = _network()
    _store_state(network)

    if change == 'coordinates':
        network.balls.at[0, 'loc'] = (0.5, 0.0, 0.0)
    elif change == 'radii':
        network.balls.at[0, 'rad'] = 1.25
    else:
        network.settings['net_type'] = 'pow'

    assert _load_cached_vertex_state(network) is None
    assert network.interface_vertex_cache_stats['reason'] == 'dependency_mismatch'


def test_native_interface_vertex_state_round_trip_keeps_telemetry_runtime_only(tmp_path):
    network = _network()
    _store_state(network)
    archive = save_network(network, tmp_path / 'interface-vertex-state.vpy')

    restored = load_network(archive)

    assert restored.interface_vertex_state == network.interface_vertex_state
    assert restored.interface_vertex_state_dependency == network.interface_vertex_state_dependency
    assert not hasattr(restored, 'cached_interface_vertex_state')
    assert not hasattr(restored, 'interface_vertex_cache_stats')
    assert _load_cached_vertex_state(restored) is not None
    assert restored.interface_vertex_cache_stats['hits'] == 1


def test_repeated_archive_loads_report_independent_runtime_hits(tmp_path):
    network = _network()
    _store_state(network)
    archive = save_network(network, tmp_path / 'repeated-interface-vertex-state.vpy')

    first = load_network(archive)
    second = load_network(archive)

    for restored in (first, second):
        assert not hasattr(restored, 'interface_vertex_cache_stats')
        assert _load_cached_vertex_state(restored) is not None
        assert restored.interface_vertex_cache_stats['hits'] == 1
        assert restored.interface_vertex_cache_stats['misses'] == 0
