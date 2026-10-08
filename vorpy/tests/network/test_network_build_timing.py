import json

import pandas as pd

import vorpy.src.network.network as network_module
from vorpy.src.network.network import Network


def _network(*, timing):
    network = Network.__new__(Network)
    network.group_name = 'timing-test'
    network.completion_kind = 'network'
    network.settings = {
        'build_type': 'all',
        'net_type': 'pow',
        'verbose': False,
        'timing': timing,
    }
    network.metrics = {'start': 0.0}
    network.box = {}
    network.sys = None
    network.group = [0]
    network.balls = pd.DataFrame({'num': [0], 'complete': [False]})
    network.verts = pd.DataFrame({'vdub': [0]})
    network.edges = pd.DataFrame(index=[])
    network.surfs = pd.DataFrame(index=[])
    network.update_progress = lambda *args, **kwargs: None
    network.connect = lambda: None
    network.build_edges = lambda **kwargs: None
    network.build_surfaces = lambda *args, **kwargs: None
    return network


def test_network_build_timing_is_opt_in_and_json_serializable(monkeypatch, capsys):
    monkeypatch.delenv('VORPY_TIMING', raising=False)
    monkeypatch.setattr(network_module, 'analyze_geometry_only', lambda net: None)
    network = _network(timing=True)

    Network.build(network, calculate_curvature=False)

    payload = network.build_timing
    assert payload['schema_version'] == 1
    assert payload['unit'] == 'seconds'
    assert set(payload['stages']) == {
        'topology', 'edge_generation', 'surface_construction', 'completeness_analysis',
    }
    assert set(payload['breakdown_seconds']) == {
        'geometric_predicates', 'surface_mesh_generation', 'surface_curvature',
        'intermediate_conversions', 'completeness_analysis',
    }
    assert json.dumps(payload, allow_nan=False)

    line = next(
        line for line in capsys.readouterr().out.splitlines()
        if line.startswith('VORPY_BUILD_TIMING ')
    )
    assert json.loads(line.removeprefix('VORPY_BUILD_TIMING ')) == payload


def test_network_build_timing_does_not_change_default_output(monkeypatch, capsys):
    monkeypatch.delenv('VORPY_TIMING', raising=False)
    monkeypatch.setattr(network_module, 'analyze_geometry_only', lambda net: None)
    network = _network(timing=False)

    Network.build(network, calculate_curvature=False)

    assert not hasattr(network, 'build_timing')
    assert 'VORPY_BUILD_TIMING ' not in capsys.readouterr().out
