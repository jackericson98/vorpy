import pandas as pd

from vorpy.src.command.vpy_cmnd import execution_profile_cli_options
from vorpy.src.network.network import Network


def test_execution_profile_option_is_opt_in_and_validated():
    args, profile = execution_profile_cli_options(
        ['input.pdb', '--profile', 'interface', '-i']
    )
    assert args == ['input.pdb', '-i']
    assert profile == 'interface'


def test_execution_profile_option_rejects_unknown_profile():
    try:
        execution_profile_cli_options(['--profile', 'fast'])
    except SystemExit as error:
        assert 'geometry, interface, analysis, or full' in str(error)
    else:
        raise AssertionError('Unknown execution profile was accepted')


def test_geometry_only_network_build_skips_scientific_stages(monkeypatch):
    network = Network.__new__(Network)
    network.group_name = 'synthetic'
    network.completion_kind = 'network'
    network.settings = {
        'build_type': 'all', 'net_type': 'aw', 'verbose': False,
    }
    network.metrics = {'start': 0.0}
    network.box = {}
    network.group = [0]
    network.balls = pd.DataFrame({'complete': [False]})
    network.verts = pd.DataFrame({'vdub': [0]})
    network.edges = pd.DataFrame(index=[])
    network.surfs = pd.DataFrame(index=[])
    network.update_progress = lambda *args, **kwargs: None
    events = []

    network.connect = lambda: events.append('connect')
    network.build_edges = lambda **kwargs: events.append(('edges', kwargs))
    network.build_surfaces = lambda *args, **kwargs: events.append(('surfaces', kwargs))
    network.analyze = lambda: events.append('analyze')
    monkeypatch.setattr(
        'vorpy.src.network.network.analyze_geometry_only',
        lambda value: events.append(('geometry_only', value)),
    )

    Network.build(network, calculate_curvature=False)

    assert ('edges', {'compute_curvature': False}) in events
    assert ('surfaces', {'calculate_curvature': False}) in events
    assert ('geometry_only', network) in events
    assert 'analyze' not in events
