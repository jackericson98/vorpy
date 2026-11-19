from importlib import import_module
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest


@pytest.mark.parametrize('interface, timeout', [(None, None), (({0}, {1}), 0.05)])
def test_disconnected_group_reseed_is_not_discarded_by_a_wall_clock_cutoff(monkeypatch, interface, timeout):
    module = import_module('vorpy.src.network.find_net_verts')
    net = SimpleNamespace(group=[0, 1], iface_grps=interface,
        balls=pd.DataFrame({'num': range(5), 'loc': [np.array([i, 0., 0.]) for i in range(5)],
                            'rad': [1.] * 5}),
        settings={'net_type': 'aw', 'max_vert': 5., 'foam_box': None}, metrics={'start': 0.},
        box={'verts': []}, update_progress=lambda *args: None)
    calls = []

    def find(**kwargs):
        calls.append(kwargs)
        if kwargs.get('b0') is None:
            return ([[0, 2, 3, 4]], [[0., 0., 0.]], [1.], [[None] * 3], [None], [1], [[]] * 5)
        if kwargs['seed_timeout'] is not None and interface is None:
            return None  # Reproduce a valid second component lost under the old 50-ms budget.
        return ([[0, 2, 3, 4], [1, 2, 3, 4]], [[0., 0., 0.], [1., 0., 0.]],
                [1., 1.], [[None] * 3] * 2, [None, None], [], [[]] * 5)

    monkeypatch.setattr(module, 'find_verts', find)
    monkeypatch.setattr(module, '_load_cached_vertex_state', lambda _: None)
    monkeypatch.setattr(module, 'get_balls', lambda *args, **kwargs: None)
    monkeypatch.setattr(module, 'box_search', lambda *args: None)
    monkeypatch.setattr(module, '_get_interface_reseed_candidates', lambda net, candidates: candidates)
    monkeypatch.setattr(module, '_store_native_vertex_state', lambda *args: None)
    monkeypatch.setattr(module, 'write_verts', lambda _: None)
    module.find_net_verts(net)
    assert len(net.verts) == 2
    assert calls[1]['b0'] == 1 and calls[1]['seed_timeout'] == timeout
