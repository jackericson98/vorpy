import importlib
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from vorpy.src.output.curvature_colors import curvature_color_limit, export_color_cache
from vorpy.src.output.draw import draw_joint
from vorpy.src.output.output import _group_export_cache


def test_export_color_scales_are_refreshed_between_plans():
    net = SimpleNamespace(
        surfs=pd.DataFrame([{"balls": [0, 1], "int_mean_curv_by_ball": {0: 2., 1: -2.}}]),
        edges=pd.DataFrame(), verts=pd.DataFrame(),
    )
    with export_color_cache(net):
        assert curvature_color_limit(net, "int_mean_curv", [0]) == 2.
        net.surfs.at[0, "int_mean_curv_by_ball"] = {0: 7., 1: -7.}
        assert curvature_color_limit(net, "int_mean_curv", [0, 0]) == 2.
        assert curvature_color_limit(net, "int_mean_curv", [1]) == 7.
    assert not hasattr(net, "_export_color_limits")
    assert curvature_color_limit(net, "int_mean_curv", [0]) == 7.


def test_export_summary_is_reused_and_cache_cleared_on_failure(monkeypatch):
    module = importlib.import_module("vorpy.src.group.group")
    group = SimpleNamespace(net=SimpleNamespace())
    calls = []
    monkeypatch.setattr(module, "get_info", lambda grp: calls.append(grp))
    with pytest.raises(RuntimeError):
        with _group_export_cache(group):
            module.Group.get_info(group)
            module.Group.get_info(group)
            assert len(calls) == 1
            raise RuntimeError("Export failed")
    assert not hasattr(group, "_export_info_ready")
    assert not hasattr(group.net, "_export_color_limits")
    module.Group.get_info(group)
    assert len(calls) == 2


@pytest.mark.parametrize("subdivisions", [0, 1, 2])
def test_joint_template_preserves_radius_translation_and_independent_results(subdivisions):
    points, triangles = draw_joint([0., 0., 0.], radius=1., subdivisions=subdivisions)
    moved, moved_triangles = draw_joint([2., 3., 4.], radius=0.5, subdivisions=subdivisions)
    assert np.asarray(moved) == pytest.approx(np.asarray(points) * 0.5 + [2., 3., 4.])
    assert np.linalg.norm(points, axis=1) == pytest.approx(np.ones(len(points)))
    assert len(triangles) == 20 * 4 ** subdivisions
    assert triangles == moved_triangles
    triangles[0][0] = -1
    assert moved_triangles[0][0] >= 0
