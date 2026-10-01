"""Regression on the actual VorPy union, including its blind axial floor."""
import json
from importlib.util import find_spec
from types import SimpleNamespace

import numpy as np
import pytest

from vorpy.src.geometry.topology_diagnostic import run_deep_pocket_diagnostic


# The real build exceeds the suite's 30-second timeout on some machines.
if find_spec("pytest_timeout") is not None:
    pytestmark = pytest.mark.timeout(120)


def test_measured_deep_pocket_is_one_closed_genus_zero_boundary(tmp_path):
    args = SimpleNamespace(output=tmp_path, body_radius=10., pocket_mouth_radius=3.2,
                           pocket_depth=14., pocket_neck_radius=2.1,
                           outer_count=48, pocket_layers=4, ring_count=8)
    run_deep_pocket_diagnostic(args)
    report = json.loads((tmp_path / "report.json").read_text())
    assert report["acceptance_errors"] == []
    assert report["selected_cell_components"] == report["boundary_components"] == 1
    assert len(report["component_reports"]) == 1
    assert report["inferred_genus"] == 0
    assert report["euler_characteristic"] == 2
    assert np.isclose(report["integrated_gaussian_curvature"], 4*np.pi, rtol=0, atol=1e-8)
    assert report["complete_selected_cells"] == report["selected_cells"]
    assert all(report[k] for k in ("complete", "closed", "manifold", "orientable"))
    underside, floor = report["axis_boundary_intersections_z"]
    assert underside < floor < 0
    assert report["measured_axial_floor_thickness"] > 2
    assert (tmp_path / "presentation_mesh.json").exists()
    assert 'cgo_transparency", 0' in (tmp_path / "presentation.py").read_text()
