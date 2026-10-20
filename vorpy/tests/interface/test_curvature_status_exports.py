from io import StringIO
from types import SimpleNamespace

import pandas as pd

from vorpy.src.interface.export import (
    summarize_direct_interface_orientation,
    write_surface_statistics,
    write_interface_orientation_statistics,
)
from vorpy.src.interface.water import summarize_water_group
from vorpy.src.output.logs import edge_curvature_log_values


def _interface(curvature_calculated):
    balls = pd.DataFrame([
        {"num": 0, "system_num": 0, "loc": (0.0, 0.0, 0.0), "rad": 1.0},
        {"num": 1, "system_num": 1, "loc": (1.0, 0.0, 0.0), "rad": 1.0},
    ])
    network = SimpleNamespace(
        curvature_calculated=curvature_calculated,
        balls=balls,
    )
    return SimpleNamespace(
        net=network,
        group1_indices=[0],
        group2_indices=[1],
        group1=SimpleNamespace(name="A"),
        group2=SimpleNamespace(name="B"),
    )


def _surface(**values):
    row = {
        "sa": 2.0,
        "balls": (0, 1),
        "flat": True,
        "com": (0.0, 0.0, 0.0),
        "func": (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0),
        "int_mean_curv": 0.0,
        "int_mean_curv_sq": 0.0,
        "int_gauss_curv": 0.0,
    }
    row.update(values)
    return pd.DataFrame([row])


def test_interface_summary_distinguishes_unavailable_from_zero():
    missing = summarize_direct_interface_orientation(
        _interface(False), _surface(), include_surface_energy=True
    )
    assert missing["curvature_status"] == "NOT_CALCULATED"
    assert missing["q"] is None
    assert missing["x"] is None
    assert missing["surf_energy"] is None

    zero = summarize_direct_interface_orientation(
        _interface(True), _surface(), include_surface_energy=True
    )
    assert zero["curvature_status"] == "CERTIFIED"
    assert zero["q"] == 0.0
    assert zero["x"] == 0.0
    assert zero["surf_energy"] == 0.0


def test_interface_log_marks_missing_curvature_as_unavailable():
    info = StringIO()
    write_interface_orientation_statistics(info, _interface(False), _surface())
    text = info.getvalue()
    assert "Integrated K dA (orientation invariant): Not available" in text
    assert "Integrated H² dA (orientation invariant): Not available" in text


def test_surface_statistics_do_not_sum_partial_curvature_rows():
    info = StringIO()
    write_surface_statistics(
        info,
        "Interface surfaces",
        _surface().drop(columns=["int_gauss_curv"]),
    )
    text = info.getvalue()
    assert "Integrated Gaussian curvature (∫K dA): Not available" in text


def test_water_summary_does_not_invent_curvature_zero():
    group = SimpleNamespace(
        name="water",
        net=SimpleNamespace(
            curvature_calculated=False,
            verts=[],
            edges=[],
            surfs=[],
        ),
        sa=3.0,
        vol=1.0,
        ball_ndxs=[],
        oriented_int_mean_curv=0.0,
        int_mean_curv_sq=0.0,
        int_gauss_curv=0.0,
        gaussian_curvature_status="NOT_CALCULATED",
        convex_sa=0.0,
        concave_sa=0.0,
        flat_sa=3.0,
        curved_unoriented_sa=0.0,
    )

    summary = summarize_water_group(group)
    assert summary["curvature_status"] == "NOT_CALCULATED"
    assert summary["int_gauss_curv"] is None
    assert summary["surf_energy"] is None


def test_log_curvature_aggregation_leaves_missing_terms_blank():
    values = edge_curvature_log_values(
        balls=(0, 1, 2), mean=None, gaussian=None, faces=None,
        rounding=lambda value: value,
    )
    assert values
    assert all(value == "" for value in values)
