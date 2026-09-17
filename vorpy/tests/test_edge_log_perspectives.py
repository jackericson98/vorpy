import numpy as np
import pandas as pd

from vorpy.src.output.logs import EDGE_CURVATURE_HEADERS, edge_curvature_log_values
from vorpy.src.inputs.net import _standardize_log_geometry_columns


def test_edge_log_perspectives_round_trip_without_aggregate_or_invented_zero():
    balls = [10, 20, 30]
    mean = {10: 1., 20: 2., 30: 3.}
    gaussian = {10: .3, 20: .5, 30: .6}
    faces = {(10, 20): .1, (10, 30): .2, (20, 10): .1,
             (20, 30): .4, (30, 10): .2, (30, 20): .4}
    values = edge_curvature_log_values(balls, mean, gaussian, faces, lambda x: x)
    assert values[:6] == [1., 2., 3., .3, .5, .6]
    assert len(values) == len(EDGE_CURVATURE_HEADERS)
    row = dict(zip(EDGE_CURVATURE_HEADERS, values))
    row.update({f"Ball {i + 1}": ball for i, ball in enumerate(balls)})
    _, _, loaded, _ = _standardize_log_geometry_columns(
        pd.DataFrame(), pd.DataFrame(), pd.DataFrame([row]), pd.DataFrame(),
    )
    assert loaded.iloc[0]["int_mean_curv_by_ball"] == mean
    assert loaded.iloc[0]["int_gauss_curv_by_generator"] == gaussian
    assert loaded.iloc[0]["int_gauss_curv_by_face"] == faces
    missing = edge_curvature_log_values(balls, mean, {10: np.nan}, None, lambda x: x)
    assert missing[3:] == [""] * 9
