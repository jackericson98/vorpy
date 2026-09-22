import json

import pandas as pd

from vorpy.src.analyze.nonpolar_interface import (
    export_geometry_tables,
    summarize_surface_features,
)


def tables():
    return {
        "surfaces": pd.DataFrame([
            {"nonpolar_surface_label": "nonpolar", "interface_label": 1,
             "surface_area": 2., "mean_curvature": .5},
            {"nonpolar_surface_label": "nonpolar", "interface_label": 0,
             "surface_area": 4., "mean_curvature": .1},
            {"nonpolar_surface_label": "mixed", "interface_label": 1,
             "surface_area": 1., "mean_curvature": .2},
        ]),
        "atoms": pd.DataFrame([{"atom_index": 1, "interface_area": 2.}]),
        "metadata": {"system_id": "fixture", "frame_id": 0},
    }


def test_export_writes_two_tables_and_metadata(tmp_path):
    paths = export_geometry_tables(tables(), tmp_path)
    assert set(paths) == {"surfaces", "atoms", "metadata"}
    assert paths["surfaces"].is_file()
    assert pd.read_csv(paths["surfaces"]).shape == (3, 4)
    assert json.loads(paths["metadata"].read_text())["system_id"] == "fixture"


def test_summary_reports_distribution_statistics():
    summary = summarize_surface_features(tables()["surfaces"])
    row = summary[(summary["nonpolar_surface_label"] == "nonpolar")
                  & (summary["interface_label"] == 1)
                  & (summary["feature"] == "surface_area")].iloc[0]
    assert row["count"] == 1
    assert row["median"] == 2.
    assert row["iqr"] == 0.
