"""Descriptive, system-aware summaries for the first geometry comparison."""

from __future__ import annotations

import pandas as pd


DEFAULT_FEATURES = (
    "surface_area",
    "mean_curvature",
    "average_mean_curvature",
    "gaussian_curvature",
    "average_gaussian_curvature",
    "integrated_mean_curvature",
    "integrated_mean_curvature_squared",
    "integrated_gaussian_curvature",
)


def summarize_surface_features(surfaces, features=DEFAULT_FEATURES):
    """Return count/mean/median/std/IQR summaries by chemistry and label.

    The grouping keys are retained in the output. This is descriptive only;
    it does not treat rows from one molecular system as independent evidence.
    """
    if not isinstance(surfaces, pd.DataFrame):
        raise TypeError("surfaces must be a pandas DataFrame")
    if surfaces.empty:
        return pd.DataFrame(columns=[
            "nonpolar_surface_label", "interface_label", "feature",
            "count", "mean", "median", "std", "q1", "q3", "iqr",
        ])
    rows = []
    grouping = ["nonpolar_surface_label", "interface_label"]
    for keys, group in surfaces.groupby(grouping, dropna=False, sort=True):
        for feature in features:
            if feature not in group:
                continue
            values = pd.to_numeric(group[feature], errors="coerce").dropna()
            if values.empty:
                continue
            q1, q3 = values.quantile([0.25, 0.75])
            rows.append({
                grouping[0]: keys[0], grouping[1]: keys[1],
                "feature": feature, "count": int(values.size),
                "mean": float(values.mean()), "median": float(values.median()),
                "std": float(values.std(ddof=1)) if values.size > 1 else 0.0,
                "q1": float(q1), "q3": float(q3), "iqr": float(q3 - q1),
            })
    return pd.DataFrame(rows)
