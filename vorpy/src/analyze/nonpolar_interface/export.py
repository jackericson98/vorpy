"""Stable file export for nonpolar/interface geometry analysis tables."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def _json_value(value):
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def export_geometry_tables(tables, destination, stem="nonpolar_interface_geometry", parquet=False):
    """Write surface/atom CSV files and JSON metadata, returning their paths.

    Parquet output is opt-in because pandas' parquet engines are not a VorPy
    dependency. If requested but unavailable, the function raises a clear
    error instead of silently changing the requested format.
    """
    required = {"surfaces", "atoms", "metadata"}
    missing = required.difference(tables)
    if missing:
        raise ValueError(f"Analysis tables are missing: {sorted(missing)}")
    output = Path(destination)
    output.mkdir(parents=True, exist_ok=True)
    paths = {}
    for name in ("surfaces", "atoms"):
        frame = tables[name]
        suffix = ".parquet" if parquet else ".csv"
        path = output / f"{stem}_{name}{suffix}"
        if parquet:
            try:
                frame.to_parquet(path, index=False)
            except ImportError as error:
                raise RuntimeError(
                    "Parquet export requires an installed pandas parquet engine "
                    "such as pyarrow or fastparquet."
                ) from error
        else:
            frame.to_csv(path, index=False)
        paths[name] = path
    metadata_path = output / f"{stem}_metadata.json"
    metadata_path.write_text(
        json.dumps(_json_value(tables["metadata"]), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    paths["metadata"] = metadata_path
    return paths
