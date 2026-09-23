"""Shared viewer/export scalar normalization; no rendering dependency."""
import numpy as np
from vorpy.workbench.domain import GeometryLayer


class DisplayColors:
    def __init__(self, layers):
        self._layer_definitions = {layer.name: layer for layer in layers}

    def _layer_scalar_limit(self, scheme: str, interpretation: str) -> float:
        values = []
        for layer in self._layer_definitions.values():
            if layer.interpretation != interpretation or scheme not in layer.cell_scalars:
                continue
            array = np.asarray(layer.cell_scalars[scheme], dtype=float).ravel()
            finite = array[np.isfinite(array)]
            if len(finite):
                values.extend(np.abs(finite))
        if not values:
            return 1.0
        limit = float(np.percentile(np.asarray(values, dtype=float), 98.0))
        return limit if np.isfinite(limit) and limit > 0.0 else 1.0

    def _display_scalars(self, layer: GeometryLayer, raw: np.ndarray) -> tuple[np.ndarray, tuple[float, float]]:
        values = np.asarray(raw, dtype=float).reshape(-1)
        finite = values[np.isfinite(values)]
        replacement = float(np.median(finite)) if len(finite) else 0.0
        values = np.nan_to_num(values, nan=replacement, posinf=replacement, neginf=replacement)
        if layer.color_scheme == "inside_outside":
            return np.clip(values, 0.0, 1.0), (0.0, 1.0)

        if layer.color_scheme == "distance":
            low, high = np.percentile(values, (2.0, 98.0)) if len(values) else (0.0, 1.0)
            if np.isclose(low, high):
                padding = max(abs(float(low)) * 0.01, 1e-9)
                low, high = low - padding, high + padding
            return values, (float(low), float(high))

        limit = self._layer_scalar_limit(layer.color_scheme, layer.interpretation)
        ratio = np.clip(values / limit, -1.0, 1.0)
        if layer.interpretation == "magnitude":
            ratio = np.abs(ratio)
        scale_gains = {
            "log10": 10.0,
            "log100": 100.0,
            "signed_log": 1000.0,  # compatibility with older projects
            "log1000": 1000.0,
            "log10000": 10000.0,
            "log100000": 100000.0,
        }
        gain = scale_gains.get(layer.scale_mode)
        if gain is not None:
            ratio = (
                np.sign(ratio)
                * np.log1p(gain * np.abs(ratio))
                / np.log1p(gain)
            )
        if layer.interpretation == "magnitude":
            ratio = 0.5 + 0.5 * np.abs(ratio)
            return ratio, (0.0, 1.0)
        return ratio, (-1.0, 1.0)

