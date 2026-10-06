"""Adapters from frozen scheme-specific curvature results to common fields."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CurvatureSummary:
    scheme: str
    selection_mode: str
    area_A2: float | None
    smooth_H_A: float | None
    raw_signed_edge_A_rad: float | None
    raw_unsigned_edge_A_rad: float | None
    positive_edge_A_rad: float | None
    negative_edge_A_rad: float | None
    conventional_edge_H_A: float | None
    combined_H_A: float | None
    H_per_area_inv_A: float | None
    unsigned_edge_per_area_rad_per_A: float | None
    cancellation_fraction: float | None
    orientation: str | None
    status: str = "available"
    source: str = ""


def _summary(scheme, selection_mode, area, smooth, signed, unsigned,
             positive, negative, orientation, source):
    edge_half = None if signed is None else 0.5 * float(signed)
    combined = None if smooth is None or edge_half is None else float(smooth) + edge_half
    area_value = None if area is None else float(area)
    h_area = (combined / area_value if combined is not None and area_value is not None
              and area_value > 0 else None)
    unsigned_area = (float(unsigned) / area_value if unsigned is not None
                     and area_value is not None and area_value > 0 else None)
    cancellation = (1.0 - abs(float(signed)) / float(unsigned)
                    if signed is not None and unsigned is not None and float(unsigned) > 0 else None)
    return CurvatureSummary(
        scheme, selection_mode, area_value,
        None if smooth is None else float(smooth),
        None if signed is None else float(signed),
        None if unsigned is None else float(unsigned),
        None if positive is None else float(positive),
        None if negative is None else float(negative),
        edge_half, combined, h_area, unsigned_area, cancellation,
        orientation, "available", source,
    )


def aw_curvature_summary(result, *, area_A2=None):
    """Adapt an existing orientation-aware AW result; perform no integrations."""
    area = result.interface_area if area_A2 is None else area_A2
    edge = result.edge_curvature
    selection = "full" if result.selection_mode == "raw" else result.selection_mode
    return _summary(
        "aw", selection, area, result.surface_curvature,
        edge["signed"], edge["unsigned"], edge["positive"], edge["negative"],
        "interface normal A-to-B; reversal flips signed smooth/edge/combined terms",
        "vorpy.src.analyze.aw_interface_curvature.AWInterfaceCurvatureResult",
    )


def power_curvature_summary(result, *, area_A2=None):
    """Adapt existing Power curvature totals without changing its selection."""
    area = area_A2
    totals = result.edge_totals
    return _summary(
        "power", "cazals_alpha0_power", area, result.smooth_surface_curvature,
        totals["signed"], totals["unsigned"], totals["positive"], totals["negative"],
        "Cazals sign: ABB positive, AAB negative; input groups A/B fixed",
        "vorpy.src.analyze.power_interface_curvature.PowerInterfaceCurvatureResult",
    )
