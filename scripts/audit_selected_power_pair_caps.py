"""Read-only two-ball audit for the selected dry-2KAI Power contacts."""

from __future__ import annotations

import csv
import json
import math
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def distance(atoms: dict[int, dict[str, str]], a: int, b: int) -> float:
    return math.sqrt(sum((float(atoms[a][key]) - float(atoms[b][key])) ** 2 for key in ("x", "y", "z")))


def cap(radius: float, partner_radius: float, d: float) -> tuple[float, float, str]:
    tol = 1e-8
    if d <= tol:
        if partner_radius > radius + tol:
            return 2 * radius, 4 * math.pi * radius * radius, "contained_full"
        if radius > partner_radius + tol:
            return 0.0, 0.0, "contained_empty"
        return 2 * radius, 4 * math.pi * radius * radius, "coincident_full"
    if d > radius + partner_radius + tol:
        return 0.0, 0.0, "empty"
    if abs(d - (radius + partner_radius)) <= tol:
        return 0.0, 0.0, "tangent"
    if partner_radius >= radius + d - tol:
        return 2 * radius, 4 * math.pi * radius * radius, "contained_full"
    if radius >= partner_radius + d - tol:
        return 0.0, 0.0, "contained_empty"
    # The retained cap is on the partner-facing side of the radical plane.
    # Its height is radius - plane_distance_from_this_center.
    height = (partner_radius * partner_radius - (d - radius) ** 2) / (2 * d)
    height = max(0.0, min(2 * radius, height))
    return height, 2 * math.pi * radius * height, "partial"


def classify(radius: float, partner_radius: float, d: float) -> str:
    tol = 1e-8
    if d <= tol:
        return "coincident_equal" if abs(radius - partner_radius) <= tol else (
            "partner_contains" if partner_radius > radius else "self_contains"
        )
    if d > radius + partner_radius + tol:
        return "disjoint"
    if abs(d - (radius + partner_radius)) <= tol:
        return "externally_tangent"
    if abs(d - abs(radius - partner_radius)) <= tol:
        return "internally_tangent"
    if d < abs(radius - partner_radius) - tol:
        return "partner_contains" if partner_radius > radius else "self_contains"
    return "partial_overlap"


def stats(values: list[float]) -> dict[str, float | int | None]:
    values = sorted(values)
    if not values:
        return {"n": 0, "min": None, "q25": None, "median": None, "q75": None, "max": None, "mean": None}

    def quantile(position: float) -> float:
        index = (len(values) - 1) * position
        lower = math.floor(index)
        upper = math.ceil(index)
        return values[lower] + (values[upper] - values[lower]) * (index - lower)

    return {
        "n": len(values),
        "min": values[0],
        "q25": quantile(0.25),
        "median": quantile(0.5),
        "q75": quantile(0.75),
        "max": values[-1],
        "mean": statistics.fmean(values),
    }


def main() -> None:
    with (ROOT / "2kai_cazals_atoms.csv").open(newline="", encoding="utf-8") as handle:
        atoms = {int(row["index"]): row for row in csv.DictReader(handle)}
    with (ROOT / "2kai_alpha0_edges.csv").open(newline="", encoding="utf-8") as handle:
        selected_rows = list(csv.DictReader(handle))
    selected = {
        tuple(sorted((int(row["facet_1_i"]), int(row["facet_1_j"]))))
        for row in selected_rows
    }
    selected |= {
        tuple(sorted((int(row["facet_2_i"]), int(row["facet_2_j"]))))
        for row in selected_rows
    }

    records = []
    for left, right in sorted(selected):
        if atoms[left]["group"] == atoms[right]["group"]:
            continue
        a, b = (left, right) if atoms[left]["group"] == "A" else (right, left)
        radius_a = float(atoms[a]["expanded_radius_A"])
        radius_b = float(atoms[b]["expanded_radius_A"])
        d = distance(atoms, a, b)
        height_a, area_a, cap_state_a = cap(radius_a, radius_b, d)
        height_b, area_b, cap_state_b = cap(radius_b, radius_a, d)
        records.append({
            "a": a,
            "b": b,
            "d": d,
            "radius_a": radius_a,
            "radius_b": radius_b,
            "class": classify(radius_a, radius_b, d),
            "depth": radius_a + radius_b - d,
            "height_a": height_a,
            "height_b": height_b,
            "area_a": area_a,
            "area_b": area_b,
            "cap_state_a": cap_state_a,
            "cap_state_b": cap_state_b,
        })

    class_counts = {}
    for record in records:
        class_counts[record["class"]] = class_counts.get(record["class"], 0) + 1
    all_cross_counts = {}
    group_a = [index for index, row in atoms.items() if row["group"] == "A"]
    group_b = [index for index, row in atoms.items() if row["group"] == "B"]
    for a in group_a:
        for b in group_b:
            state = classify(
                float(atoms[a]["expanded_radius_A"]),
                float(atoms[b]["expanded_radius_A"]),
                distance(atoms, a, b),
            )
            all_cross_counts[state] = all_cross_counts.get(state, 0) + 1
    summary = {
        "selected_contacts": len(records),
        "class_counts": class_counts,
        "all_cross_group_pair_counts": all_cross_counts,
        "nonempty_a_caps": sum(record["area_a"] > 1e-12 for record in records),
        "nonempty_b_caps": sum(record["area_b"] > 1e-12 for record in records),
        "overlap_depth_Ra_plus_Rb_minus_d_A": stats([record["depth"] for record in records]),
        "distance_A": stats([record["d"] for record in records]),
        "cap_height_A": stats([record["height_a"] for record in records]),
        "cap_height_B": stats([record["height_b"] for record in records]),
        "pairwise_cap_area_A_A2": stats([record["area_a"] for record in records]),
        "pairwise_cap_area_B_A2": stats([record["area_b"] for record in records]),
        "pairwise_cap_area_total_A2": stats([record["area_a"] + record["area_b"] for record in records]),
        "depth_bins_A": {
            "<=0": sum(record["depth"] <= 0 for record in records),
            "0-1": sum(0 < record["depth"] < 1 for record in records),
            "1-2": sum(1 <= record["depth"] < 2 for record in records),
            "2-3": sum(2 <= record["depth"] < 3 for record in records),
            ">=3": sum(record["depth"] >= 3 for record in records),
        },
    }
    output = ROOT / "output" / "dual_side_incidence_audit" / "2kai_power" / "selected_power_pair_caps.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
