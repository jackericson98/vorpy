"""Independent topology audit for the frozen 2KAI contact surfaces.

This script reads the persisted curved geometry and does not change production
surface construction or area values.
"""

from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SURFACE_DIR = ROOT / "output" / "dual_side_incidence_audit" / "2kai_power" / "molecular_contact_surface"
TWO_PI = 2.0 * math.pi


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def scale(a, value):
    return (a[0] * value, a[1] * value, a[2] * value)


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def norm(a):
    return math.sqrt(dot(a, a))


def distance(a, b):
    return norm(sub(a, b))


def unit(a):
    return scale(a, 1.0 / norm(a))


def stable_id(row):
    return (
        f"{row['index']}|pdb={row['pdb_serial']}|group={row['group']}|"
        f"{row['chain']}:{row['residue_name']}{row['residue_number']}:{row['atom_name']}"
    )


def load_atoms():
    with (ROOT / "2kai_cazals_atoms.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    atoms = {}
    by_index = {}
    for row in rows:
        atom = {
            "id": stable_id(row),
            "index": int(row["index"]),
            "group": row["group"],
            "center": (float(row["x"]), float(row["y"]), float(row["z"])),
            "radius": float(row["expanded_radius_A"]),
        }
        atoms[atom["id"]] = atom
        by_index[atom["index"]] = atom
    return atoms, by_index


def load_selected_partners(by_index):
    partners = defaultdict(list)
    with (ROOT / "2kai_alpha0_edges.csv").open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            for left, right in ((int(row["facet_1_i"]), int(row["facet_1_j"])),
                                (int(row["facet_2_i"]), int(row["facet_2_j"]))):
                first, second = by_index[left], by_index[right]
                if first["group"] == second["group"]:
                    continue
                partners[first["id"]].append(second)
                partners[second["id"]].append(first)
    return {key: tuple({item["id"]: item for item in values}.values()) for key, values in partners.items()}


def merge_intervals(intervals):
    values = sorted((left, right) for left, right in intervals if right - left > 1e-10)
    merged = []
    for left, right in values:
        if not merged or left > merged[-1][1] + 1e-8:
            merged.append([left, right])
        else:
            merged[-1][1] = max(merged[-1][1], right)
    return tuple(tuple(value) for value in merged)


def intersection_intervals(first, second):
    return merge_intervals(
        (max(left_a, left_b), min(right_a, right_b))
        for left_a, right_a in first
        for left_b, right_b in second
        if min(right_a, right_b) - max(left_a, left_b) > 1e-8
    )


def split_interval(start, delta):
    if abs(delta) >= TWO_PI - 1e-8:
        return ((0.0, TWO_PI),)
    if delta >= 0.0:
        left, right = start % TWO_PI, (start + delta) % TWO_PI
    else:
        left, right = (start + delta) % TWO_PI, start % TWO_PI
    if right >= left:
        return ((left, right),)
    return ((left, TWO_PI), (0.0, right))


def circle_geometry(first, second):
    displacement = sub(second["center"], first["center"])
    separation = norm(displacement)
    normal = scale(displacement, 1.0 / separation)
    offset = (
        first["radius"] ** 2 - second["radius"] ** 2 + separation ** 2
    ) / (2.0 * separation)
    center = add(first["center"], scale(normal, offset))
    radius_squared = max(0.0, first["radius"] ** 2 - offset ** 2)
    reference = min(
        ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
        key=lambda value: abs(dot(value, normal)),
    )
    basis_1 = unit(cross(reference, normal))
    basis_2 = cross(normal, basis_1)
    return center, math.sqrt(radius_squared), basis_1, basis_2


def arc_interval(arc, atoms):
    first_id, second_id = sorted((arc["carrier_atom_id"], arc["source_atom_ids"][0]))
    first, second = atoms[first_id], atoms[second_id]
    center, circle_radius, basis_1, basis_2 = circle_geometry(first, second)
    if circle_radius <= 1e-12:
        return ()

    def theta(point):
        direction = scale(sub(point, center), 1.0 / circle_radius)
        return math.atan2(dot(direction, basis_2), dot(direction, basis_1))

    start = theta(arc["start_xyz"])
    end = theta(arc["end_xyz"])
    positive = (end - start) % TWO_PI
    negative = positive - TWO_PI
    target = abs(arc["angular_span"])
    delta = positive if abs(abs(positive) - target) <= abs(abs(negative) - target) else negative
    return split_interval(start, delta)


def endpoint_gap(left, right):
    direct = max(distance(left["start_xyz"], right["start_xyz"]), distance(left["end_xyz"], right["end_xyz"]))
    reverse = max(distance(left["start_xyz"], right["end_xyz"]), distance(left["end_xyz"], right["start_xyz"]))
    return min(direct, reverse)


def arc_midpoint(arc, atoms):
    carrier = atoms[arc["carrier_atom_id"]]
    normal = tuple(arc["circle_normal"])
    reference = min(
        ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
        key=lambda value: abs(dot(value, normal)),
    )
    basis_1 = unit(cross(reference, normal))
    basis_2 = cross(normal, basis_1)
    direction = unit(sub(arc["start_xyz"], carrier["center"]))
    theta = math.atan2(dot(direction, basis_2), dot(direction, basis_1))
    radial = math.sqrt(max(0.0, 1.0 - arc["circle_offset"] ** 2))
    unit_direction = add(
        scale(normal, arc["circle_offset"]),
        scale(
            add(
                scale(basis_1, math.cos(theta + arc["angular_span"] / 2.0)),
                scale(basis_2, math.sin(theta + arc["angular_span"] / 2.0)),
            ),
            radial,
        ),
    )
    return add(carrier["center"], scale(unit_direction, carrier["radius"]))


def validate_arc(arc, atoms, selected_partners):
    point = arc_midpoint(arc, atoms)
    carrier = atoms[arc["carrier_atom_id"]]
    carrier_residual = abs(distance(point, carrier["center"]) - carrier["radius"])
    source_residuals = {
        atom_id: abs(distance(point, atoms[atom_id]["center"]) - atoms[atom_id]["radius"])
        for atom_id in arc["source_atom_ids"]
    }
    self_occluded = any(
        other["group"] == carrier["group"]
        and other["id"] != carrier["id"]
        and distance(point, other["center"]) < other["radius"] - 1e-7
        for other in atoms.values()
    )
    partner_member = any(
        distance(point, partner["center"]) <= partner["radius"] + 1e-7
        for partner in selected_partners.get(carrier["id"], ())
    )
    valid = (
        carrier_residual <= 1e-8
        and max(source_residuals.values(), default=0.0) <= 1e-8
        and not self_occluded
        and (not arc["source_contact_ids"] or partner_member)
    )
    return {
        "valid": valid,
        "midpoint_xyz": point,
        "carrier_sphere_residual_A": carrier_residual,
        "source_sphere_residuals_A": source_residuals,
        "same_group_strict_occlusion": self_occluded,
        "selected_partner_membership": partner_member,
    }


class UnionFind:
    def __init__(self, size):
        self.parent = list(range(size))

    def find(self, value):
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, first, second):
        first, second = self.find(first), self.find(second)
        if first != second:
            self.parent[second] = first


def current_graph(surface):
    patches = surface["patches"]
    arcs = surface["arcs"]
    patch_index = {patch["patch_id"]: index for index, patch in enumerate(patches)}
    arc_patch = {
        arc_id: patch_index[patch["patch_id"]]
        for patch in patches
        for arc_id in patch["boundary_arc_ids"]
    }
    self_arcs = [
        arc for arc in arcs
        if not arc["source_contact_ids"] and len(arc["source_atom_ids"]) == 1
    ]
    matched = set()
    endpoint_pairs = 0
    for left_index, left in enumerate(self_arcs):
        for right in self_arcs[left_index + 1:]:
            if left["carrier_atom_id"] == right["carrier_atom_id"]:
                continue
            if left["source_atom_ids"][0] != right["carrier_atom_id"]:
                continue
            if right["source_atom_ids"][0] != left["carrier_atom_id"]:
                continue
            if endpoint_gap(left, right) <= 1e-7:
                endpoint_pairs += 1
                matched.update((left["arc_id"], right["arc_id"]))

    union = UnionFind(len(patches))
    for left_index, left in enumerate(self_arcs):
        if left["arc_id"] not in matched:
            continue
        for right in self_arcs[left_index + 1:]:
            if right["arc_id"] not in matched:
                continue
            if left["carrier_atom_id"] == right["carrier_atom_id"]:
                continue
            if left["source_atom_ids"][0] != right["carrier_atom_id"]:
                continue
            if right["source_atom_ids"][0] != left["carrier_atom_id"]:
                continue
            if endpoint_gap(left, right) <= 1e-7:
                union.union(arc_patch[left["arc_id"]], arc_patch[right["arc_id"]])

    groups = defaultdict(list)
    for arc in arcs:
        if arc["arc_id"] in matched:
            continue
        groups[union.find(arc_patch[arc["arc_id"]])].append(arc)

    bad = {}
    for root, component_arcs in groups.items():
        outgoing = defaultdict(list)
        incoming = defaultdict(list)
        for arc in component_arcs:
            if arc["start_vertex_id"] is not None:
                outgoing[arc["start_vertex_id"]].append(arc)
            if arc["end_vertex_id"] is not None:
                incoming[arc["end_vertex_id"]].append(arc)
        for vertex_id in set(outgoing) | set(incoming):
            out_degree = len(outgoing[vertex_id])
            in_degree = len(incoming[vertex_id])
            if out_degree != 1 or in_degree != 1:
                bad[vertex_id] = {
                    "root": root,
                    "outgoing": outgoing[vertex_id],
                    "incoming": incoming[vertex_id],
                    "out_degree": out_degree,
                    "in_degree": in_degree,
                }
    return {
        "patches": patches,
        "arcs": arcs,
        "junctions": surface["junctions"],
        "arc_patch": arc_patch,
        "self_arcs": self_arcs,
        "matched": matched,
        "endpoint_pairs": endpoint_pairs,
        "union": union,
        "groups": groups,
        "bad": bad,
    }


def seam_audit(surface, graph, atoms):
    by_pair = defaultdict(lambda: defaultdict(list))
    for arc in graph["self_arcs"]:
        key = tuple(sorted((arc["carrier_atom_id"], arc["source_atom_ids"][0])))
        by_pair[key][arc["carrier_atom_id"]].append(arc)

    internal_arcs = set()
    union = UnionFind(len(graph["patches"]))
    pair_reports = []
    for pair, sides in by_pair.items():
        side_values = list(sides.values())
        if len(side_values) != 2:
            pair_reports.append({"pair": pair, "status": "missing_counterpart"})
            continue
        first, second = side_values
        first_intervals = merge_intervals(
            interval for arc in first for interval in arc_interval(arc, atoms)
        )
        second_intervals = merge_intervals(
            interval for arc in second for interval in arc_interval(arc, atoms)
        )
        overlap = intersection_intervals(first_intervals, second_intervals)
        endpoints = sorted({
            endpoint
            for arc in first + second
            for interval in arc_interval(arc, atoms)
            for endpoint in interval
        })
        coverage_mismatch = 0
        for index, left in enumerate(endpoints):
            right = endpoints[index + 1] if index + 1 < len(endpoints) else left + TWO_PI
            midpoint = ((left + right) / 2.0) % TWO_PI
            first_cover = any(interval[0] - 1e-8 <= midpoint <= interval[1] + 1e-8 for interval in first_intervals)
            second_cover = any(interval[0] - 1e-8 <= midpoint <= interval[1] + 1e-8 for interval in second_intervals)
            coverage_mismatch += first_cover != second_cover
        equal_coverage = coverage_mismatch == 0
        if equal_coverage:
            for arc in first + second:
                internal_arcs.add(arc["arc_id"])
            for left in first:
                for right in second:
                    if intersection_intervals(arc_interval(left, atoms), arc_interval(right, atoms)):
                        union.union(graph["arc_patch"][left["arc_id"]], graph["arc_patch"][right["arc_id"]])
        pair_reports.append({
            "pair": pair,
            "left_arcs": len(first),
            "right_arcs": len(second),
            "left_right_count_equal": len(first) == len(second),
            "coverage_equal": equal_coverage,
            "geometric_seam_pieces": len(overlap),
            "current_endpoint_pairs": sum(
                1
                for left in first
                for right in second
                if endpoint_gap(left, right) <= 1e-7
            ),
        })

    return {
        "pair_reports": pair_reports,
        "internal_arcs": internal_arcs,
        "union": union,
    }


def corrected_boundary_graph(surface, graph, seams):
    patches = graph["patches"]
    corrected_union = seams["union"]
    remaining = [arc for arc in graph["arcs"] if arc["arc_id"] not in seams["internal_arcs"]]
    outgoing = defaultdict(list)
    incoming = defaultdict(list)
    for arc in remaining:
        if arc["start_vertex_id"] is not None:
            outgoing[arc["start_vertex_id"]].append(arc)
        if arc["end_vertex_id"] is not None:
            incoming[arc["end_vertex_id"]].append(arc)
    bad = {
        vertex: (len(outgoing[vertex]), len(incoming[vertex]))
        for vertex in set(outgoing) | set(incoming)
        if len(outgoing[vertex]) != 1 or len(incoming[vertex]) != 1
    }
    loops = 0
    visited = set()
    by_start = {vertex: values[0] for vertex, values in outgoing.items() if len(values) == 1}
    for arc in remaining:
        if arc["arc_id"] in visited or arc["start_vertex_id"] is None:
            continue
        loops += 1
        current = arc
        while current["arc_id"] not in visited:
            visited.add(current["arc_id"])
            current = by_start.get(current["end_vertex_id"])
            if current is None:
                break
    components = len({corrected_union.find(index) for index in range(len(patches))})
    return {
        "remaining_arcs": len(remaining),
        "boundary_arc_count": sum(bool(arc["source_contact_ids"]) for arc in remaining),
        "boundary_degree_defects": bad,
        "boundary_loops": loops if not bad else None,
        "components": components,
    }


def classify(surface, graph, atoms, selected_partners):
    junctions = {junction["vertex_id"]: junction for junction in surface["junctions"]}
    arcs = {arc["arc_id"]: arc for arc in surface["arcs"]}
    patch_by_arc = graph["arc_patch"]
    categories = Counter()
    representatives = {}
    max_residual = 0.0
    strict_occlusion = 0
    partner_misses = 0
    records = []
    for vertex_id, defect in graph["bad"].items():
        junction = junctions[vertex_id]
        incident = {
            arc_id
            for arc_id in junction["incident_arc_ids"]
            if arc_id in arcs
        }
        sphere_ids = set(junction["source_carrier_atom_ids"] or (junction["carrier_atom_id"],))
        contact_ids = set()
        for arc_id in incident:
            arc = arcs[arc_id]
            sphere_ids.update(arc["source_atom_ids"])
            contact_ids.update(arc["source_contact_ids"])
        arc_validations = {
            arc_id: validate_arc(arcs[arc_id], atoms, selected_partners)
            for arc_id in sorted(incident)
        }
        residuals = {
            atom_id: abs(distance(junction["xyz"], atoms[atom_id]["center"]) - atoms[atom_id]["radius"])
            for atom_id in sphere_ids
            if atom_id in atoms
        }
        max_residual = max(max_residual, max(residuals.values(), default=0.0))
        carrier_ids = set(junction["source_carrier_atom_ids"] or (junction["carrier_atom_id"],))
        for carrier_id in carrier_ids:
            carrier = atoms[carrier_id]
            if any(
                other["group"] == carrier["group"]
                and other["id"] != carrier_id
                and distance(junction["xyz"], other["center"]) < other["radius"] - 1e-8
                for other in atoms.values()
            ):
                strict_occlusion += 1
        if not any(
            distance(junction["xyz"], partner["center"]) <= partner["radius"] + 1e-7
            for carrier_id in carrier_ids
            for partner in selected_partners.get(carrier_id, ())
        ):
            partner_misses += 1

        sphere_count = len(sphere_ids)
        if sphere_count == 2 and len(contact_ids) == 0:
            category = "two-sphere-seam-unmatched"
        elif sphere_count == 3 and len(contact_ids) == 0:
            category = "three-sphere-same-group-junction"
        elif sphere_count == 3 and contact_ids:
            category = "three-sphere-mixed-boundary-junction"
        else:
            category = "other-or-unresolved"
        record = {
            "vertex_id": vertex_id,
            "xyz": junction["xyz"],
            "incident_carrier_atom_ids": sorted(carrier_ids),
            "incident_patch_ids": sorted({surface_patch for surface_patch in (
                graph["patches"][patch_by_arc[arc_id]]["patch_id"]
                for arc_id in incident
                if arc_id in patch_by_arc
            )}),
            "incident_arc_ids": sorted(incident),
            "out_degree": defect["out_degree"],
            "in_degree": defect["in_degree"],
            "sphere_ids": sorted(sphere_ids),
            "sphere_residuals_A": residuals,
            "selected_contact_ids": sorted(contact_ids),
            "stored_source_circle_ids": sorted(junction["source_circle_ids"]),
            "all_incident_arcs_are_nonzero": all(abs(arcs[arc_id]["angular_span"]) > 1e-10 for arc_id in incident),
            "true_multi_sphere_intersection": max(residuals.values(), default=math.inf) <= 1e-8,
            "should_pair_across_carriers": len(carrier_ids) > 1,
            "arc_validations": arc_validations,
            "invalid_incident_arc_count": sum(not value["valid"] for value in arc_validations.values()),
        }
        categories[category] += 1
        records.append(record)
        representatives.setdefault(category, record)
    coincident_pairs = sum(
        distance(left["xyz"], right["xyz"]) <= 1e-7
        for index, left in enumerate(surface["junctions"])
        for right in surface["junctions"][index + 1:]
    )
    return {
        "categories": dict(categories),
        "vertices": sorted(records, key=lambda value: value["vertex_id"]),
        "representatives": representatives,
        "max_sphere_residual_A": max_residual,
        "strict_self_occluded_carrier_incidences": strict_occlusion,
        "vertices_without_selected_partner_membership": partner_misses,
        "coincident_stored_junction_pairs": coincident_pairs,
    }


def audit_side(side, atoms, by_index, selected_partners):
    surface = json.loads((SURFACE_DIR / f"surface_{side}.json").read_text(encoding="utf-8"))
    graph = current_graph(surface)
    seams = seam_audit(surface, graph, atoms)
    corrected = corrected_boundary_graph(surface, graph, seams)
    classified = classify(surface, graph, atoms, selected_partners)
    pair_reports = seams["pair_reports"]
    return {
        "side": side,
        "stored_summary": surface["summary"],
        "diagnostic_vertex_count": len(graph["bad"]),
        "current_endpoint_matched_seam_pairs": graph["endpoint_pairs"],
        "current_matched_seam_arcs": len(graph["matched"]),
        "self_seam_arc_count": len(graph["self_arcs"]),
        "sphere_pair_groups": len(pair_reports),
        "sphere_pair_groups_with_unequal_segmentation": sum(not item.get("left_right_count_equal", True) for item in pair_reports),
        "sphere_pair_groups_with_missing_counterpart": sum(item.get("status") == "missing_counterpart" for item in pair_reports),
        "sphere_pair_groups_with_coverage_mismatch": sum(not item.get("coverage_equal", True) for item in pair_reports if "coverage_equal" in item),
        "geometric_internal_seam_arcs": len(seams["internal_arcs"]),
        "geometric_internal_seam_pieces": sum(item.get("geometric_seam_pieces", 0) for item in pair_reports),
        "endpoint_matched_pairs_with_geometric_overlap": sum(
            item.get("current_endpoint_pairs", 0) for item in pair_reports
        ),
        "corrected_boundary": corrected,
        "classification": classified,
    }


def main():
    atoms, by_index = load_atoms()
    selected_partners = load_selected_partners(by_index)
    report = {
        "definition": "read-only audit of the persisted selected partner-restricted union-of-balls surfaces",
        "tolerances": {"endpoint_A": 1e-7, "sphere_residual_A": 1e-8},
        "A": audit_side("A", atoms, by_index, selected_partners),
        "B": audit_side("B", atoms, by_index, selected_partners),
    }
    output = SURFACE_DIR / "topology_audit.json"
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
