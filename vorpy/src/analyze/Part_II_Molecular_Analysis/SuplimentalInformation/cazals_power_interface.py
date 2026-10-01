"""
cazals_power_interface.py

Standalone Cazals-style bicolor power-interface analysis using GUDHI.

Repository location:
    vorpy/src/analyze/Part_II_Molecular_Analysis/SuplimentalInformation/cazals_power_interface.py

Inputs
------
1. Existing NPZ:
       points, radii, group_a, group_b

2. PDB directly:
       --group-a-chain A --group-b-chain B

For PDB input the script uses:
    vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation.cazals_pdb.prepare_cazals_pdb

which assigns Cazals/Chothia base radii, expands each by the 1.4 A probe,
and sends R_i = r_i + 1.4 to the weighted power construction.

Examples
--------
NPZ:
    python -m vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation.cazals_power_interface system.npz

2KAI:
    python -m vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation.cazals_power_interface \
        vorpy/data/2kai.pdb \
        --group-a-chain A \
        --group-b-chain B \
        --strict-heavy-fallbacks \
        --audit-csv 2kai_cazals_atoms.csv \
        --npz-out 2kai_cazals_input.npz \
        --edge-csv 2kai_interface_edges.csv

Important
---------
This script computes the geometric bicolor power-interface statistic. Exact
reproduction of a published Cazals value additionally requires matching the
paper's alpha restriction and any additional interface-selection condition.
Those restrictions should be validated explicitly rather than silently
approximated here.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

try:
    import gudhi
except ImportError as exc:
    raise ImportError(
        "GUDHI is required. Install it in the active environment with:\n"
        "    pip install gudhi"
    ) from exc

from vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation.cazals_pdb import (
    export_audit_csv,
    export_npz,
    prepare_cazals_pdb,
)


@dataclass
class PowerVertex:
    tetrahedron: tuple[int, int, int, int]
    position: np.ndarray
    power_value: float
    residual: float


@dataclass
class InterfaceEdgeMeasurement:
    triangle: tuple[int, int, int]
    facet_1: tuple[int, int]
    facet_2: tuple[int, int]
    tetrahedron_1: tuple[int, int, int, int]
    tetrahedron_2: tuple[int, int, int, int]
    vertex_1: np.ndarray
    vertex_2: np.ndarray
    length: float
    beta_rad: float
    beta_deg: float
    abs_beta_rad: float
    abs_beta_deg: float
    length_times_beta: float
    length_times_abs_beta: float


@dataclass
class InterfaceResult:
    num_input_atoms: int
    num_group_a: int
    num_group_b: int
    active_vertices: set[int]
    hidden_vertices: set[int]
    regular_edges: set[tuple[int, int]]
    regular_triangles: set[tuple[int, int, int]]
    regular_tetrahedra: set[tuple[int, int, int, int]]
    bicolor_edges: set[tuple[int, int]]
    power_vertices: dict[tuple[int, int, int, int], PowerVertex]
    measurements: list[InterfaceEdgeMeasurement]
    skipped_boundary_triangles: int
    skipped_nonmanifold_triangles: int
    skipped_wrong_bicolor_count: int
    skipped_degenerate_edges: int
    skipped_power_vertices: int
    total_interface_edge_length: float
    signed_integrated_angle: float
    unsigned_integrated_angle: float
    mean_signed_angle_rad: float
    mean_signed_angle_deg: float
    mean_unsigned_angle_rad: float
    mean_unsigned_angle_deg: float


def canonical_simplex(indices: Iterable[int]) -> tuple[int, ...]:
    return tuple(sorted(int(i) for i in indices))


def normalize(v: np.ndarray, tol: float = 1e-12) -> np.ndarray:
    n = float(np.linalg.norm(v))
    if n <= tol:
        raise ValueError("Cannot normalize a zero-length vector.")
    return v / n


def validate_inputs(points, radii, group_a, group_b):
    points = np.asarray(points, dtype=float)
    radii = np.asarray(radii, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError(f"points must have shape (N,3), got {points.shape}")
    n = len(points)
    if radii.shape != (n,):
        raise ValueError(f"radii must have shape ({n},), got {radii.shape}")
    if not np.all(np.isfinite(points)) or not np.all(np.isfinite(radii)):
        raise ValueError("points/radii contain non-finite values")
    if np.any(radii <= 0):
        raise ValueError("All radii must be positive")

    a = {int(i) for i in group_a}
    b = {int(i) for i in group_b}
    if not a or not b:
        raise ValueError("Both groups must be non-empty")
    if a & b:
        raise ValueError(f"Groups overlap at indices {sorted(a & b)[:20]}")
    bad = sorted(i for i in a | b if i < 0 or i >= n)
    if bad:
        raise IndexError(f"Group indices outside [0,{n-1}]: {bad[:20]}")
    return points, radii, a, b


def build_regular_complex(points: np.ndarray, radii: np.ndarray):
    """
    Build a weighted GUDHI alpha complex and extract its regular simplices.

    GUDHI API compatibility is checked explicitly because weighted-alpha support
    differs between releases.
    """
    weights = np.asarray(radii, dtype=float) ** 2
    try:
        ac = gudhi.AlphaComplex(
            points=points.tolist(),
            weights=weights.tolist(),
            precision="safe",
        )
    except TypeError as exc:
        raise RuntimeError(
            "This installed GUDHI version does not accept `weights=` in "
            "gudhi.AlphaComplex. Do not fall back to an unweighted alpha "
            "complex: that would change the Cazals geometry. Install/use a "
            "GUDHI release/API that supports weighted alpha complexes, or "
            "adapt this function to GUDHI's weighted DelaunayComplex API "
            "after checking that version's documentation."
        ) from exc

    st = ac.create_simplex_tree()
    vertices, edges, triangles, tetrahedra = set(), set(), set(), set()

    for simplex, _filtration in st.get_simplices():
        s = canonical_simplex(simplex)
        if len(s) == 1:
            vertices.add(s[0])
        elif len(s) == 2:
            edges.add(s)
        elif len(s) == 3:
            triangles.add(s)
        elif len(s) == 4:
            tetrahedra.add(s)
        elif len(s) > 4:
            raise RuntimeError(f"Unexpected >3D simplex: {s}")

    return st, vertices, edges, triangles, tetrahedra


def calculate_power_vertex(tetrahedron, points, weights, condition_limit=1e12):
    tet = canonical_simplex(tetrahedron)
    i = tet[0]
    p0, w0 = points[i], weights[i]
    A, b = [], []

    for j in tet[1:]:
        pj, wj = points[j], weights[j]
        A.append(2.0 * (pj - p0))
        b.append(np.dot(pj, pj) - np.dot(p0, p0) - wj + w0)

    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float)
    condition = np.linalg.cond(A)
    if not np.isfinite(condition) or condition > condition_limit:
        raise np.linalg.LinAlgError(
            f"Ill-conditioned tetrahedron {tet}, condition={condition:.6g}"
        )

    x = np.linalg.solve(A, b)
    values = np.asarray(
        [np.dot(x - points[j], x - points[j]) - weights[j] for j in tet]
    )
    value = float(np.mean(values))
    residual = float(np.max(np.abs(values - value)))
    return PowerVertex(tet, x, value, residual)


def is_bicolor_edge(edge, group_a, group_b):
    i, j = edge
    return ((i in group_a and j in group_b) or
            (j in group_a and i in group_b))


def oriented_bicolor_facet_normal(edge, points, group_a, group_b):
    i, j = edge
    if i in group_a and j in group_b:
        a, b = i, j
    elif j in group_a and i in group_b:
        a, b = j, i
    else:
        raise ValueError(f"{edge} is not bicolor")
    return normalize(points[b] - points[a])


def bicolor_edges_of_triangle(triangle, group_a, group_b):
    i, j, k = triangle
    edges = [
        canonical_simplex((i, j)),
        canonical_simplex((i, k)),
        canonical_simplex((j, k)),
    ]
    return [e for e in edges if is_bicolor_edge(e, group_a, group_b)]


def build_triangle_to_tetrahedra(tetrahedra):
    result = defaultdict(list)
    for i, j, k, l in tetrahedra:
        for face in (
            canonical_simplex((i, j, k)),
            canonical_simplex((i, j, l)),
            canonical_simplex((i, k, l)),
            canonical_simplex((j, k, l)),
        ):
            result[face].append((i, j, k, l))
    return dict(result)


def signed_interface_angle(tangent, normal_1, normal_2):
    t = normalize(tangent)
    n1 = normalize(normal_1)
    n2 = normalize(normal_2)
    sin_term = float(np.dot(t, np.cross(n1, n2)))
    cos_term = float(np.clip(np.dot(n1, n2), -1.0, 1.0))
    return float(np.arctan2(sin_term, cos_term))


def orient_power_edge_tangent(triangle, v1, v2, points):
    i, j, k = canonical_simplex(triangle)
    tri_normal = normalize(np.cross(points[j] - points[i], points[k] - points[i]))
    edge = v2 - v1
    if np.dot(edge, tri_normal) < 0:
        return -edge, v2, v1
    return edge, v1, v2


def analyze_bicolor_power_interface(
    points,
    radii,
    group_a,
    group_b,
    *,
    power_residual_tolerance=1e-8,
    length_tolerance=1e-10,
    verbose=True,
):
    points, radii, group_a, group_b = validate_inputs(
        points, radii, group_a, group_b
    )
    weights = radii ** 2

    st, vertices, edges, triangles, tetrahedra = build_regular_complex(
        points, radii
    )
    del st

    hidden = set(range(len(points))) - vertices
    bicolor = {e for e in edges if is_bicolor_edge(e, group_a, group_b)}

    power_vertices = {}
    skipped_pv = 0
    for tet in tetrahedra:
        try:
            pv = calculate_power_vertex(tet, points, weights)
        except np.linalg.LinAlgError:
            skipped_pv += 1
            continue
        if pv.residual > power_residual_tolerance:
            skipped_pv += 1
            continue
        power_vertices[tet] = pv

    tri_to_tets = build_triangle_to_tetrahedra(tetrahedra)

    measurements = []
    skipped_boundary = skipped_nonmanifold = 0
    skipped_wrong_bicolor = skipped_degenerate = 0

    for tri in sorted(triangles):
        bi_edges = bicolor_edges_of_triangle(tri, group_a, group_b)
        if len(bi_edges) != 2:
            skipped_wrong_bicolor += 1
            continue

        incident = tri_to_tets.get(tri, [])
        if len(incident) == 1:
            skipped_boundary += 1
            continue
        if len(incident) != 2:
            skipped_nonmanifold += 1
            continue

        tet1, tet2 = incident
        if tet1 not in power_vertices or tet2 not in power_vertices:
            skipped_degenerate += 1
            continue

        raw1 = power_vertices[tet1].position
        raw2 = power_vertices[tet2].position
        try:
            tangent, v1, v2 = orient_power_edge_tangent(
                tri, raw1, raw2, points
            )
        except ValueError:
            skipped_degenerate += 1
            continue

        length = float(np.linalg.norm(tangent))
        if length <= length_tolerance:
            skipped_degenerate += 1
            continue

        facet1, facet2 = sorted(bi_edges)
        try:
            n1 = oriented_bicolor_facet_normal(
                facet1, points, group_a, group_b
            )
            n2 = oriented_bicolor_facet_normal(
                facet2, points, group_a, group_b
            )
        except ValueError:
            skipped_degenerate += 1
            continue

        beta = signed_interface_angle(tangent, n1, n2)
        beta_deg = float(np.degrees(beta))
        measurements.append(
            InterfaceEdgeMeasurement(
                triangle=tri,
                facet_1=facet1,
                facet_2=facet2,
                tetrahedron_1=tet1,
                tetrahedron_2=tet2,
                vertex_1=np.asarray(v1),
                vertex_2=np.asarray(v2),
                length=length,
                beta_rad=beta,
                beta_deg=beta_deg,
                abs_beta_rad=abs(beta),
                abs_beta_deg=abs(beta_deg),
                length_times_beta=length * beta,
                length_times_abs_beta=length * abs(beta),
            )
        )

    total_length = float(sum(m.length for m in measurements))
    signed = float(sum(m.length_times_beta for m in measurements))
    unsigned = float(sum(m.length_times_abs_beta for m in measurements))

    if total_length:
        mean_signed = signed / total_length
        mean_unsigned = unsigned / total_length
    else:
        mean_signed = mean_unsigned = float("nan")

    result = InterfaceResult(
        num_input_atoms=len(points),
        num_group_a=len(group_a),
        num_group_b=len(group_b),
        active_vertices=vertices,
        hidden_vertices=hidden,
        regular_edges=edges,
        regular_triangles=triangles,
        regular_tetrahedra=tetrahedra,
        bicolor_edges=bicolor,
        power_vertices=power_vertices,
        measurements=measurements,
        skipped_boundary_triangles=skipped_boundary,
        skipped_nonmanifold_triangles=skipped_nonmanifold,
        skipped_wrong_bicolor_count=skipped_wrong_bicolor,
        skipped_degenerate_edges=skipped_degenerate,
        skipped_power_vertices=skipped_pv,
        total_interface_edge_length=total_length,
        signed_integrated_angle=signed,
        unsigned_integrated_angle=unsigned,
        mean_signed_angle_rad=mean_signed,
        mean_signed_angle_deg=float(np.degrees(mean_signed)),
        mean_unsigned_angle_rad=mean_unsigned,
        mean_unsigned_angle_deg=float(np.degrees(mean_unsigned)),
    )
    if verbose:
        print_summary(result)
    return result


def print_summary(r):
    pos = sum(m.beta_rad > 0 for m in r.measurements)
    neg = sum(m.beta_rad < 0 for m in r.measurements)
    zero = len(r.measurements) - pos - neg
    print()
    print("=" * 68)
    print("BICOLOR POWER-INTERFACE CURVATURE")
    print("=" * 68)
    print(f"Input atoms:                     {r.num_input_atoms:,}")
    print(f"Group A / B:                     {r.num_group_a:,} / {r.num_group_b:,}")
    print(f"Hidden weighted vertices:        {len(r.hidden_vertices):,}")
    print(f"Regular edges:                   {len(r.regular_edges):,}")
    print(f"Regular triangles:               {len(r.regular_triangles):,}")
    print(f"Regular tetrahedra:              {len(r.regular_tetrahedra):,}")
    print(f"Bicolor regular edges/facets:    {len(r.bicolor_edges):,}")
    print(f"Finite interface edges:          {len(r.measurements):,}")
    print("-" * 68)
    print(f"Total interface edge length:     {r.total_interface_edge_length:.10f} A")
    print(f"Sum(length * beta):              {r.signed_integrated_angle:.10f} A rad")
    print(f"Sum(length * |beta|):            {r.unsigned_integrated_angle:.10f} A rad")
    print(f"Mean signed angle:               {r.mean_signed_angle_deg:.10f} deg")
    print(f"Mean unsigned angle:             {r.mean_unsigned_angle_deg:.10f} deg")
    print(f"Positive / negative / zero:      {pos:,} / {neg:,} / {zero:,}")
    print("-" * 68)
    print(f"Boundary triangles skipped:      {r.skipped_boundary_triangles:,}")
    print(f"Nonmanifold triangles skipped:   {r.skipped_nonmanifold_triangles:,}")
    print(f"Degenerate edges skipped:        {r.skipped_degenerate_edges:,}")
    print(f"Power vertices unresolved:       {r.skipped_power_vertices:,}")
    print("=" * 68)
    print()


def export_measurements_csv(result, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "triangle_i", "triangle_j", "triangle_k",
        "facet_1_i", "facet_1_j", "facet_2_i", "facet_2_j",
        "tetrahedron_1", "tetrahedron_2",
        "v1_x", "v1_y", "v1_z", "v2_x", "v2_y", "v2_z",
        "length_A", "beta_rad", "beta_deg", "abs_beta_deg",
        "length_times_beta_A_rad", "length_times_abs_beta_A_rad",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for m in result.measurements:
            writer.writerow(
                {
                    "triangle_i": m.triangle[0],
                    "triangle_j": m.triangle[1],
                    "triangle_k": m.triangle[2],
                    "facet_1_i": m.facet_1[0],
                    "facet_1_j": m.facet_1[1],
                    "facet_2_i": m.facet_2[0],
                    "facet_2_j": m.facet_2[1],
                    "tetrahedron_1": " ".join(map(str, m.tetrahedron_1)),
                    "tetrahedron_2": " ".join(map(str, m.tetrahedron_2)),
                    "v1_x": m.vertex_1[0],
                    "v1_y": m.vertex_1[1],
                    "v1_z": m.vertex_1[2],
                    "v2_x": m.vertex_2[0],
                    "v2_y": m.vertex_2[1],
                    "v2_z": m.vertex_2[2],
                    "length_A": m.length,
                    "beta_rad": m.beta_rad,
                    "beta_deg": m.beta_deg,
                    "abs_beta_deg": m.abs_beta_deg,
                    "length_times_beta_A_rad": m.length_times_beta,
                    "length_times_abs_beta_A_rad": m.length_times_abs_beta,
                }
            )
    return path


def load_npz(path):
    with np.load(path) as data:
        required = {"points", "radii", "group_a", "group_b"}
        missing = required - set(data.files)
        if missing:
            raise KeyError(f"{path} missing arrays: {sorted(missing)}")
        return (
            np.asarray(data["points"], dtype=float),
            np.asarray(data["radii"], dtype=float),
            np.asarray(data["group_a"], dtype=int),
            np.asarray(data["group_b"], dtype=int),
        )


def _chain_set(values):
    result = set()
    for value in values or []:
        for piece in value.split(","):
            piece = piece.strip()
            if piece:
                result.add(piece)
    return result


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("input", type=Path, help=".npz or .pdb input")
    p.add_argument("--group-a-chain", action="append", default=[])
    p.add_argument("--group-b-chain", action="append", default=[])
    p.add_argument("--include-hydrogens", action="store_true")
    p.add_argument("--include-hetatm", action="store_true")
    p.add_argument("--strict-heavy-fallbacks", action="store_true")
    p.add_argument("--probe-radius", type=float, default=1.40)
    p.add_argument("--audit-csv", type=Path)
    p.add_argument("--npz-out", type=Path)
    p.add_argument("--edge-csv", type=Path)
    return p.parse_args()


def main():
    args = parse_args()
    suffix = args.input.suffix.lower()

    if suffix == ".npz":
        points, radii, group_a, group_b = load_npz(args.input)

    elif suffix == ".pdb":
        a = _chain_set(args.group_a_chain)
        b = _chain_set(args.group_b_chain)
        if not a or not b:
            raise SystemExit(
                "PDB input requires --group-a-chain and --group-b-chain."
            )

        prepared = prepare_cazals_pdb(
            args.input,
            a,
            b,
            probe_radius=args.probe_radius,
            include_hydrogens=args.include_hydrogens,
            include_hetatm=args.include_hetatm,
            strict_heavy_fallbacks=args.strict_heavy_fallbacks,
            verbose=True,
        )
        points = prepared.points
        radii = prepared.expanded_radii
        group_a = prepared.group_a
        group_b = prepared.group_b

        if args.audit_csv:
            print(f"Audit CSV: {export_audit_csv(prepared, args.audit_csv)}")
        if args.npz_out:
            print(f"NPZ: {export_npz(prepared, args.npz_out)}")
    else:
        raise SystemExit("Input must be .pdb or .npz")

    result = analyze_bicolor_power_interface(
        points, radii, group_a, group_b, verbose=True
    )

    if args.edge_csv:
        print(f"Interface-edge CSV: {export_measurements_csv(result, args.edge_csv)}")


if __name__ == "__main__":
    main()
