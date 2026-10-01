"""
Diagnostic comparison of unrestricted and Cazals-style 2KAI interface curvature.

Place this file beside:
    cazals_power_interface.py

It deliberately reuses the existing power/regular-complex implementation and
does not alter the validated baseline script.

Reported variants
-----------------
1. Existing unrestricted normal/tangent statistic (baseline)
2. Unrestricted Cazals triangle-angle/sign statistic
3. alpha=0 selected Cazals statistic

Condition beta (M=5) is NOT implemented here. It requires the exact largest
orthogonal-ball construction for each retained bicolor facet and should not be
approximated merely to move the result toward the published ~17 degrees.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation import (
    cazals_power_interface as base,
)
from vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation.cazals_pdb import (
    prepare_cazals_pdb,
)


ALPHA_TOL = 1.0e-12


@dataclass
class DiagnosticMeasurement:
    triangle: tuple[int, int, int]
    singleton: int
    singleton_group: str
    facet_1: tuple[int, int]
    facet_2: tuple[int, int]
    length: float
    beta_rad: float
    beta_deg: float
    contribution: float
    facet_1_filtration: float
    facet_2_filtration: float


def _chain_set(values):
    result = set()
    for value in values or []:
        for piece in value.split(","):
            piece = piece.strip()
            if piece:
                result.add(piece)
    return result


def _simplex_filtrations(st):
    """Return canonical simplex -> GUDHI filtration value."""
    values = {}
    for simplex, filtration in st.get_filtration():
        values[base.canonical_simplex(simplex)] = float(filtration)
    return values


def _singleton_and_sign(triangle, group_a, group_b):
    """
    Cazals convention used for this diagnostic:

        A-B-B -> positive, angle at A
        A-A-B -> negative, angle at B

    Returns (singleton_index, sign, singleton_group).
    """
    tri = tuple(triangle)
    a = [i for i in tri if i in group_a]
    b = [i for i in tri if i in group_b]

    if len(a) == 1 and len(b) == 2:
        return a[0], +1.0, "A"
    if len(a) == 2 and len(b) == 1:
        return b[0], -1.0, "B"
    raise ValueError(f"Triangle {tri} is not a 1:2 bicolor triangle")


def _triangle_angle_at(singleton, triangle, points):
    others = [i for i in triangle if i != singleton]
    if len(others) != 2:
        raise ValueError(f"Bad triangle/singleton: {triangle}, {singleton}")

    u = points[others[0]] - points[singleton]
    v = points[others[1]] - points[singleton]
    nu = float(np.linalg.norm(u))
    nv = float(np.linalg.norm(v))
    if nu <= 1e-14 or nv <= 1e-14:
        raise ValueError("Degenerate Delaunay triangle")

    cosine = float(np.clip(np.dot(u, v) / (nu * nv), -1.0, 1.0))
    return float(np.arccos(cosine))


def _calculate_power_vertices(tetrahedra, points, weights, residual_tol=1e-8):
    power_vertices = {}
    unresolved = 0
    for tet in tetrahedra:
        try:
            pv = base.calculate_power_vertex(tet, points, weights)
        except np.linalg.LinAlgError:
            unresolved += 1
            continue
        if pv.residual > residual_tol:
            unresolved += 1
            continue
        power_vertices[tet] = pv
    return power_vertices, unresolved


def _measure(
    points,
    radii,
    group_a,
    group_b,
    *,
    alpha_zero=False,
    alpha_tol=ALPHA_TOL,
):
    points, radii, group_a, group_b = base.validate_inputs(
        points, radii, group_a, group_b
    )
    weights = radii ** 2

    st, vertices, edges, triangles, tetrahedra = base.build_regular_complex(
        points, radii
    )
    filtration = _simplex_filtrations(st)

    bicolor = {e for e in edges if base.is_bicolor_edge(e, group_a, group_b)}
    if alpha_zero:
        retained_bicolor = {
            e for e in bicolor
            if filtration.get(e, float("inf")) <= alpha_tol
        }
    else:
        retained_bicolor = set(bicolor)

    power_vertices, unresolved_pv = _calculate_power_vertices(
        tetrahedra, points, weights
    )
    tri_to_tets = base.build_triangle_to_tetrahedra(tetrahedra)

    measurements = []
    skipped_boundary = 0
    skipped_nonmanifold = 0
    skipped_wrong_selected_count = 0
    skipped_degenerate = 0

    for tri in sorted(triangles):
        all_bicolor = base.bicolor_edges_of_triangle(tri, group_a, group_b)
        if len(all_bicolor) != 2:
            continue

        # The interface is the set of retained bicolor Voronoi facets.
        # An interior interface edge contributes only when BOTH of the two
        # bicolor facets incident to this regular triangle survive selection.
        selected = [e for e in all_bicolor if e in retained_bicolor]
        if len(selected) != 2:
            skipped_wrong_selected_count += 1
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

        v1 = power_vertices[tet1].position
        v2 = power_vertices[tet2].position
        length = float(np.linalg.norm(v2 - v1))
        if length <= 1e-10:
            skipped_degenerate += 1
            continue

        try:
            singleton, sign, singleton_group = _singleton_and_sign(
                tri, group_a, group_b
            )
            magnitude = _triangle_angle_at(singleton, tri, points)
        except ValueError:
            skipped_degenerate += 1
            continue

        beta = sign * magnitude
        facet1, facet2 = sorted(selected)
        measurements.append(
            DiagnosticMeasurement(
                triangle=tri,
                singleton=singleton,
                singleton_group=singleton_group,
                facet_1=facet1,
                facet_2=facet2,
                length=length,
                beta_rad=beta,
                beta_deg=float(np.degrees(beta)),
                contribution=length * beta,
                facet_1_filtration=filtration[facet1],
                facet_2_filtration=filtration[facet2],
            )
        )

    return {
        "vertices": vertices,
        "edges": edges,
        "triangles": triangles,
        "tetrahedra": tetrahedra,
        "all_bicolor": bicolor,
        "retained_bicolor": retained_bicolor,
        "measurements": measurements,
        "unresolved_power_vertices": unresolved_pv,
        "skipped_boundary": skipped_boundary,
        "skipped_nonmanifold": skipped_nonmanifold,
        "skipped_selection": skipped_wrong_selected_count,
        "skipped_degenerate": skipped_degenerate,
        "filtration": filtration,
    }


def _summary(name, data):
    ms = data["measurements"]
    total_length = float(sum(m.length for m in ms))
    positive = float(sum(m.contribution for m in ms if m.beta_rad > 0))
    negative = float(sum(m.contribution for m in ms if m.beta_rad < 0))
    signed = positive + negative
    unsigned = float(sum(abs(m.contribution) for m in ms))

    signed_deg = (
        float(np.degrees(signed / total_length)) if total_length else float("nan")
    )
    unsigned_deg = (
        float(np.degrees(unsigned / total_length)) if total_length else float("nan")
    )

    pos_n = sum(m.beta_rad > 0 for m in ms)
    neg_n = sum(m.beta_rad < 0 for m in ms)
    zero_n = len(ms) - pos_n - neg_n

    print()
    print("=" * 72)
    print(name)
    print("=" * 72)
    print(f"Bicolor regular facets total:      {len(data['all_bicolor']):,}")
    print(f"Bicolor facets retained:           {len(data['retained_bicolor']):,}")
    print(f"Interior interface edges:          {len(ms):,}")
    print(f"Positive / negative / zero:        {pos_n:,} / {neg_n:,} / {zero_n:,}")
    print("-" * 72)
    print(f"Total interface edge length:       {total_length:.10f} A")
    print(f"Positive sum(length*beta):         {positive:.10f} A rad")
    print(f"Negative sum(length*beta):         {negative:.10f} A rad")
    print(f"Signed sum(length*beta):           {signed:.10f} A rad")
    print(f"Unsigned sum(length*|beta|):       {unsigned:.10f} A rad")
    print(f"s_H signed:                        {signed_deg:.10f} deg")
    print(f"s_H unsigned:                      {unsigned_deg:.10f} deg")
    print("-" * 72)
    print(f"Selection-rejected triangles:      {data['skipped_selection']:,}")
    print(f"Boundary triangles skipped:        {data['skipped_boundary']:,}")
    print(f"Nonmanifold triangles skipped:     {data['skipped_nonmanifold']:,}")
    print(f"Degenerate/unresolved skipped:     {data['skipped_degenerate']:,}")
    print(f"Power vertices unresolved:         {data['unresolved_power_vertices']:,}")
    print("=" * 72)

    return {
        "total_length": total_length,
        "positive": positive,
        "negative": negative,
        "signed": signed,
        "unsigned": unsigned,
        "signed_deg": signed_deg,
        "unsigned_deg": unsigned_deg,
    }


def export_csv(data, path):
    import csv

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "triangle_i", "triangle_j", "triangle_k",
        "singleton", "singleton_group",
        "facet_1_i", "facet_1_j", "facet_2_i", "facet_2_j",
        "facet_1_filtration", "facet_2_filtration",
        "length_A", "beta_rad", "beta_deg",
        "length_times_beta_A_rad",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for m in data["measurements"]:
            writer.writerow({
                "triangle_i": m.triangle[0],
                "triangle_j": m.triangle[1],
                "triangle_k": m.triangle[2],
                "singleton": m.singleton,
                "singleton_group": m.singleton_group,
                "facet_1_i": m.facet_1[0],
                "facet_1_j": m.facet_1[1],
                "facet_2_i": m.facet_2[0],
                "facet_2_j": m.facet_2[1],
                "facet_1_filtration": m.facet_1_filtration,
                "facet_2_filtration": m.facet_2_filtration,
                "length_A": m.length,
                "beta_rad": m.beta_rad,
                "beta_deg": m.beta_deg,
                "length_times_beta_A_rad": m.contribution,
            })
    return path


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("input", type=Path, help=".pdb or .npz")
    p.add_argument("--group-a-chain", action="append", default=[])
    p.add_argument("--group-b-chain", action="append", default=[])
    p.add_argument("--strict-heavy-fallbacks", action="store_true")
    p.add_argument("--include-hydrogens", action="store_true")
    p.add_argument("--include-hetatm", action="store_true")
    p.add_argument("--probe-radius", type=float, default=1.40)
    p.add_argument("--alpha-tol", type=float, default=ALPHA_TOL)
    p.add_argument("--alpha-edge-csv", type=Path)
    return p.parse_args()


def main():
    args = parse_args()

    if args.input.suffix.lower() == ".pdb":
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

    elif args.input.suffix.lower() == ".npz":
        points, radii, group_a, group_b = base.load_npz(args.input)

    else:
        raise SystemExit("Input must be .pdb or .npz")

    print()
    print("#" * 72)
    print("EXISTING BASELINE IMPLEMENTATION")
    print("#" * 72)
    baseline = base.analyze_bicolor_power_interface(
        points, radii, group_a, group_b, verbose=True
    )

    unrestricted = _measure(
        points, radii, group_a, group_b,
        alpha_zero=False,
        alpha_tol=args.alpha_tol,
    )
    _summary("CAZALS SIGN — UNRESTRICTED REGULAR INTERFACE", unrestricted)

    alpha0 = _measure(
        points, radii, group_a, group_b,
        alpha_zero=True,
        alpha_tol=args.alpha_tol,
    )
    _summary("CAZALS SIGN + CONDITION ALPHA (alpha = 0)", alpha0)

    print()
    print("=" * 72)
    print("COMPARISON")
    print("=" * 72)
    print(f"Original unrestricted baseline:    {baseline.mean_signed_angle_deg:.10f} deg")

    def signed_deg(data):
        ms = data["measurements"]
        L = sum(m.length for m in ms)
        C = sum(m.contribution for m in ms)
        return float(np.degrees(C / L)) if L else float("nan")

    print(f"Cazals sign, unrestricted:         {signed_deg(unrestricted):.10f} deg")
    print(f"Cazals sign + alpha=0:             {signed_deg(alpha0):.10f} deg")
    print("Published 2KAI benchmark:          approximately 17 deg")
    print()
    print("Condition beta (M=5): NOT IMPLEMENTED.")
    print("Do not interpret any remaining difference as an error until that")
    print("published large-facet rejection has been reproduced exactly.")
    print("=" * 72)

    if args.alpha_edge_csv:
        print(f"alpha=0 edge CSV: {export_csv(alpha0, args.alpha_edge_csv)}")


if __name__ == "__main__":
    main()
