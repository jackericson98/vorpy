r"""
cazals_alpha0_audit.py

Audit the weighted GUDHI alpha complex used for the Cazals-style 2KAI
protein-protein interface calculation.

This is intentionally diagnostic. It does NOT claim that a GUDHI simplex-tree
filtration value by itself supplies clipped/restricted dual Voronoi geometry.

Place beside:
    cazals_power_interface.py
    cazals_pdb.py

Typical 2KAI run
----------------
python -m vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation.cazals_alpha0_audit ^
    vorpy\data\2kai.pdb ^
    --group-a-chain A --group-a-chain B ^
    --group-b-chain I ^
    --strict-heavy-fallbacks ^
    --out-dir 2kai_cazals_alpha0_audit

Outputs
-------
simplex_filtration.csv
bicolor_facets.csv
interface_edge_candidates.csv

Terminology
-----------
Regular edge (i,j)       <-> full power/Voronoi facet F_ij
Regular triangle (i,j,k) <-> full power/Voronoi edge E_ijk
Regular tetrahedron      <-> power/Voronoi vertex

At alpha=0 this script classifies simplices by the filtration values returned
by GUDHI's weighted AlphaComplex.  For a candidate interface edge it reports
the length of the FULL unrestricted dual power edge between its two incident
power vertices.  That full length is an audit quantity, not automatically the
length of a clipped/restricted alpha-shape interface edge.

Condition beta (M=5) is intentionally not implemented here.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np

from vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation import (
    cazals_power_interface as base,
)
from vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation.cazals_pdb import (
    prepare_cazals_pdb,
)


DEFAULT_ALPHA = 0.0
DEFAULT_TOL = 1.0e-12


def _chain_set(values):
    result = set()
    for value in values or []:
        for piece in value.split(","):
            piece = piece.strip()
            if piece:
                result.add(piece)
    return result


def _simplex_type(simplex):
    return {
        1: "vertex",
        2: "edge",
        3: "triangle",
        4: "tetrahedron",
    }.get(len(simplex), f"{len(simplex)-1}simplex")


def _filtration_map(st):
    out = {}
    for simplex, filtration in st.get_filtration():
        out[base.canonical_simplex(simplex)] = float(filtration)
    return out


def _in_alpha(filtration_value, alpha, tol):
    return filtration_value <= alpha + tol


def _group_label(i, group_a, group_b):
    if i in group_a:
        return "A"
    if i in group_b:
        return "B"
    return "?"


def _composition(simplex, group_a, group_b):
    return "".join(sorted(_group_label(i, group_a, group_b) for i in simplex))


def _triangle_singleton_and_beta(tri, points, group_a, group_b):
    """
    Diagnostic Cazals-style triangle angle:
      ABB -> + angle at A
      AAB -> - angle at B

    The sign convention is exposed explicitly so group reversal can be audited.
    """
    a = [i for i in tri if i in group_a]
    b = [i for i in tri if i in group_b]

    if len(a) == 1 and len(b) == 2:
        singleton = a[0]
        sign = +1.0
        singleton_group = "A"
    elif len(a) == 2 and len(b) == 1:
        singleton = b[0]
        sign = -1.0
        singleton_group = "B"
    else:
        raise ValueError("triangle is not 1:2 bicolor")

    others = [i for i in tri if i != singleton]
    u = points[others[0]] - points[singleton]
    v = points[others[1]] - points[singleton]
    nu = float(np.linalg.norm(u))
    nv = float(np.linalg.norm(v))
    if nu <= 1e-14 or nv <= 1e-14:
        raise ValueError("degenerate triangle")

    c = float(np.clip(np.dot(u, v) / (nu * nv), -1.0, 1.0))
    magnitude = float(np.arccos(c))
    return singleton, singleton_group, sign * magnitude


def _calculate_power_vertices(tetrahedra, points, weights, residual_tol):
    pvs = {}
    unresolved = {}
    for tet in tetrahedra:
        try:
            pv = base.calculate_power_vertex(tet, points, weights)
        except np.linalg.LinAlgError as exc:
            unresolved[tet] = f"solve:{exc}"
            continue
        if pv.residual > residual_tol:
            unresolved[tet] = f"residual:{pv.residual:.16g}"
            continue
        pvs[tet] = pv
    return pvs, unresolved


def _triangle_to_tets(tetrahedra):
    result = defaultdict(list)
    for tet in tetrahedra:
        for tri in combinations(tet, 3):
            result[base.canonical_simplex(tri)].append(tet)
    return dict(result)


def _edge_to_triangles(triangles):
    result = defaultdict(list)
    for tri in triangles:
        for edge in combinations(tri, 2):
            result[base.canonical_simplex(edge)].append(tri)
    return dict(result)


@dataclass
class Candidate:
    triangle: tuple[int, int, int]
    composition: str
    singleton: int
    singleton_group: str
    beta_rad: float
    beta_deg: float
    bicolor_facet_1: tuple[int, int]
    bicolor_facet_2: tuple[int, int]
    facet_1_filtration: float
    facet_2_filtration: float
    triangle_filtration: float
    facet_1_alpha: bool
    facet_2_alpha: bool
    triangle_alpha: bool
    incident_tetrahedra: tuple
    finite_full_dual_edge: bool
    full_dual_length: float | None
    full_length_times_beta: float | None
    status: str


def analyze(points, radii, group_a, group_b, alpha=0.0, tol=1e-12,
            residual_tol=1e-8):
    points, radii, group_a, group_b = base.validate_inputs(
        points, radii, group_a, group_b
    )
    weights = radii ** 2

    st, vertices, edges, triangles, tetrahedra = base.build_regular_complex(
        points, radii
    )
    filt = _filtration_map(st)

    pvs, unresolved_pvs = _calculate_power_vertices(
        tetrahedra, points, weights, residual_tol
    )
    tri_to_tets = _triangle_to_tets(tetrahedra)
    edge_to_tris = _edge_to_triangles(triangles)

    alpha_vertices = {v for v in vertices if _in_alpha(filt[(v,)], alpha, tol)}
    alpha_edges = {e for e in edges if _in_alpha(filt[e], alpha, tol)}
    alpha_triangles = {t for t in triangles if _in_alpha(filt[t], alpha, tol)}
    alpha_tetrahedra = {t for t in tetrahedra if _in_alpha(filt[t], alpha, tol)}

    bicolor_edges = {e for e in edges if base.is_bicolor_edge(e, group_a, group_b)}
    alpha_bicolor_edges = bicolor_edges & alpha_edges

    candidates = []
    for tri in sorted(triangles):
        bi = base.bicolor_edges_of_triangle(tri, group_a, group_b)
        if len(bi) != 2:
            continue

        try:
            singleton, singleton_group, beta = _triangle_singleton_and_beta(
                tri, points, group_a, group_b
            )
        except ValueError:
            continue

        facet1, facet2 = sorted(bi)
        f1 = filt[facet1]
        f2 = filt[facet2]
        tf = filt[tri]
        f1a = _in_alpha(f1, alpha, tol)
        f2a = _in_alpha(f2, alpha, tol)
        ta = _in_alpha(tf, alpha, tol)

        incident = tuple(tri_to_tets.get(tri, ()))
        finite = False
        full_len = None
        contribution = None
        status_parts = []

        if len(incident) == 1:
            status_parts.append("unbounded_full_power_edge")
        elif len(incident) != 2:
            status_parts.append(f"nonmanifold_incidence_{len(incident)}")
        elif incident[0] not in pvs or incident[1] not in pvs:
            status_parts.append("unresolved_power_vertex")
        else:
            finite = True
            full_len = float(
                np.linalg.norm(
                    pvs[incident[1]].position - pvs[incident[0]].position
                )
            )
            contribution = full_len * beta
            status_parts.append("finite_full_power_edge")

        if f1a and f2a:
            status_parts.append("both_bicolor_facets_alpha")
        elif f1a or f2a:
            status_parts.append("one_bicolor_facet_alpha")
        else:
            status_parts.append("no_bicolor_facets_alpha")

        if ta:
            status_parts.append("triangle_alpha")
        else:
            status_parts.append("triangle_not_alpha")

        candidates.append(
            Candidate(
                triangle=tri,
                composition=_composition(tri, group_a, group_b),
                singleton=singleton,
                singleton_group=singleton_group,
                beta_rad=beta,
                beta_deg=float(np.degrees(beta)),
                bicolor_facet_1=facet1,
                bicolor_facet_2=facet2,
                facet_1_filtration=f1,
                facet_2_filtration=f2,
                triangle_filtration=tf,
                facet_1_alpha=f1a,
                facet_2_alpha=f2a,
                triangle_alpha=ta,
                incident_tetrahedra=incident,
                finite_full_dual_edge=finite,
                full_dual_length=full_len,
                full_length_times_beta=contribution,
                status=";".join(status_parts),
            )
        )

    return {
        "st": st,
        "filtration": filt,
        "vertices": vertices,
        "edges": edges,
        "triangles": triangles,
        "tetrahedra": tetrahedra,
        "alpha_vertices": alpha_vertices,
        "alpha_edges": alpha_edges,
        "alpha_triangles": alpha_triangles,
        "alpha_tetrahedra": alpha_tetrahedra,
        "bicolor_edges": bicolor_edges,
        "alpha_bicolor_edges": alpha_bicolor_edges,
        "candidates": candidates,
        "power_vertices": pvs,
        "unresolved_power_vertices": unresolved_pvs,
        "edge_to_triangles": edge_to_tris,
        "group_a": group_a,
        "group_b": group_b,
        "points": points,
        "radii": radii,
        "alpha": alpha,
        "tol": tol,
    }


def _fmt_simplex(s):
    return " ".join(map(str, s))


def export_simplex_filtration(data, out_dir):
    path = out_dir / "simplex_filtration.csv"
    filt = data["filtration"]
    with path.open("w", newline="", encoding="utf-8") as h:
        fields = [
            "dimension", "type", "simplex", "composition",
            "filtration", "in_alpha",
        ]
        w = csv.DictWriter(h, fieldnames=fields)
        w.writeheader()
        for simplex, value in sorted(
            filt.items(), key=lambda kv: (len(kv[0]), kv[0])
        ):
            w.writerow({
                "dimension": len(simplex) - 1,
                "type": _simplex_type(simplex),
                "simplex": _fmt_simplex(simplex),
                "composition": _composition(
                    simplex, data["group_a"], data["group_b"]
                ),
                "filtration": f"{value:.17g}",
                "in_alpha": int(
                    _in_alpha(value, data["alpha"], data["tol"])
                ),
            })
    return path


def export_bicolor_facets(data, out_dir):
    path = out_dir / "bicolor_facets.csv"
    filt = data["filtration"]
    with path.open("w", newline="", encoding="utf-8") as h:
        fields = [
            "i", "j", "filtration", "in_alpha",
            "incident_regular_triangles",
            "incident_alpha_triangles",
        ]
        w = csv.DictWriter(h, fieldnames=fields)
        w.writeheader()
        for edge in sorted(data["bicolor_edges"]):
            tris = data["edge_to_triangles"].get(edge, [])
            alpha_tris = [t for t in tris if t in data["alpha_triangles"]]
            w.writerow({
                "i": edge[0],
                "j": edge[1],
                "filtration": f"{filt[edge]:.17g}",
                "in_alpha": int(edge in data["alpha_bicolor_edges"]),
                "incident_regular_triangles": "|".join(
                    _fmt_simplex(t) for t in tris
                ),
                "incident_alpha_triangles": "|".join(
                    _fmt_simplex(t) for t in alpha_tris
                ),
            })
    return path


def export_candidates(data, out_dir):
    path = out_dir / "interface_edge_candidates.csv"
    with path.open("w", newline="", encoding="utf-8") as h:
        fields = [
            "triangle", "composition", "singleton", "singleton_group",
            "beta_rad", "beta_deg",
            "facet_1", "facet_2",
            "facet_1_filtration", "facet_2_filtration",
            "triangle_filtration",
            "facet_1_alpha", "facet_2_alpha", "triangle_alpha",
            "incident_tetrahedra",
            "finite_full_dual_edge",
            "full_dual_length_A",
            "full_length_times_beta_A_rad",
            "status",
        ]
        w = csv.DictWriter(h, fieldnames=fields)
        w.writeheader()
        for c in data["candidates"]:
            w.writerow({
                "triangle": _fmt_simplex(c.triangle),
                "composition": c.composition,
                "singleton": c.singleton,
                "singleton_group": c.singleton_group,
                "beta_rad": f"{c.beta_rad:.17g}",
                "beta_deg": f"{c.beta_deg:.17g}",
                "facet_1": _fmt_simplex(c.bicolor_facet_1),
                "facet_2": _fmt_simplex(c.bicolor_facet_2),
                "facet_1_filtration": f"{c.facet_1_filtration:.17g}",
                "facet_2_filtration": f"{c.facet_2_filtration:.17g}",
                "triangle_filtration": f"{c.triangle_filtration:.17g}",
                "facet_1_alpha": int(c.facet_1_alpha),
                "facet_2_alpha": int(c.facet_2_alpha),
                "triangle_alpha": int(c.triangle_alpha),
                "incident_tetrahedra": "|".join(
                    _fmt_simplex(t) for t in c.incident_tetrahedra
                ),
                "finite_full_dual_edge": int(c.finite_full_dual_edge),
                "full_dual_length_A": (
                    "" if c.full_dual_length is None
                    else f"{c.full_dual_length:.17g}"
                ),
                "full_length_times_beta_A_rad": (
                    "" if c.full_length_times_beta is None
                    else f"{c.full_length_times_beta:.17g}"
                ),
                "status": c.status,
            })
    return path


def _euler(v, e, f, t):
    return len(v) - len(e) + len(f) - len(t)


def print_summary(data):
    alpha = data["alpha"]
    candidates = data["candidates"]

    both_facets = [
        c for c in candidates if c.facet_1_alpha and c.facet_2_alpha
    ]
    triangle_alpha = [c for c in candidates if c.triangle_alpha]
    finite_both = [
        c for c in both_facets if c.finite_full_dual_edge
    ]

    def diagnostic_stats(items):
        usable = [
            c for c in items
            if c.full_dual_length is not None
            and c.full_length_times_beta is not None
        ]
        L = float(sum(c.full_dual_length for c in usable))
        C = float(sum(c.full_length_times_beta for c in usable))
        U = float(sum(abs(c.full_length_times_beta) for c in usable))
        pos = float(sum(
            c.full_length_times_beta for c in usable
            if c.beta_rad > 0
        ))
        neg = float(sum(
            c.full_length_times_beta for c in usable
            if c.beta_rad < 0
        ))
        return usable, L, C, U, pos, neg

    usable, L, C, U, pos, neg = diagnostic_stats(finite_both)

    print()
    print("=" * 78)
    print("GUDHI WEIGHTED ALPHA=0 / CAZALS INTERFACE AUDIT")
    print("=" * 78)
    print(f"Alpha:                              {alpha:.12g}")
    print(f"Input atoms:                        {len(data['points']):,}")
    print(f"Group A / B:                        {len(data['group_a']):,} / {len(data['group_b']):,}")
    print("-" * 78)
    print("FULL REGULAR COMPLEX")
    print(f"  vertices:                         {len(data['vertices']):,}")
    print(f"  edges:                            {len(data['edges']):,}")
    print(f"  triangles:                        {len(data['triangles']):,}")
    print(f"  tetrahedra:                       {len(data['tetrahedra']):,}")
    print(f"  Euler characteristic:             {_euler(data['vertices'], data['edges'], data['triangles'], data['tetrahedra']):,}")
    print(f"  bicolor regular edges/facets:     {len(data['bicolor_edges']):,}")
    print("-" * 78)
    print("GUDHI ALPHA COMPLEX")
    print(f"  vertices:                         {len(data['alpha_vertices']):,}")
    print(f"  edges:                            {len(data['alpha_edges']):,}")
    print(f"  triangles:                        {len(data['alpha_triangles']):,}")
    print(f"  tetrahedra:                       {len(data['alpha_tetrahedra']):,}")
    print(f"  Euler characteristic:             {_euler(data['alpha_vertices'], data['alpha_edges'], data['alpha_triangles'], data['alpha_tetrahedra']):,}")
    print(f"  bicolor edges/facets retained:    {len(data['alpha_bicolor_edges']):,}")
    print("-" * 78)
    print("BICOLOR TRIANGLE / DUAL-EDGE AUDIT")
    print(f"  all 1:2 bicolor triangles:        {len(candidates):,}")
    print(f"  both bicolor facets in alpha:     {len(both_facets):,}")
    print(f"  triangle itself in alpha:         {len(triangle_alpha):,}")
    print(f"  finite full dual edges (both):    {len(finite_both):,}")
    print(f"  unresolved power vertices:        {len(data['unresolved_power_vertices']):,}")
    print("-" * 78)
    print("UNCLIPPED FULL-DUAL LENGTH DIAGNOSTIC")
    print("  WARNING: these lengths are full power-edge lengths.")
    print("  They are NOT asserted to be restricted alpha-interface lengths.")
    print(f"  full dual length sum:              {L:.10f} A")
    print(f"  positive sum(l*beta):              {pos:.10f} A rad")
    print(f"  negative sum(l*beta):              {neg:.10f} A rad")
    print(f"  signed sum(l*beta):                {C:.10f} A rad")
    print(f"  unsigned sum(l*|beta|):            {U:.10f} A rad")
    if L:
        print(f"  signed normalized diagnostic:      {np.degrees(C/L):.10f} deg")
        print(f"  unsigned normalized diagnostic:    {np.degrees(U/L):.10f} deg")
    print("-" * 78)
    print("NOT YET CLAIMED / IMPLEMENTED")
    print("  * clipped/restricted dual Voronoi edge lengths")
    print("  * exact Cazals condition beta (M=5)")
    print("  * final published s_H reproduction")
    print("=" * 78)
    print()


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("input", type=Path, help=".pdb or .npz")
    p.add_argument("--group-a-chain", action="append", default=[])
    p.add_argument("--group-b-chain", action="append", default=[])
    p.add_argument("--strict-heavy-fallbacks", action="store_true")
    p.add_argument("--include-hydrogens", action="store_true")
    p.add_argument("--include-hetatm", action="store_true")
    p.add_argument("--probe-radius", type=float, default=1.40)
    p.add_argument("--alpha", type=float, default=DEFAULT_ALPHA)
    p.add_argument("--tol", type=float, default=DEFAULT_TOL)
    p.add_argument("--out-dir", type=Path, default=Path("cazals_alpha0_audit"))
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

    data = analyze(
        points,
        radii,
        group_a,
        group_b,
        alpha=args.alpha,
        tol=args.tol,
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    p1 = export_simplex_filtration(data, args.out_dir)
    p2 = export_bicolor_facets(data, args.out_dir)
    p3 = export_candidates(data, args.out_dir)

    print_summary(data)
    print(f"Simplex filtration CSV:        {p1}")
    print(f"Bicolor facet CSV:             {p2}")
    print(f"Interface-edge candidate CSV:  {p3}")


if __name__ == "__main__":
    main()
