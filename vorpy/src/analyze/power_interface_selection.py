"""Selection views of solved Power geometry; no curvature/geometry formulas.

The frozen Cazals result is never mutated. Full Power removes both alpha and M
filters, retaining bounded bicolor facets. Existing measured finite dual edges
are included precisely when both incident bicolor facets belong to the view.
Open-interface boundary edges are not assigned an invented turning curvature.
"""
from dataclasses import replace
from collections import defaultdict
from itertools import combinations

from vorpy.src.analyze.cazals_validation import power_facet_polygon_area


def power_selection_views(result):
    # Index incidence once. The existing polygon helper only needs tetrahedra
    # containing this facet's generator pair. Supplying that exact local subset
    # avoids repeatedly scanning the entire solvated complex, without changing
    # vertex coordinates, triangle incidence at the facet, or area arithmetic.
    pairs = {f.generator_ids for f in result.facets}
    incident = defaultdict(set)
    for tet in result.full_simplices[3]:
        for pair in combinations(tet, 2):
            if pair in pairs:
                incident[pair].add(tet)
    polygons = {}
    for facet in result.facets:
        local = replace(result, full_simplices={**result.full_simplices, 3: incident[facet.generator_ids]})
        polygons[facet.generator_ids] = power_facet_polygon_area(local, facet)
    bounded = {ids for ids, p in polygons.items() if p["status"] == "finite_bounded"}
    alpha = {f.generator_ids for f in result.facets if f.alpha_zero_selected} & bounded
    selections = {"full_power": bounded, "alpha0_power_before_M": alpha}
    views = {"cazals_alpha0_power": result}
    for method, selected in selections.items():
        facets = [replace(f, final_selected=f.generator_ids in selected,
                          status="selected" if f.generator_ids in selected else "excluded",
                          reason="" if f.generator_ids in selected else
                          polygons[f.generator_ids].get("reason") or "facet fails alpha=0 selection")
                  for f in result.facets]
        edges = []
        for edge in result.edges:
            inside = all(pair in selected for pair in edge.bicolor_facets)
            included = inside and edge.geometry_status == "finite"
            edges.append(replace(edge, status="included" if included else "excluded",
                                 reason="" if included else
                                 "one or both facets not selected" if not inside else edge.geometry_status))
        views[method] = replace(result, facets=facets, edges=edges)
    return views, polygons
