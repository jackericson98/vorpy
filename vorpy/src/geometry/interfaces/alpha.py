"""Map scheme-native alpha-selected dual edges to solved Voronoi surfaces."""

from __future__ import annotations

import numpy as np

from .builders import build_full_interface
from .curvature import _summary
from .model import InterfaceEdge, InterfaceRepresentation, InterfaceVertex


def build_alpha_interface(network, dual, filtration, alpha, group_a, group_b, *,
                          calculate_curvature=True, aw_options=None):
    """Build a physical interface selected by a scheme-specific filtration.

    Alpha simplices decide which generator pairs survive. Every area, boundary,
    and curvature value is then obtained from the same solved Voronoi-side
    network. The filtration object is queried as-is; Cazals M selection is
    never applied here.
    """
    scheme = str((getattr(network, "settings", None) or {}).get("net_type", "aw"))
    if filtration.network is not network or dual.network is not network:
        raise ValueError("Network, dual, and filtration must reference the same solve.")
    if filtration.scheme != scheme or dual.scheme != scheme:
        raise ValueError("Network, dual, and filtration schemes must match.")
    group_a, group_b = frozenset(map(int, group_a)), frozenset(map(int, group_b))
    base = build_full_interface(network, dual, group_a, group_b)
    selected_edges = filtration.interface_at(alpha, group_a, group_b)
    selected_pairs = {tuple(sorted(row.generator_tuple)) for row in selected_edges}
    pair_surfaces = {}
    for surface in base.surfaces:
        pair_surfaces.setdefault(tuple(sorted(surface.generator_ids)), []).append(surface)

    mapping_rows = []
    selected_surface_ids = set()
    for pair in sorted(selected_pairs):
        surfaces = pair_surfaces.get(pair, [])
        supported = [surface for surface in surfaces if surface.supported]
        unique = len(supported) == 1
        surface = supported[0] if unique else (surfaces[0] if len(surfaces) == 1 else None)
        mapped = surface is not None
        eligible = bool(
            surface and surface.supported and surface.bounded and surface.complete
            and surface.area_A2 is not None
        )
        reason = "" if eligible else (
            "no physical bicolor surface exists" if not surfaces else
            "multiple supported physical surfaces map to one dual edge" if len(supported) > 1 else
            "mapped surface is incomplete, unsupported, or geometrically excluded"
        )
        atoms = [_stable_atom(network, atom_id, group_a, group_b, filtration) for atom_id in pair]
        birth = filtration.simplex_birth(pair)
        mapping_rows.append({
            "scheme": scheme,
            "alpha_native": float(alpha),
            "alpha_units": filtration.metadata.get("native_alpha_units", "unknown"),
            "generator_i": pair[0], "generator_j": pair[1],
            "stable_atom_i": atoms[0]["stable_atom_id"],
            "stable_atom_j": atoms[1]["stable_atom_id"],
            "group_i": atoms[0]["group"], "group_j": atoms[1]["group"],
            "dual_edge_id": birth.simplex_id if birth else f"1:{pair[0]},{pair[1]}",
            "surface_id": surface.feature_id if surface else None,
            "mapped": mapped,
            "surface_complete": bool(surface and surface.complete),
            "surface_supported": bool(surface and surface.supported),
            "surface_cells_complete": bool(surface and surface.generator_cells_complete),
            "curvature_eligible": bool(surface and surface.selected),
            "selected": eligible,
            "exclusion_reason": reason,
        })
        if eligible:
            selected_surface_ids.add(surface.feature_id)

    surfaces = []
    for surface in base.surfaces:
        selected = surface.feature_id in selected_surface_ids
        atom_i, atom_j = [
            _stable_atom(network, atom_id, group_a, group_b, filtration)
            for atom_id in surface.generator_ids
        ]
        surfaces.append(type(surface)(
            surface.feature_id, surface.generator_ids, surface.dual_simplex_id,
            surface.area_A2, surface.bounded, surface.complete,
            surface.generator_cells_complete, surface.supported, selected,
            "selected" if selected else "excluded_alpha" if surface.selected else surface.status,
            (("selected by alpha; AW curvature excluded because participating cells are incomplete"
              if not surface.generator_cells_complete else "")
             if selected else "dual edge is not present at the requested alpha"
             if surface.selected else surface.reason),
            surface.boundary_edge_ids, surface.boundary_vertex_ids,
            surface.smooth_H_A, atom_i["stable_atom_id"], atom_j["stable_atom_id"],
            atom_i["group"], atom_j["group"],
        ))

    alpha_edge_records = []
    incidence_rows = []
    for edge in base.edges:
        selected_incident = tuple(
            surface_id for surface_id in edge.incident_candidate_surface_ids
            if surface_id in selected_surface_ids
        )
        count = len(selected_incident)
        triangle = tuple(sorted(edge.generator_ids))
        triangle_birth = filtration.simplex_birth(triangle)
        present = bool(triangle_birth and triangle_birth.supported
                       and triangle_birth.filtration_birth is not None
                       and triangle_birth.filtration_birth <= float(alpha) + filtration.tolerance)
        tri_status = (
            "present_at_alpha" if present else
            "unresolved" if triangle_birth is None or not triangle_birth.supported
            or triangle_birth.filtration_birth is None else "resolved_absent_at_alpha"
        )
        if not edge.bounded or not edge.complete or not edge.supported:
            classification = "unresolved"
        elif count > 2:
            classification = "nonmanifold"
        elif count == 2:
            classification = "interior_selected"
        elif count == 1:
            classification = "boundary_selected"
        else:
            classification = "excluded"
        reason = {
            "interior_selected": "physical edge is shared by two alpha-selected surfaces",
            "boundary_selected": "physical edge bounds one alpha-selected surface; no boundary curvature term",
            "excluded": "physical edge has no alpha-selected surface incident",
            "nonmanifold": "more than two alpha-selected surfaces meet this physical edge",
        }.get(classification, edge.reason)
        updated = InterfaceEdge(
            edge.feature_id, edge.generator_ids, edge.dual_simplex_id,
            edge.incident_candidate_surface_ids, selected_incident,
            edge.bounded, edge.complete, edge.supported, classification,
            reason, edge.endpoint_vertex_ids,
            classification=classification,
            dual_triangle_present_at_alpha=present,
            dual_triangle_birth=(triangle_birth.filtration_birth if triangle_birth else None),
            dual_triangle_status=tri_status,
            curvature_status=("boundary_no_term" if classification == "boundary_selected"
                              else "not_computed"),
        )
        alpha_edge_records.append(updated)
        if count:
            incidence_rows.append({
                "scheme": scheme, "alpha_native": float(alpha),
                "alpha_units": filtration.metadata.get("native_alpha_units", "unknown"),
                "physical_edge_id": edge.feature_id,
                "generator_tuple": triangle,
                "selected_surface_ids": selected_incident,
                "dual_triangle_tuple": triangle,
                "dual_triangle_present_at_alpha": present,
                "dual_triangle_birth": triangle_birth.filtration_birth if triangle_birth else None,
                "dual_triangle_status": tri_status,
                "classification": classification,
            })

    selected_vertex_records = _selected_vertices(
        network, dual, base.surfaces, selected_surface_ids
    )

    result = InterfaceRepresentation(
        scheme, "alpha", float(alpha), group_a, group_b, surfaces,
        alpha_edge_records, selected_vertex_records,
        metadata={
            "alpha_native": float(alpha),
            "alpha_units": filtration.metadata.get("native_alpha_units", "unknown"),
            "alpha_convention": filtration.metadata.get("filtration_kind", "scheme-native"),
            "physical_geometry_source": "same solved network Voronoi surfaces, edges, and vertices",
            "selection_source": "same-scheme alpha filtration, bicolor dual 1-simplices",
            "cazals_alpha0_m5_applied": False,
            "cazals_selection_is_separate": True,
            "selected_pair_count": len(selected_pairs),
            "mapped_pair_count": sum(row["mapped"] for row in mapping_rows),
            "selected_surface_count": len(selected_surface_ids),
            "boundary_curvature_policy": "no boundary curvature term; edge turning is integrated only when two selected patches meet",
            "topology_note": "incidence data may be incomplete when higher-dimensional alpha births are unresolved",
        },
        selection_mappings=mapping_rows,
        alpha_incidence_audit=incidence_rows,
    )
    triangle_counts = {}
    for row in incidence_rows:
        status = row["dual_triangle_status"]
        triangle_counts[status] = triangle_counts.get(status, 0) + 1
    result.metadata["physical_interface_edge_count"] = len(incidence_rows)
    result.metadata["dual_triangle_status_counts"] = triangle_counts
    if calculate_curvature:
        result.curvature = _curvature_for_selected_interface(
            result, filtration, float(alpha), aw_options or {}
        )
    return result


def _curvature_for_selected_interface(interface, filtration, alpha, aw_options):
    scheme = interface.scheme
    selected_ids = {surface.feature_id for surface in interface.selected_surfaces}
    area = interface.area_A2
    if scheme == "aw":
        from vorpy.src.analyze.aw_interface_curvature import (
            analyze_aw_interface_curvature,
        )

        full = analyze_aw_interface_curvature(
            filtration.network, interface.group_a, interface.group_b,
            selection_mode="supported", **aw_options,
        )
        surface_curvature = {
            record.surface_id: record.integrated_mean_curvature
            for record in full.selected_surfaces
            if record.integrated_mean_curvature is not None
        }
        for index, surface in enumerate(interface.surfaces):
            if surface.feature_id in selected_ids:
                value = surface_curvature.get(surface.feature_id)
                interface.surfaces[index] = type(surface)(
                    surface.feature_id, surface.generator_ids, surface.dual_simplex_id,
                    surface.area_A2, surface.bounded, surface.complete,
                    surface.generator_cells_complete, surface.supported, surface.selected,
                    surface.status, surface.reason, surface.boundary_edge_ids,
                    surface.boundary_vertex_ids, value, surface.stable_atom_i,
                    surface.stable_atom_j, surface.group_i, surface.group_j,
                )
        missing_surface_curvature = selected_ids - set(surface_curvature)
        smooth = (None if missing_surface_curvature else
                  sum(surface_curvature[sid] for sid in selected_ids))
        edge_by_id = {edge.edge_id: edge for edge in full.selected_edges
                      if set(edge.surface_ids) <= selected_ids}
        signed = unsigned = positive = negative = 0.0
        for i, edge in enumerate(interface.edges):
            measured = edge_by_id.get(edge.feature_id)
            if edge.classification != "interior_selected":
                continue
            if measured is None:
                interface.edges[i] = _with_curvature_status(edge, "unresolved_curvature")
                continue
            signed += measured.signed_contribution or 0.0
            unsigned += measured.unsigned_contribution or 0.0
            positive += measured.positive_contribution or 0.0
            negative += measured.negative_contribution or 0.0
            interface.edges[i] = _with_edge_measure(edge, measured.length,
                measured.signed_contribution, measured.unsigned_contribution,
                None)
        summary = _summary(
            "aw", "alpha", area, smooth, signed, unsigned, positive, negative,
            "interface normal A-to-B; reversal flips signed smooth/edge/combined terms",
            "existing orientation-aware AW curvature analyzer; selected surface/edge records filtered by alpha pair mapping",
        )
        if missing_surface_curvature:
            from dataclasses import replace
            summary = replace(
                summary, status="partial",
                source=summary.source +
                f"; smooth term unresolved for {len(missing_surface_curvature)} selected surfaces with incomplete AW cells",
            )
        return summary

    for index, surface in enumerate(interface.surfaces):
        if surface.selected:
            interface.surfaces[index] = _surface_with_smooth(surface, 0.0)
    return _planar_curvature(interface, filtration, alpha, area)


def _planar_curvature(interface, filtration, alpha, area):
    """Reuse the validated Power edge measurement engine for planar facets."""
    from vorpy.src.analyze.power_interface_curvature import (
        analyze_power_regular_complex,
    )

    network = filtration.network
    ids = [int(value) for value in network.balls.index]
    position = {generator_id: index for index, generator_id in enumerate(ids)}
    points = np.asarray(network.balls.loc[ids, "loc"].tolist(), dtype=float)
    radii = (np.asarray(network.balls.loc[ids, "rad"], dtype=float)
             if interface.scheme == "pow" else np.ones(len(ids), dtype=float))
    full, births = {dimension: set() for dimension in range(4)}, {}
    for dimension, records in filtration.records.items():
        for key, record in records.items():
            compact = tuple(sorted(position[int(value)] for value in key))
            full[dimension].add(compact)
            if record.filtration_birth is not None:
                births[compact] = float(record.filtration_birth)
    group_a = {position[value] for value in interface.group_a}
    group_b = {position[value] for value in interface.group_b}
    # A large finite M disables only the Cazals size rejection. Alpha pair
    # selection still comes from the supplied scheme-specific filtration.
    measured = analyze_power_regular_complex(
        points, radii, group_a, group_b, full, births, alpha=alpha,
        condition_beta_m=1e300,
    )
    selected_pairs = {tuple(sorted(position[value] for value in pair))
                      for pair in (surface.generator_ids for surface in interface.selected_surfaces)}
    edges_by_triangle = {}
    for record in measured.edges:
        if set(record.bicolor_facets) <= selected_pairs:
            edges_by_triangle.setdefault(tuple(record.regular_triangle), record)
    signed = unsigned = positive = negative = 0.0
    for index, edge in enumerate(interface.edges):
        if edge.classification != "interior_selected":
            continue
        compact_tri = tuple(sorted(position[value] for value in edge.generator_ids))
        measurement = edges_by_triangle.get(compact_tri)
        if measurement is None or measurement.status != "included":
            interface.edges[index] = _with_curvature_status(edge, "unresolved_curvature")
            continue
        s = float(measurement.length_beta or 0.0)
        u = float(measurement.length_abs_beta or 0.0)
        signed += s
        unsigned += u
        positive += max(s, 0.0)
        negative += min(s, 0.0)
        beta = (s / measurement.length) if measurement.length else None
        interface.edges[index] = _with_edge_measure(edge, measurement.length, s, u, beta)
    orientation = "Cazals AAB/ABB sign: ABB positive and AAB negative; group labels fixed"
    if interface.scheme == "prm":
        orientation += "; uniform unit weights invoke the validated planar Power engine and preserve ordinary Voronoi geometry"
    return _summary(
        "power" if interface.scheme == "pow" else "prm", "alpha", area, 0.0,
        signed, unsigned, positive, negative, orientation,
        "existing Power planar edge geometry/curvature analyzer; M rejection disabled for generic alpha interface",
    )


def _with_edge_measure(edge, length, signed, unsigned, beta):
    return type(edge)(
        edge.feature_id, edge.generator_ids, edge.dual_simplex_id,
        edge.incident_candidate_surface_ids, edge.incident_selected_surface_ids,
        edge.bounded, edge.complete, edge.supported, edge.incidence_status,
        edge.reason, edge.endpoint_vertex_ids, edge.classification,
        float(length) if length is not None else None,
        float(beta) if beta is not None else None,
        float(signed) if signed is not None else None,
        float(unsigned) if unsigned is not None else None,
        edge.dual_triangle_present_at_alpha, edge.dual_triangle_birth,
        edge.dual_triangle_status, "integrated",
    )


def _surface_with_smooth(surface, value):
    return type(surface)(
        surface.feature_id, surface.generator_ids, surface.dual_simplex_id,
        surface.area_A2, surface.bounded, surface.complete,
        surface.generator_cells_complete, surface.supported, surface.selected,
        surface.status, surface.reason, surface.boundary_edge_ids,
        surface.boundary_vertex_ids, value, surface.stable_atom_i,
        surface.stable_atom_j, surface.group_i, surface.group_j,
    )


def _with_curvature_status(edge, status):
    return type(edge)(
        edge.feature_id, edge.generator_ids, edge.dual_simplex_id,
        edge.incident_candidate_surface_ids, edge.incident_selected_surface_ids,
        edge.bounded, edge.complete, edge.supported, edge.incidence_status,
        edge.reason, edge.endpoint_vertex_ids, edge.classification, edge.length_A,
        edge.signed_turn_or_beta, edge.signed_integral_A_rad,
        edge.unsigned_integral_A_rad, edge.dual_triangle_present_at_alpha,
        edge.dual_triangle_birth, edge.dual_triangle_status, status,
    )


def _stable_atom(network, generator_id, group_a, group_b, filtration):
    row = network.balls.loc[generator_id]
    try:
        from vorpy.src.geometry.aw_alpha.complex import _network_atom_metadata
        source = _network_atom_metadata(network, generator_id, row)
    except (AttributeError, KeyError, TypeError, ValueError):
        source = {}
    def first(names):
        for name in names:
            value = row.get(name)
            if value is not None and str(value).strip() not in {"", "nan", "None"}:
                return str(value).strip()
            value = source.get(name)
            if value is not None and str(value).strip() not in {"", "nan", "None"}:
                return str(value).strip()
        return ""
    chain = first(("chain_name", "chain"))
    resnum = first(("auth_seq_id", "res_seq", "residue_number"))
    icode = first(("pdb_ins_code", "insertion_code", "icode"))
    resname = first(("auth_comp_id", "res_name", "residue_name"))
    atom = first(("auth_atom_id", "name", "atom_name"))
    override = filtration.metadata.get("stable_atom_metadata", {}).get(generator_id)
    if override:
        if isinstance(override, dict):
            chain = str(override.get("chain", chain))
            resnum = str(override.get("residue_number", resnum))
            icode = str(override.get("insertion_code", icode))
            resname = str(override.get("residue_name", resname))
            atom = str(override.get("atom_name", atom))
            stable = str(override.get("stable_atom_id", ""))
        else:
            chain, resnum, icode, resname, atom = map(str, override)
            stable = "|".join((chain, resnum + icode, resname, atom))
    else:
        stable = (str(source.get("stable_id", "")) or first(("stable_id",))
                  or f"{chain}|{resnum}{icode}|{resname}|{atom}")
    group = "A" if generator_id in group_a else "B" if generator_id in group_b else ""
    return {"stable_atom_id": stable, "chain": chain, "residue_number": resnum,
            "insertion_code": icode, "residue_name": resname, "atom_name": atom,
            "group": group}


def _selected_vertices(network, dual, candidate_surfaces, selected_surface_ids):
    candidate_by_vertex = {}
    for surface in candidate_surfaces:
        for vertex_id in surface.boundary_vertex_ids:
            candidate_by_vertex.setdefault(vertex_id, set()).add(surface.feature_id)
    output = []
    for vertex_id, candidate_ids in sorted(candidate_by_vertex.items()):
        selected = tuple(sorted(candidate_ids & selected_surface_ids))
        if not selected or vertex_id not in network.verts.index:
            continue
        row = network.verts.loc[vertex_id]
        try:
            generator_ids = tuple(sorted(int(value) for value in row.get("balls", ())))
        except (TypeError, ValueError):
            generator_ids = ()
        simplex = dual.simplex(3, generator_ids) if len(generator_ids) == 4 else None
        feature = next((item for item in (simplex.primal_features if simplex else ())
                        if item.kind == "vert" and item.feature_id == int(vertex_id)), None)
        bounded = bool(feature and feature.bounded)
        supported = bool(simplex and simplex.supported and feature)
        try:
            xyz = tuple(float(value) for value in row.get("loc"))
            if len(xyz) != 3 or not np.all(np.isfinite(xyz)):
                xyz = None
        except (TypeError, ValueError):
            xyz = None
        output.append(InterfaceVertex(
            int(vertex_id), generator_ids, tuple(sorted(candidate_ids)), selected,
            bounded, supported,
            "supported" if bounded and supported else "unresolved_or_incomplete",
            xyz,
        ))
    return output
