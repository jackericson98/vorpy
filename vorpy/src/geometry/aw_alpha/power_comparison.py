"""Diagnostic Power alpha=0 versus experimental AW alpha=0 pair overlap."""

from __future__ import annotations

import numpy as np

from vorpy.src.analyze.cazals_validation import power_facet_polygon_area
from vorpy.src.analyze.power_interface_curvature import (
    analyze_power_interface_curvature,
)

from .complex import AWAlphaFiltration, _network_atom_metadata


def compare_alpha0_pairs(
    network, group_a, group_b, *, name=None, tolerance=1e-6, filtration=None
):
    """Compare alpha=0 bicolor generator pairs using the same centers/radii.

    The Power branch is the existing weighted GUDHI alpha implementation,
    with the AW network's exact current per-generator radii supplied as its
    power radii. M=5 is not used to determine the Power alpha-zero pair set.
    """
    if not hasattr(network, "balls") or not hasattr(network, "surfs"):
        raise ValueError("Pair comparison requires an existing solved AW network.")
    ids = [int(value) for value in network.balls.index]
    locations = np.asarray([network.balls.loc[i, "loc"] for i in ids], dtype=float)
    radii = np.asarray([network.balls.loc[i, "rad"] for i in ids], dtype=float)
    index_pos = {generator_id: position for position, generator_id in enumerate(ids)}
    network_group_a = set(map(int, group_a))
    network_group_b = set(map(int, group_b))
    group_a = {index_pos[i] for i in network_group_a}
    group_b = {index_pos[i] for i in network_group_b}
    power = analyze_power_interface_curvature(
        locations, radii, group_a, group_b, alpha=0.0, condition_beta_m=5.0
    )
    power_pairs = {
        tuple(sorted(facet.generator_ids)) for facet in power.alpha_zero_facets
    }
    power_area = {}
    for facet in power.alpha_zero_facets:
        pair = tuple(sorted(facet.generator_ids))
        if not (
            (pair[0] in group_a and pair[1] in group_b)
            or (pair[1] in group_a and pair[0] in group_b)
        ):
            continue
        area = power_facet_polygon_area(power, facet)
        power_area[pair] = (
            area.get("area_A2") if area.get("status") == "finite_bounded" else None
        )

    # Restore network-local stable generator IDs after the existing analyzer's
    # compact positional remapping.
    power_pairs = {
        tuple(sorted((ids[pair[0]], ids[pair[1]])))
        for pair in power_pairs
        if pair[0] in group_a
        and pair[1] in group_b
        or pair[1] in group_a
        and pair[0] in group_b
    }
    power_area = {
        tuple(sorted((ids[pair[0]], ids[pair[1]]))): area
        for pair, area in power_area.items()
    }
    aw = (
        filtration or AWAlphaFiltration(network, name=name, tolerance=tolerance)
    ).calculate_births(max_dimension=1)
    aw_rows, aw_summary = aw.bicolor_surface_records(
        0.0, network_group_a, network_group_b, max_dimension=1
    )
    aw_pairs = {tuple(sorted(row["generator_tuple"])) for row in aw_rows}
    # A pair can map to disconnected primal surface records; pair counts and
    # overlap are unique generator tuples, while mapped surface counts remain
    # separate in the AW summary.
    shared = power_pairs & aw_pairs
    union = power_pairs | aw_pairs
    pair_rows = []
    aw_area = {}
    for row in aw_rows:
        pair = tuple(sorted(row["generator_tuple"]))
        aw_area[pair] = aw_area.get(pair, 0.0) + (row["surface_area_A2"] or 0.0)
    for pair in sorted(union):
        atom_i = _network_atom_metadata(network, pair[0], network.balls.loc[pair[0]])
        atom_j = _network_atom_metadata(network, pair[1], network.balls.loc[pair[1]])
        pair_rows.append(
            {
                "system": str(name or getattr(network, "group_name", "Network")),
                "generator_tuple": pair,
                **{f"atom_i_{key}": value for key, value in atom_i.items()},
                **{f"atom_j_{key}": value for key, value in atom_j.items()},
                "in_power_alpha0": pair in power_pairs,
                "in_aw_alpha0": pair in aw_pairs,
                "power_facet_area_A2": power_area.get(pair),
                "aw_surface_area_A2": aw_area.get(pair),
                "shared": pair in shared,
            }
        )
    summary = {
        "system": str(name or getattr(network, "group_name", "Network")),
        "N_power": len(power_pairs),
        "N_aw": len(aw_pairs),
        "N_shared": len(shared),
        "N_power_only": len(power_pairs - aw_pairs),
        "N_aw_only": len(aw_pairs - power_pairs),
        "N_union": len(union),
        "Jaccard": len(shared) / len(union) if union else 1.0,
        "power_area_finite_pair_count": sum(
            value is not None for value in power_area.values()
        ),
        "power_area_unresolved_pair_count": sum(
            value is None for value in power_area.values()
        ),
        "aw_surface_count": aw_summary["mapped_surface_count"],
        "aw_total_area_A2": aw_summary["total_surface_area_A2"],
        "same_coordinates_and_radii": True,
        "power_method": "existing weighted GUDHI Power alpha=0 bicolor facets (no M=5 selection)",
        "aw_method": "experimental AW restricted-cell nerve alpha=0",
    }
    return pair_rows, summary
