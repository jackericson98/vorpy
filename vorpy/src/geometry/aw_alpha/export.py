"""CSV/JSON exports for the experimental AW restricted-cell nerve."""

from __future__ import annotations

import csv
import json
from pathlib import Path


def _json(value):
    return json.dumps(value, separators=(",", ":"), allow_nan=False, default=str)


def _simplex_row(record):
    return {
        "dimension": record.dimension,
        "generator_tuple": _json(record.generator_tuple),
        "geometric_birth_A": record.geometric_birth,
        "filtration_birth_A": record.filtration_birth,
        "birth_source": record.birth_source,
        "closure_adjustment_A": record.closure_adjustment,
        "supported": record.supported,
        "residual_A": record.residual,
        "notes": record.notes,
    }


def _write_csv(path, rows, fields):
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def export_aw_alpha(
    filtration,
    destination,
    *,
    alpha_values=(0.0,),
    group_a=None,
    group_b=None,
    radius_convention="network-supplied radii",
):
    filtration.calculate_births()
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    records = filtration.all_records()
    simplex_fields = (
        "dimension",
        "generator_tuple",
        "geometric_birth_A",
        "filtration_birth_A",
        "birth_source",
        "closure_adjustment_A",
        "supported",
        "residual_A",
        "notes",
    )
    _write_csv(
        destination / "aw_alpha_simplices.csv",
        [_simplex_row(record) for record in records],
        simplex_fields,
    )

    audit_rows = []
    for record in records:
        audit_rows.append(
            {
                **_simplex_row(record),
                "simplex_id": record.simplex_id,
                "diagnostics_json": _json(record.diagnostics),
            }
        )
    _write_csv(
        destination / "aw_alpha_birth_audit.csv",
        audit_rows,
        ("simplex_id", *simplex_fields, "diagnostics_json"),
    )
    unresolved = [
        row
        for row in audit_rows
        if row["geometric_birth_A"] is None
        or row["filtration_birth_A"] is None
        or not row["supported"]
    ]
    _write_csv(
        destination / "unresolved_births.csv",
        unresolved,
        ("simplex_id", *simplex_fields, "diagnostics_json"),
    )

    alpha0 = filtration.alpha_complex(0.0)
    _write_csv(
        destination / "aw_alpha0_simplices.csv",
        [_simplex_row(record) for record in alpha0.records],
        simplex_fields,
    )
    surface_rows, surface_summary = (
        [],
        {
            "pair_count": 0,
            "mapped_surface_count": 0,
            "interface_atom_count": 0,
            "total_surface_area_A2": 0.0,
            "connected_component_count": 0,
            "unsupported_birth_count": 0,
        },
    )
    if group_a is not None and group_b is not None:
        surface_rows, surface_summary = filtration.bicolor_surface_records(
            0.0, group_a, group_b
        )
    surface_fields = (
        "generator_tuple",
        "surface_id",
        *(
            f"atom_{side}_{field}"
            for side in ("i", "j")
            for field in (
                "atom_index",
                "generator_id",
                "stable_id",
                "metadata_status",
                "chain",
                "residue_number",
                "insertion_code",
                "residue_name",
                "atom_name",
                "radius_A",
                "group",
            )
        ),
        "geometric_birth_A",
        "filtration_birth_A",
        "birth_source",
        "surface_area_A2",
        "surface_integrated_mean_curvature_A",
        "surface_mean_curvature",
        "feature_complete",
        "supported",
        "confidence_status",
        "primal_feature",
    )
    _write_csv(
        destination / "aw_alpha0_bicolor_surfaces.csv", surface_rows, surface_fields
    )

    closure_issues = [
        issue for issue in filtration.issues if issue["kind"] == "closure_adjustment"
    ]
    metadata = {
        "experimental": True,
        "name": filtration.name,
        "definition": "experimental additively weighted restricted-cell nerve",
        "distance": "d_i(x)=||x-p_i||-r_i",
        "sublevel_region": "V_i^AW intersect B(p_i,r_i+alpha)",
        "filtration_units": "angstrom additive radial expansion",
        "alpha_zero_interpretation": "network-supplied spheres without further expansion",
        "radius_convention": radius_convention,
        "power_or_gudhi_alpha_used": False,
        "cazals_M5_applied": False,
        "groups_A": sorted(map(int, group_a or ())),
        "groups_B": sorted(map(int, group_b or ())),
        "tolerance_A": filtration.tolerance,
        "alpha_values_reported": [float(value) for value in alpha_values],
        "simplex_counts_by_dimension": {
            str(d): len(filtration.records[d]) for d in range(4)
        },
        "unresolved_birth_count": len(unresolved),
        "filtration_order_violations": [
            issue
            for issue in filtration.validate_filtration()
            if issue["kind"] == "monotonicity_violation"
        ],
        "closure_adjustments_over_tolerance": closure_issues,
        "nerve_good_cover_or_homotopy_equivalence_claimed": False,
        "topology_by_alpha": {},
        "alpha0_bicolor_summary": surface_summary,
    }
    for value in alpha_values:
        subcomplex = filtration.alpha_complex(value)
        metadata["topology_by_alpha"][str(float(value))] = {
            "counts": {str(d): n for d, n in subcomplex.counts.items()},
            "euler_characteristic": subcomplex.euler_characteristic,
            "blocked_simplices": len(subcomplex.blocked),
        }
    (destination / "aw_alpha_metadata.json").write_text(
        json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )

    summary = filtration.summary()
    summary += "\nTopology snapshots\n" + "\n".join(
        f"alpha={float(value):g} A: counts={metadata['topology_by_alpha'][str(float(value))]['counts']}, "
        f"chi={metadata['topology_by_alpha'][str(float(value))]['euler_characteristic']}"
        for value in alpha_values
    )
    summary += (
        f"\nalpha=0 bicolor surfaces: pairs={surface_summary['pair_count']}, "
        f"mapped surfaces={surface_summary['mapped_surface_count']}, "
        f"area={surface_summary['total_surface_area_A2']:.9g} A^2, "
        f"interface atoms={surface_summary['interface_atom_count']}, "
        f"components={surface_summary['connected_component_count']}, "
        f"unsupported pair births={surface_summary['unsupported_birth_count']}\n"
    )
    if closure_issues:
        summary += "\nCLOSURE ADJUSTMENTS ABOVE TOLERANCE\n" + "\n".join(
            f"{issue['simplex']}: +{issue['adjustment_A']:.9g} A"
            for issue in closure_issues
        )
    (destination / "aw_alpha_summary.txt").write_text(summary + "\n", encoding="utf-8")
    return {
        "metadata": metadata,
        "surface_rows": surface_rows,
        "surface_summary": surface_summary,
        "unresolved": unresolved,
    }
