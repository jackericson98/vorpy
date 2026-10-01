"""CSV and text exports for Apollonius and AW alpha complexes."""

import csv
import json
from pathlib import Path


def _json(value):
    return json.dumps(value, separators=(",", ":"), allow_nan=False)


def _row(simplex):
    return {
        "simplex_id": simplex.simplex_id,
        "generator_ids": _json(simplex.generator_ids),
        "dimension": simplex.dimension,
        "alpha_birth": "" if simplex.alpha_birth is None else f"{simplex.alpha_birth:.12g}",
        "alpha_status": simplex.alpha_status,
        "alpha_reason": simplex.diagnostics.get("alpha", {}).get("reason", ""),
        "complete": simplex.complete,
        "bounded": simplex.bounded,
        "status": simplex.status,
        "num_primal_components": simplex.num_primal_components,
        "primal_feature_ids": _json(simplex.primal_feature_ids),
        "generator_atom_numbers": _json(simplex.diagnostics.get("atom_numbers", ())),
        "boundary_generator_flags": _json(simplex.diagnostics.get("boundary_generator_flags", ())),
        "generator_coordinates": _json(simplex.generator_coordinates),
        "generator_radii": _json(simplex.generator_radii),
    }


def _write_simplex_csv(path, simplices):
    fields = (
        "simplex_id", "generator_ids", "dimension", "alpha_birth", "alpha_status", "alpha_reason",
        "complete", "bounded", "status", "num_primal_components", "primal_feature_ids",
        "generator_atom_numbers", "generator_coordinates", "generator_radii",
        "boundary_generator_flags",
    )
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for simplex in simplices:
            writer.writerow(_row(simplex))


def export_apollonius(complex_, destination, alpha_complex=None):
    """Write per-dimension dual tables, a summary, and optional alpha outputs."""
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    filenames = {0: "apollonius_vertices.csv", 1: "apollonius_edges.csv",
                 2: "apollonius_triangles.csv", 3: "apollonius_tetrahedra.csv"}
    paths = {}
    for dimension, filename in filenames.items():
        path = destination / filename
        _write_simplex_csv(path, complex_.simplices[dimension].values())
        paths[dimension] = path

    summary_path = destination / "apollonius_summary.txt"
    summary_path.write_text(complex_.summary(), encoding="utf-8")
    paths["summary"] = summary_path

    # Optional PyMOL CGO script: its displayed spheres are generator centers,
    # and its cylinders are the straight abstract dual edges. They are not
    # exported as molecular atoms or chemical bonds.
    generator_points = [simplex.generator_coordinates[0]
                        for simplex in complex_.simplices[0].values()]
    edge_points = [simplex.generator_coordinates
                   for simplex in complex_.simplices[1].values()]
    visualization_path = destination / "apollonius_dual.py"
    script = [
        "# Run inside PyMOL with: run apollonius_dual.py",
        "from pymol import cmd",
        "from pymol.cgo import COLOR, CYLINDER, SPHERE",
        f"generators = {generator_points!r}",
        f"edges = {edge_points!r}",
        "objects = []",
        "for point in generators:",
        "    objects.extend([COLOR, 0.15, 0.55, 0.85, SPHERE, *point, 0.18])",
        "for edge in edges:",
        "    objects.extend([CYLINDER, *edge[0], *edge[1], 0.035, 0.95, 0.68, 0.18, 0.95, 0.68, 0.18])",
        "cmd.delete('apollonius_dual')",
        "cmd.load_cgo(objects, 'apollonius_dual')",
        "cmd.zoom('apollonius_dual')",
        "",
    ]
    visualization_path.write_text("\n".join(script), encoding="utf-8")
    paths["visualization"] = visualization_path

    if alpha_complex is not None:
        label = f"{alpha_complex.alpha:g}"
        alpha_summary = destination / f"alpha_complex_{label}_summary.txt"
        alpha_summary.write_text(alpha_complex.summary(), encoding="utf-8")
        alpha_table = destination / f"alpha_complex_{label}_simplices.csv"
        _write_simplex_csv(alpha_table, (
            simplex for table in alpha_complex.simplices.values()
            for simplex in table.values()
        ))
        paths["alpha_summary"] = alpha_summary
        paths["alpha_simplices"] = alpha_table
    return paths
