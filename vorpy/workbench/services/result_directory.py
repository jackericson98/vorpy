"""Discover a completed VorPy output directory without importing VorPy."""

from pathlib import Path

from vorpy.workbench.domain import GeometryLayer
from vorpy.workbench.services.structure_loader import load_pdb

MESH_SUFFIXES = {".off", ".ply", ".vtp"}


def load_result_directory(directory: Path):
    pdb_files = sorted(directory.rglob("*.pdb"))
    preferred = [
        path for path in pdb_files
        if not any(token in path.stem.lower() for token in ("surr", "group_atoms", "ext_atoms"))
    ]
    if not preferred and not pdb_files:
        raise ValueError("The selected directory does not contain a PDB structure")
    result = load_pdb((preferred or pdb_files)[0])
    result.name = directory.name

    for mesh_path in sorted(path for path in directory.rglob("*") if path.suffix.lower() in MESH_SUFFIXES):
        lower = mesh_path.stem.lower()
        if "surf" in lower:
            kind, color, opacity = "surfaces", "#4f9fcf", 0.45
        elif "edge" in lower:
            kind, color, opacity = "edges", "#72bde4", 1.0
        elif "vert" in lower:
            kind, color, opacity = "vertices", "#efb84f", 1.0
        else:
            kind, color, opacity = "mesh", "#9f8bd4", 0.75
        relative = mesh_path.relative_to(directory)
        is_atom_cell = "atoms" in {part.lower() for part in relative.parts}
        initially_visible = not is_atom_cell and ("shell" in lower or lower in {"edges", "verts"})
        result.layers.append(GeometryLayer(
            name=str(relative),
            kind=kind,
            source_path=mesh_path,
            color=color,
            opacity=opacity,
            visible=initially_visible,
        ))
    return result


