"""Shared prepared-mesh representation and file-format writers for VorPy."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Optional
import xml.etree.ElementTree as ET

import numpy as np


SUPPORTED_MESH_FORMATS = ("off", "ply", "vtp")


@dataclass
class MeshData:
    """Geometry prepared by a VorPy component exporter.

    Colors are stored per face in normalized RGB form. ``face_data`` contains
    optional numeric scientific fields with one value per triangle.
    """

    points: np.ndarray
    triangles: np.ndarray
    face_colors: Optional[np.ndarray] = None
    face_data: Dict[str, np.ndarray] = field(default_factory=dict)

    def __post_init__(self):
        self.points = np.asarray(self.points, dtype=float).reshape((-1, 3))
        self.triangles = np.asarray(self.triangles, dtype=np.int64).reshape((-1, 3))

        if self.face_colors is not None:
            colors = np.asarray(self.face_colors, dtype=float)
            if colors.ndim != 2 or colors.shape[0] != len(self.triangles) or colors.shape[1] < 3:
                raise ValueError("face_colors must contain one RGB or RGBA color per triangle")
            self.face_colors = np.clip(colors[:, :3], 0.0, 1.0)

        for name, values in self.face_data.items():
            array = np.asarray(values)
            if array.ndim != 1 or len(array) != len(self.triangles):
                raise ValueError(f"face_data[{name!r}] must contain one value per triangle")
            self.face_data[name] = array


def _normalize_face_colors(colors, triangle_count, part_index):
    """Return one RGB row per triangle for a single prepared mesh part.

    Surface exporters can legitimately produce no triangles for an individual
    mesh part.  Some callers may also represent a uniform color as one RGB
    triplet instead of repeating that triplet for every triangle.  Normalize
    both cases here, where the triangle count is known, so MeshData always
    receives a strict (N, 3) face-color array.
    """
    triangle_count = int(triangle_count)
    array = np.asarray(colors, dtype=float)

    if triangle_count == 0:
        if array.size not in (0, 3, 4):
            raise ValueError(
                f"color_parts[{part_index}] contains {array.size} color values "
                "for a mesh part with zero triangles"
            )
        return np.empty((0, 3), dtype=float)

    if array.size == 0:
        raise ValueError(
            f"color_parts[{part_index}] is empty but the mesh part contains "
            f"{triangle_count} triangle(s)"
        )

    # A single RGB/RGBA triplet represents one uniform color for the part.
    if array.ndim == 1 and array.size in (3, 4):
        rgb = array[:3].reshape((1, 3))
        return np.repeat(rgb, triangle_count, axis=0)

    # Accept a flat sequence containing one RGB/RGBA value per triangle.
    if array.ndim == 1:
        if array.size == triangle_count * 3:
            return array.reshape((triangle_count, 3))
        if array.size == triangle_count * 4:
            return array.reshape((triangle_count, 4))[:, :3]
        raise ValueError(
            f"color_parts[{part_index}] has shape {array.shape}; expected one "
            f"RGB/RGBA color or {triangle_count} per-triangle colors"
        )

    if array.ndim != 2 or array.shape[1] < 3:
        raise ValueError(
            f"color_parts[{part_index}] has shape {array.shape}; expected "
            "(N, 3+) face colors"
        )

    if array.shape[0] == 1 and triangle_count > 1:
        return np.repeat(array[:, :3], triangle_count, axis=0)

    if array.shape[0] != triangle_count:
        raise ValueError(
            f"color_parts[{part_index}] contains {array.shape[0]} colors for "
            f"{triangle_count} triangles"
        )

    return array[:, :3]


def combine_mesh_parts(point_parts: Iterable, triangle_parts: Iterable,
                       color_parts: Optional[Iterable] = None,
                       face_data_parts: Optional[Dict[str, Iterable]] = None) -> MeshData:
    """Combine independently indexed mesh parts into one globally indexed mesh."""
    point_parts = [np.asarray(points, dtype=float).reshape((-1, 3)) for points in point_parts]
    triangle_parts = [np.asarray(tris, dtype=np.int64).reshape((-1, 3)) for tris in triangle_parts]

    if len(point_parts) != len(triangle_parts):
        raise ValueError(
            "point_parts and triangle_parts must contain the same number of mesh parts"
        )

    triangles = []
    offset = 0
    for points, tris in zip(point_parts, triangle_parts):
        triangles.append(tris + offset)
        offset += len(points)

    points = np.concatenate(point_parts) if point_parts else np.empty((0, 3), dtype=float)
    triangles = np.concatenate(triangles) if triangles else np.empty((0, 3), dtype=np.int64)

    colors = None
    if color_parts is not None:
        color_parts = list(color_parts)
        if len(color_parts) != len(triangle_parts):
            raise ValueError(
                "color_parts and triangle_parts must contain the same number of mesh parts"
            )

        parts = [
            _normalize_face_colors(colors_part, len(tris), part_index)
            for part_index, (colors_part, tris) in enumerate(zip(color_parts, triangle_parts))
        ]
        colors = np.concatenate(parts, axis=0) if parts else np.empty((0, 3), dtype=float)

    combined_data = {}
    for name, parts in (face_data_parts or {}).items():
        arrays = [np.asarray(part) for part in parts]
        combined_data[name] = np.concatenate(arrays) if arrays else np.empty(0)

    return MeshData(points, triangles, colors, combined_data)


def write_mesh(mesh: MeshData, file_name, file_type="off", directory=None, chunk_size=10000):
    """Write prepared geometry using one of VorPy's supported mesh formats."""
    file_type = str(file_type).lower().lstrip(".")
    if file_type not in SUPPORTED_MESH_FORMATS:
        supported = ", ".join(SUPPORTED_MESH_FORMATS)
        raise ValueError(f"Unsupported mesh format {file_type!r}. Supported formats: {supported}")

    path = Path(directory or ".") / str(file_name)
    if path.suffix.lower() != f".{file_type}":
        path = Path(str(path) + f".{file_type}")
    path.parent.mkdir(parents=True, exist_ok=True)

    if file_type == "off":
        _write_off(mesh, path, chunk_size)
    elif file_type == "ply":
        _write_ply(mesh, path, chunk_size)
    else:
        _write_vtp(mesh, path)
    return path


def _write_off(mesh, path, chunk_size):
    with path.open("w", buffering=1024 * 1024, newline="\n") as output:
        output.write(f"OFF\n{len(mesh.points)} {len(mesh.triangles)} 0\n\n\n")
        _write_lines(output, (
            f"{round(float(x), 4)} {round(float(y), 4)} {round(float(z), 4)}\n"
            for x, y, z in mesh.points
        ), chunk_size)

        if mesh.face_colors is None:
            lines = (f"3 {a} {b} {c}\n" for a, b, c in mesh.triangles)
        else:
            lines = (
                f"3 {a} {b} {c} {r:g} {g:g} {b_col:g}\n"
                for (a, b, c), (r, g, b_col) in zip(mesh.triangles, mesh.face_colors)
            )
        _write_lines(output, lines, chunk_size)


def _write_ply(mesh, path, chunk_size):
    has_colors = mesh.face_colors is not None
    with path.open("w", buffering=1024 * 1024, newline="\n") as output:
        output.write("ply\nformat ascii 1.0\ncomment Generated by VorPy\n")
        output.write(f"element vertex {len(mesh.points)}\nproperty float x\nproperty float y\nproperty float z\n")
        output.write(f"element face {len(mesh.triangles)}\nproperty list uchar int vertex_indices\n")
        if has_colors:
            output.write("property uchar red\nproperty uchar green\nproperty uchar blue\n")
        output.write("end_header\n")
        _write_lines(output, (f"{x:.4f} {y:.4f} {z:.4f}\n" for x, y, z in mesh.points), chunk_size)

        if has_colors:
            rgb = np.rint(mesh.face_colors * 255).astype(np.uint8)
            lines = (
                f"3 {a} {b} {c} {r} {g} {b_col}\n"
                for (a, b, c), (r, g, b_col) in zip(mesh.triangles, rgb)
            )
        else:
            lines = (f"3 {a} {b} {c}\n" for a, b, c in mesh.triangles)
        _write_lines(output, lines, chunk_size)


def _write_vtp(mesh, path):
    root = ET.Element("VTKFile", type="PolyData", version="0.1", byte_order="LittleEndian")
    poly_data = ET.SubElement(root, "PolyData")
    piece = ET.SubElement(poly_data, "Piece", NumberOfPoints=str(len(mesh.points)),
                          NumberOfPolys=str(len(mesh.triangles)))

    points = ET.SubElement(piece, "Points")
    point_array = ET.SubElement(points, "DataArray", type="Float64", NumberOfComponents="3", format="ascii")
    point_array.text = "\n" + " ".join(f"{value:.8g}" for value in mesh.points.ravel()) + "\n"

    polys = ET.SubElement(piece, "Polys")
    connectivity = ET.SubElement(polys, "DataArray", type="Int64", Name="connectivity", format="ascii")
    connectivity.text = "\n" + " ".join(str(value) for value in mesh.triangles.ravel()) + "\n"
    offsets = ET.SubElement(polys, "DataArray", type="Int64", Name="offsets", format="ascii")
    offsets.text = "\n" + " ".join(str(value) for value in range(3, 3 * len(mesh.triangles) + 1, 3)) + "\n"

    cell_data = ET.SubElement(piece, "CellData")
    if mesh.face_colors is not None:
        color_array = ET.SubElement(cell_data, "DataArray", type="UInt8", Name="RGB",
                                    NumberOfComponents="3", format="ascii")
        rgb = np.rint(mesh.face_colors * 255).astype(np.uint8)
        color_array.text = "\n" + " ".join(str(value) for value in rgb.ravel()) + "\n"
    for name, values in mesh.face_data.items():
        data_type = "Int64" if np.issubdtype(values.dtype, np.integer) else "Float64"
        array = ET.SubElement(cell_data, "DataArray", type=data_type, Name=str(name), format="ascii")
        array.text = "\n" + " ".join(str(value) for value in values) + "\n"

    ET.indent(root, space="  ")
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)


def _write_lines(output, lines, chunk_size):
    buffer = []
    for line in lines:
        buffer.append(line)
        if len(buffer) >= chunk_size:
            output.write("".join(buffer))
            buffer.clear()
    if buffer:
        output.write("".join(buffer))
