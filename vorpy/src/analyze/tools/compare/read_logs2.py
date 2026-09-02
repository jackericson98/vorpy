"""Backward-compatible VorPy log reader with modern header-driven parsing.

Key improvement
---------------
Modern VorPy logs contain additional atom/surface curvature and representative
surface-energy fields. Older read_logs2 versions parsed rows positionally,
which silently shifted fields after newly inserted columns.

This reader uses each section's CSV header to map values by column name.
Older logs remain supported through aliases and positional fallbacks.
"""

import ast
import csv
import re
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# Generic converters
# ---------------------------------------------------------------------------

def sort_bool(value):
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() == "true"


def _strip_numpy_wrappers(value):
    """Turn strings like np.float64(1.2) into 1.2 before literal_eval."""
    value = str(value)
    value = re.sub(
        r"(?:np\.)?(?:float16|float32|float64|int16|int32|int64)\(([^()]*)\)",
        r"\1",
        value,
    )
    return value


def parse_string_lists(value, apply_type=float):
    """Parse scalar/list/nested-list CSV strings robustly, including negatives."""
    if value is None:
        return []

    if isinstance(value, (list, tuple)):
        def convert(obj):
            if isinstance(obj, (list, tuple)):
                return [convert(x) for x in obj]
            return apply_type(obj)
        return convert(value)

    text = str(value).strip()
    if text == "":
        return []

    text = _strip_numpy_wrappers(text)

    try:
        parsed = ast.literal_eval(text)
    except (ValueError, SyntaxError):
        # Conservative compatibility fallback.
        numbers = re.findall(
            r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?",
            text,
        )
        return [apply_type(x) for x in numbers]

    def convert(obj):
        if isinstance(obj, (list, tuple)):
            return [convert(x) for x in obj]
        return apply_type(obj)

    if isinstance(parsed, (list, tuple)):
        return convert(parsed)

    return [apply_type(parsed)]


def parse_string_lists_int(value, apply_type=int):
    return parse_string_lists(value, apply_type=apply_type)


def _to_int(value, default=0):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def _to_float(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _to_str(value):
    return "" if value is None else str(value)


def _norm_header(value):
    """Normalize a CSV header enough to survive spelling/case variations."""
    value = str(value).strip().lower()
    value = value.replace("_", " ")
    value = value.replace("-", " ")
    value = re.sub(r"\s+", " ", value)
    return value


def _row_dict(headers, row):
    """Map a CSV row by header, padding missing trailing fields with ''."""
    return {
        _norm_header(header): (row[i] if i < len(row) else "")
        for i, header in enumerate(headers)
    }


def _get(mapped, *names, default=""):
    for name in names:
        key = _norm_header(name)
        if key in mapped and mapped[key] != "":
            return mapped[key]
    return default


# ---------------------------------------------------------------------------
# Atom parsing
# ---------------------------------------------------------------------------

ATOM_CONVERTERS = {
    "Index": _to_int,
    "Name": _to_str,
    "Residue": _to_str,
    "Residue Sequence": _to_int,
    "Chain": _to_str,
    "Mass": _to_float,
    "X": _to_float,
    "Y": _to_float,
    "Z": _to_float,
    "Radius": _to_float,
    "Volume": _to_float,
    "Van Der Waals Volume": _to_float,
    "Surface Area": _to_float,
    "Complete Cell?": sort_bool,

    # Modern curvature fields
    "Maximum Mean Curvature": _to_float,
    "Average Mean Surface Curvature": _to_float,
    "Maximum Gaussian Curvature": _to_float,
    "Average Gaussian Surface Curvature": _to_float,
    "Integrated Mean Curvature": _to_float,
    "Integrated Mean Curvature Squared": _to_float,
    "Integrated Gaussian Curvature": _to_float,
    "Representative Surface Energy": _to_float,

    # Older curvature aliases retained below too
    "Maximum Curvature": _to_float,
    "Average Surface Curvature": _to_float,

    "Sphericity": _to_float,
    "Isometric Quotient": _to_float,
    "Inner Ball?": sort_bool,
    "Number of Neighbors": _to_int,
    "Closest Neighbor": _to_int,
    "Closest Neighbor Distance": _to_float,
    "Layer Distance Average": parse_string_lists,
    "Layer Distance RMSD": parse_string_lists,
    "Minimum Point Distance": _to_float,
    "Maximum Point Distance": _to_float,
    "Number of Overlaps": _to_int,
    "Contact Area": _to_float,
    "Non-Overlap Volume": _to_float,
    "Overlap Volume": _to_float,
    "Center of Mass": parse_string_lists,
    "Moment of Inertia Tensor": parse_string_lists,
    "Bounding Box": parse_string_lists,
    "Neighbors": parse_string_lists_int,
}


ATOM_ALIASES = {
    "Van Der Waals Volume": ["VDW Volume"],
    "Non-Overlap Volume": ["Non - Overlap Volume", "Non Overlap Volume"],
    "Neighbors": ["neighbors"],
}


def read_atom(atom_line, headers=None):
    """Parse an atom row, preferably from its actual CSV header."""
    if headers is None:
        # Legacy fallback layouts.
        modern_legacy_headers = [
            "Index", "Name", "Residue", "Residue Sequence", "Chain", "Mass",
            "X", "Y", "Z", "Radius", "Volume", "Van Der Waals Volume",
            "Surface Area", "Complete Cell?", "Maximum Mean Curvature",
            "Average Mean Surface Curvature", "Maximum Gaussian Curvature",
            "Average Gaussian Surface Curvature", "Sphericity",
            "Isometric Quotient", "Inner Ball?", "Number of Neighbors",
            "Closest Neighbor", "Closest Neighbor Distance",
            "Layer Distance Average", "Layer Distance RMSD",
            "Minimum Point Distance", "Maximum Point Distance",
            "Number of Overlaps", "Contact Area", "Non-Overlap Volume",
            "Overlap Volume", "Center of Mass", "Moment of Inertia Tensor",
            "Bounding Box", "Neighbors",
        ]

        old_headers = [
            "Index", "Name", "Residue", "Residue Sequence", "Chain", "Mass",
            "X", "Y", "Z", "Radius", "Volume", "Van Der Waals Volume",
            "Surface Area", "Complete Cell?", "Maximum Curvature",
            "Average Surface Curvature", "Sphericity", "Isometric Quotient",
            "Inner Ball?", "Number of Neighbors", "Closest Neighbor",
            "Closest Neighbor Distance", "Layer Distance Average",
            "Layer Distance RMSD", "Minimum Point Distance",
            "Maximum Point Distance", "Number of Overlaps", "Contact Area",
            "Non-Overlap Volume", "Overlap Volume", "Center of Mass",
            "Moment of Inertia Tensor", "Bounding Box", "Neighbors",
        ]

        headers = (
            modern_legacy_headers
            if len(atom_line) >= len(modern_legacy_headers)
            else old_headers
        )

    mapped = _row_dict(headers, atom_line)
    atom = {}

    for canonical, converter in ATOM_CONVERTERS.items():
        aliases = ATOM_ALIASES.get(canonical, [])
        raw = _get(mapped, canonical, *aliases, default="")
        if raw == "":
            continue
        atom[canonical] = converter(raw)

    # Friendly aliases expected by some older analysis code.
    if "Name" in atom:
        atom["name"] = atom["Name"]

    if "Radius" in atom:
        atom["rad"] = atom["Radius"]

    if all(k in atom for k in ("X", "Y", "Z")):
        atom["loc"] = [atom["X"], atom["Y"], atom["Z"]]

    # Normalize old curvature names into modern names where possible.
    if "Maximum Mean Curvature" not in atom and "Maximum Curvature" in atom:
        atom["Maximum Mean Curvature"] = atom["Maximum Curvature"]

    if (
        "Average Mean Surface Curvature" not in atom
        and "Average Surface Curvature" in atom
    ):
        atom["Average Mean Surface Curvature"] = atom["Average Surface Curvature"]

    return atom


# ---------------------------------------------------------------------------
# Surface parsing
# ---------------------------------------------------------------------------

def read_surf(surf_line, headers=None):
    """
    Parse old and modern surface rows.

    Modern 3.5.x example header:
        Index
        Ball 1
        Ball 2
        Surface Area
        Mean Curvature
        Average Mean Curvature
        Gaussian Curvature
        Average Gaussian Curvature
        Integrated Mean Curvature
        Integrated Mean Curvature Squared
        Integrated Gaussian Curvature
        Representative Surface Energy
        Ball 1 Volume Contribution
        Ball 2 Volume Contribution
        Contact Area
        Overlap
    """
    while surf_line and surf_line[-1] == "":
        surf_line = surf_line[:-1]

    if not surf_line:
        return None

    if headers is not None:
        mapped = _row_dict(headers, surf_line)

        surf = {
            "Index": _to_int(_get(mapped, "Index")),
            "Balls": [
                _to_int(_get(mapped, "Ball 1")),
                _to_int(_get(mapped, "Ball 2")),
            ],
            "Surface Area": _to_float(_get(mapped, "Surface Area")),
        }

        optional_float_fields = {
            "Mean Curvature": ["Mean Curvature", "Curvature"],
            "Average Mean Curvature": ["Average Mean Curvature"],
            "Gauss Curvature": ["Gaussian Curvature", "Gauss Curvature"],
            "Average Gauss Curvature": [
                "Average Gaussian Curvature",
                "Average Gauss Curvature",
            ],
            "Integrated Mean Curvature": ["Integrated Mean Curvature"],
            "Integrated Mean Curvature Squared": [
                "Integrated Mean Curvature Squared"
            ],
            "Integrated Gaussian Curvature": [
                "Integrated Gaussian Curvature"
            ],
            "Representative Surface Energy": [
                "Representative Surface Energy"
            ],
            "Contact Area": ["Contact Area"],
            "Overlap": ["Overlap"],
        }

        for canonical, aliases in optional_float_fields.items():
            raw = _get(mapped, *aliases, default="")
            if raw != "":
                surf[canonical] = _to_float(raw)

        b1v = _get(
            mapped,
            "Ball 1 Volume Contribution",
            "Ball 1 Volume",
            default="",
        )
        b2v = _get(
            mapped,
            "Ball 2 Volume Contribution",
            "Ball 2 Volume",
            default="",
        )

        volumes = []
        if b1v != "":
            volumes.append(_to_float(b1v))
        if b2v != "":
            volumes.append(_to_float(b2v))
        surf["Ball Volumes"] = volumes

        surf.setdefault("Contact Area", 0.0)
        surf.setdefault("Overlap", 0.0)

        return surf

    # Positional fallback for legacy files without a captured header.
    n = len(surf_line)

    if n >= 16:
        return {
            "Index": _to_int(surf_line[0]),
            "Balls": [_to_int(surf_line[1]), _to_int(surf_line[2])],
            "Surface Area": _to_float(surf_line[3]),
            "Mean Curvature": _to_float(surf_line[4]),
            "Average Mean Curvature": _to_float(surf_line[5]),
            "Gauss Curvature": _to_float(surf_line[6]),
            "Average Gauss Curvature": _to_float(surf_line[7]),
            "Integrated Mean Curvature": _to_float(surf_line[8]),
            "Integrated Mean Curvature Squared": _to_float(surf_line[9]),
            "Integrated Gaussian Curvature": _to_float(surf_line[10]),
            "Representative Surface Energy": _to_float(surf_line[11]),
            "Ball Volumes": [
                _to_float(surf_line[12]),
                _to_float(surf_line[13]),
            ],
            "Contact Area": _to_float(surf_line[14]),
            "Overlap": _to_float(surf_line[15]),
        }

    if n == 12:
        return {
            "Index": _to_int(surf_line[0]),
            "Balls": [_to_int(surf_line[1]), _to_int(surf_line[2])],
            "Surface Area": _to_float(surf_line[3]),
            "Mean Curvature": _to_float(surf_line[4]),
            "Average Mean Curvature": _to_float(surf_line[5]),
            "Gauss Curvature": _to_float(surf_line[6]),
            "Average Gauss Curvature": _to_float(surf_line[7]),
            "Ball Volumes": [
                _to_float(surf_line[8]),
                _to_float(surf_line[9]),
            ],
            "Contact Area": _to_float(surf_line[10]),
            "Overlap": _to_float(surf_line[11]),
        }

    if n == 10:
        return {
            "Index": _to_int(surf_line[0]),
            "Balls": [_to_int(surf_line[1]), _to_int(surf_line[2])],
            "Surface Area": _to_float(surf_line[3]),
            "Mean Curvature": _to_float(surf_line[4]),
            "Gauss Curvature": _to_float(surf_line[5]),
            "Ball Volumes": [
                _to_float(surf_line[6]),
                _to_float(surf_line[7]),
            ],
            "Contact Area": _to_float(surf_line[8]),
            "Overlap": _to_float(surf_line[9]),
        }

    if n == 9:
        return {
            "Index": _to_int(surf_line[0]),
            "Balls": [_to_int(surf_line[1]), _to_int(surf_line[2])],
            "Surface Area": _to_float(surf_line[3]),
            "Mean Curvature": _to_float(surf_line[4]),
            "Gauss Curvature": _to_float(surf_line[5]),
            "Ball Volumes": [
                _to_float(surf_line[6]),
                _to_float(surf_line[7]),
            ],
            "Contact Area": _to_float(surf_line[8]),
            "Overlap": 0.0,
        }

    if n == 8:
        return {
            "Index": _to_int(surf_line[0]),
            "Balls": [_to_int(surf_line[1]), _to_int(surf_line[2])],
            "Surface Area": _to_float(surf_line[3]),
            "Curvature": _to_float(surf_line[4]),
            "Ball Volumes": [
                _to_float(surf_line[5]),
                _to_float(surf_line[6]),
            ],
            "Contact Area": _to_float(surf_line[7]),
            "Overlap": 0.0,
        }

    return None


# ---------------------------------------------------------------------------
# Edge / vertex parsing
# ---------------------------------------------------------------------------

def read_edge(edge_line, headers=None):
    if headers is not None:
        mapped = _row_dict(headers, edge_line)
        return {
            "Index": _to_int(_get(mapped, "Index")),
            "Balls": [
                _to_int(_get(mapped, "Ball 1")),
                _to_int(_get(mapped, "Ball 2")),
                _to_int(_get(mapped, "Ball 3")),
            ],
            "Length": _to_float(_get(mapped, "Length")),
        }

    return {
        "Index": _to_int(edge_line[0]),
        "Balls": [_to_int(x) for x in edge_line[1:4]],
        "Length": _to_float(edge_line[4]),
    }


def read_vert(vert_line, headers=None):
    if headers is not None:
        mapped = _row_dict(headers, vert_line)
        return {
            "Index": _to_int(_get(mapped, "Index")),
            "Balls": [
                _to_int(_get(mapped, "Ball 1")),
                _to_int(_get(mapped, "Ball 2")),
                _to_int(_get(mapped, "Ball 3")),
                _to_int(_get(mapped, "Ball 4")),
            ],
            "loc": [
                _to_float(_get(mapped, "x")),
                _to_float(_get(mapped, "y")),
                _to_float(_get(mapped, "z")),
            ],
            "rad": _to_float(_get(mapped, "r", "Radius")),
        }

    return {
        "Index": _to_int(vert_line[0]),
        "Balls": [_to_int(x) for x in vert_line[1:5]],
        "loc": [_to_float(x) for x in vert_line[5:8]],
        "rad": _to_float(vert_line[8]),
    }


# ---------------------------------------------------------------------------
# Main log reader
# ---------------------------------------------------------------------------

SECTION_NAMES = {
    "build information",
    "build informaiton",   # historical typo in current logs
    "group information",
    "atoms",
    "surfaces",
    "edges",
    "vertices",
}


def _is_section(line):
    return bool(line) and _norm_header(line[0]) in SECTION_NAMES


def _parse_build(headers, row):
    mapped = _row_dict(headers, row)

    def get(*keys, default=None):
        return _get(mapped, *keys, default=default)

    return {
        "name": get("Name", default=""),
        "location": get("Location"),
        "time": get("Completion Date", "Completion Time", "Time"),
        "network_type": get("Network Type", default=""),
        "surface_resolution": _to_float(get("Surface Resolution", default=0.0)),
        "box_size": _to_float(get("Box Size", default=0.0)),
        "max_vert": _to_float(
            get("Maximum Allowable Vertex", "Max Vert", default=0.0)
        ),
        "Total_Time": _to_float(get("Total Time", default=0.0)),
        "vert_time": _to_float(get("Vertex Time", default=0.0)),
        "connect_time": _to_float(get("Connect Time", default=0.0)),
        "surf_time": _to_float(
            get("Surface Building Time", "Surface Time", default=0.0)
        ),
        "analysis_time": _to_float(
            get("Analysis time", "Analysis Time", default=0.0)
        ),
        "max_vertex": _to_float(
            get("Maximum Found Vertex", "Max Vertex", default=0.0)
        ),
        "version": get(
            "vorPy version",
            "VorPy Version",
            "Version",
            default="< 3.2.0",
        ),
    }


def _parse_group(headers, row):
    mapped = _row_dict(headers, row)

    group = {}
    converters = {
        "Name": _to_str,
        "Volume": _to_float,
        "Surface Area": _to_float,
        "Mass": _to_float,
        "Density": _to_float,
        "Center of Mass": parse_string_lists,
        "VDW Volume": _to_float,
        "VDW Center of Mass": parse_string_lists,
        "Moment of Inertia": parse_string_lists,
        "Spatial Moment of Inertia": parse_string_lists,
    }

    for key, converter in converters.items():
        raw = _get(mapped, key, default="")
        if raw != "":
            group[key] = converter(raw)

    return group


def read_logs2(
    log_files,
    return_dict=False,
    no_sol=False,
    all_=True,
    balls=False,
    surfs=False,
    edges=False,
    verts=False,
):
    """
    Read one VorPy log or a list of logs.

    The public return structure matches the historical read_logs2 API.
    """
    one_file = isinstance(log_files, (str, bytes))
    if one_file:
        log_files = [log_files]

    file_info = {}

    for file in log_files:
        with open(file, "r", newline="", encoding="utf-8-sig") as logs:
            rows = list(csv.reader(logs))

        data = {}
        group_data = {}
        atoms = []
        surf_list = []
        edge_list = []
        vert_list = []

        current_section = None
        current_headers = None
        expect_header = False

        for line in rows:
            if not line or all(str(x).strip() == "" for x in line):
                continue

            first = _norm_header(line[0])

            if first in SECTION_NAMES:
                current_section = first
                current_headers = None
                expect_header = True
                continue

            if expect_header:
                current_headers = line
                expect_header = False
                continue

            if current_section in {"build information", "build informaiton"}:
                data = _parse_build(current_headers, line)
                current_section = None
                continue

            if current_section == "group information":
                group_data = _parse_group(current_headers, line)
                current_section = None
                continue

            if current_section == "atoms" and (all_ or balls):
                atom = read_atom(line, headers=current_headers)

                if no_sol:
                    name = atom.get("Name", "").strip().lower()
                    residue = atom.get("Residue", "").strip().lower()

                    if (
                        name in {
                            "hw1", "hw2", "ow", "h02", "h01",
                            "na", "cl", "mg", "k"
                        }
                        or residue in {
                            "sol", "hoh", "wat", "h2o",
                            "tip3", "tip3p", "spc", "spce"
                        }
                    ):
                        continue

                atoms.append(atom)
                continue

            if current_section == "surfaces" and (all_ or surfs):
                surf = read_surf(line, headers=current_headers)
                if surf is not None:
                    surf_list.append(surf)
                continue

            if current_section == "edges" and (all_ or edges):
                edge_list.append(read_edge(line, headers=current_headers))
                continue

            if current_section == "vertices" and (all_ or verts):
                vert_list.append(read_vert(line, headers=current_headers))
                continue

        file_name = data.get("name") or Path(file).stem
        unique_name = file_name
        index = 0

        while unique_name in file_info:
            unique_name = f"{file_name}{index}"
            index += 1

        if return_dict:
            file_info[unique_name] = {
                "data": data,
                "group data": group_data,
                "atoms": atoms,
                "surfs": surf_list,
                "edges": edge_list,
                "verts": vert_list,
            }
        else:
            file_info[unique_name] = {
                "data": data,
                "group data": group_data,
                "atoms": pd.DataFrame(atoms),
                "surfs": pd.DataFrame(surf_list),
                "edges": pd.DataFrame(edge_list),
                "verts": pd.DataFrame(vert_list),
            }

    if one_file:
        first_key = next(iter(file_info))
        return file_info[first_key]

    return file_info
