
import os
import sys
import datetime
import tkinter as tk
from tkinter import filedialog
from tkinter import simpledialog

# Add the project root directory to the Python path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../"))
sys.path.append(project_root)


from vorpy.src.analyze.tools.CleanData.CombineLogs2 import read_logs_sections, write_logs_sections, parse_string_lists


def parse_ndxs_spec(spec):
    """
    Parse a user spec into a set of integer atom indices.

    Supported:
    - "10-20" (inclusive)
    - "10:20" (inclusive)
    - "10,11,12"
    - "10 11 12"
    """
    if spec is None:
        return set()

    s = str(spec).strip()
    if s == "":
        return set()

    if "-" in s or ":" in s:
        sep = "-" if "-" in s else ":"
        parts = [_.strip() for _ in s.split(sep) if _.strip() != ""]
        if len(parts) != 2:
            raise ValueError(f"Bad range spec: {spec!r}")
        a = int(parts[0])
        b = int(parts[1])
        if a <= b:
            return set(range(a, b + 1))
        else:
            return set(range(b, a + 1))

    # list form
    vals = []
    for tok in s.replace(",", " ").split():
        vals.append(int(tok))

    return set(vals)


def _safe_int(x):
    try:
        return int(str(x).strip())
    except Exception:
        return None


def _safe_float(x):
    try:
        return float(str(x).strip())
    except Exception:
        return None


def _filter_atoms(atom_rows, ndxs):
    keep = []

    for row in atom_rows:
        if not row:
            continue
        idx = _safe_int(row[0])
        if idx is None:
            continue
        if idx in ndxs:
            keep.append(row)

    return keep


def _filter_surfaces(surface_rows, ndxs, keep_boundary=False):
    """
    Surface rows: [Index, Ball 1, Ball 2, ...]
    """
    keep = []

    for row in surface_rows:
        if not row:
            continue
        b1 = _safe_int(row[1]) if len(row) > 1 else None
        b2 = _safe_int(row[2]) if len(row) > 2 else None
        if b1 is None or b2 is None:
            continue

        if keep_boundary:
            if (b1 in ndxs) or (b2 in ndxs):
                keep.append(row)
        else:
            if (b1 in ndxs) and (b2 in ndxs):
                keep.append(row)

    return keep


def _filter_edges(edge_rows, ndxs):
    """
    Edge rows: [Index, Ball 1, Ball 2, Ball 3, ...]
    """
    keep = []

    for row in edge_rows:
        if not row:
            continue

        b1 = _safe_int(row[1]) if len(row) > 1 else None
        b2 = _safe_int(row[2]) if len(row) > 2 else None
        b3 = _safe_int(row[3]) if len(row) > 3 else None

        if b1 is None or b2 is None or b3 is None:
            continue

        if (b1 in ndxs) and (b2 in ndxs) and (b3 in ndxs):
            keep.append(row)

    return keep


def _filter_vertices(vertex_rows, ndxs):
    """
    Vertex rows: [Index, Ball 1, Ball 2, Ball 3, Ball 4, ...]
    """
    keep = []

    for row in vertex_rows:
        if not row:
            continue

        b1 = _safe_int(row[1]) if len(row) > 1 else None
        b2 = _safe_int(row[2]) if len(row) > 2 else None
        b3 = _safe_int(row[3]) if len(row) > 3 else None
        b4 = _safe_int(row[4]) if len(row) > 4 else None

        if b1 is None or b2 is None or b3 is None or b4 is None:
            continue

        if (b1 in ndxs) and (b2 in ndxs) and (b3 in ndxs) and (b4 in ndxs):
            keep.append(row)

    return keep


def _recompute_group_information(group_section, atom_rows):
    """
    Recompute a *minimal* group information line from the filtered atom rows.

    We keep the section formatting intact, but update:
    Volume, Surface Area, Mass, Density, Center of Mass, VDW Volume, VDW Center of Mass (if present).

    Any additional columns in the original group information row are preserved if possible.
    """
    if group_section is None:
        return

    if group_section.get("header") is None or not group_section.get("rows"):
        return

    header = group_section["header"]
    row = group_section["rows"][0]

    # Map titles to column index
    col = {name: i for i, name in enumerate(header)}

    # Atom header format (from CombineLogs): X,Y,Z at 6..8; Mass at 5; Volume at 10; VDW Vol at 11; SA at 12.
    vols = []
    sas = []
    masses = []
    vdw_vols = []
    locs = []

    for a in atom_rows:
        v = _safe_float(a[10]) if len(a) > 10 else None
        sa = _safe_float(a[12]) if len(a) > 12 else None
        m = _safe_float(a[5]) if len(a) > 5 else None
        vdw = _safe_float(a[11]) if len(a) > 11 else None

        x = _safe_float(a[6]) if len(a) > 6 else None
        y = _safe_float(a[7]) if len(a) > 7 else None
        z = _safe_float(a[8]) if len(a) > 8 else None

        if x is not None and y is not None and z is not None:
            locs.append((x, y, z))
        else:
            locs.append((0.0, 0.0, 0.0))

        vols.append(0.0 if v is None else v)
        sas.append(0.0 if sa is None else sa)
        masses.append(0.0 if m is None else m)
        vdw_vols.append(0.0 if vdw is None else vdw)

    total_vol = sum(vols)
    total_sa = sum(sas)
    total_mass = sum(masses)
    total_vdw = sum(vdw_vols)

    def weighted_com(weights):
        wsum = sum(weights)
        if wsum == 0:
            return [0.0, 0.0, 0.0]
        cx = sum(locs[i][0] * weights[i] for i in range(len(locs))) / wsum
        cy = sum(locs[i][1] * weights[i] for i in range(len(locs))) / wsum
        cz = sum(locs[i][2] * weights[i] for i in range(len(locs))) / wsum
        return [cx, cy, cz]

    com = weighted_com(vols)
    vdw_com = weighted_com(vdw_vols)

    density = (total_vdw / total_vol) if total_vol != 0 else 0.0

    # Update row if columns exist
    if "Volume" in col:
        row[col["Volume"]] = total_vol
    if "Surface Area" in col:
        row[col["Surface Area"]] = total_sa
    if "Mass" in col:
        row[col["Mass"]] = total_mass
    if "Density" in col:
        row[col["Density"]] = density
    if "Center of Mass" in col:
        row[col["Center of Mass"]] = [round(_, 6) for _ in com]
    if "VDW Volume" in col:
        row[col["VDW Volume"]] = total_vdw
    if "VDW Center of Mass" in col:
        row[col["VDW Center of Mass"]] = [round(_, 6) for _ in vdw_com]

    group_section["rows"][0] = row


def _touch_build_information(build_section, output_path):
    """
    If possible, update Location and Completion Date in build information.
    """
    if build_section is None:
        return

    if build_section.get("header") is None or not build_section.get("rows"):
        return

    header = build_section["header"]
    row = build_section["rows"][0]
    col = {name: i for i, name in enumerate(header)}

    if "Location" in col:
        row[col["Location"]] = output_path

    if "Completion Date" in col:
        row[col["Completion Date"]] = datetime.datetime.now()

    build_section["rows"][0] = row


def filter_logs_sections(sections, ndxs):
    """
    Filter parsed logs sections to only keep the atoms in ndxs, and any
    Surfaces/Edges/Vertices fully supported by those atoms.

    Notes
    -----
    - Surfaces are kept only if *both* Ball 1 and Ball 2 are in ndxs.
    - Edges are kept only if Ball 1/2/3 are all in ndxs.
    - Vertices are kept only if Ball 1/2/3/4 are all in ndxs.
    """
    out = sections.copy()

    if "Atoms" in out:
        out["Atoms"] = dict(out["Atoms"])
        out["Atoms"]["rows"] = _filter_atoms(out["Atoms"].get("rows", []), ndxs)

    if "Surfaces" in out:
        out["Surfaces"] = dict(out["Surfaces"])
        out["Surfaces"]["rows"] = _filter_surfaces(out["Surfaces"].get("rows", []), ndxs, keep_boundary=False)

    if "Edges" in out:
        out["Edges"] = dict(out["Edges"])
        out["Edges"]["rows"] = _filter_edges(out["Edges"].get("rows", []), ndxs)

    if "Vertices" in out:
        out["Vertices"] = dict(out["Vertices"])
        out["Vertices"]["rows"] = _filter_vertices(out["Vertices"].get("rows", []), ndxs)

    # Recompute group info if present
    if "group information" in out and "Atoms" in out:
        _recompute_group_information(out["group information"], out["Atoms"]["rows"])

    return out


def extract_logs(input_file=None, ndxs=None, output_dir=None, output_name=None):
    """
    Main API: extract an atom subset from a full logs file.

    Parameters
    ----------
    input_file : str | None
    ndxs : set[int] | None
    output_dir : str | None
    output_name : str | None

    Returns
    -------
    out_path : str
    """
    if input_file is None:
        input_file = filedialog.askopenfilename(title="Select a full logs CSV (e.g., Total_logs.csv)")

    if not input_file or not os.path.exists(input_file):
        raise FileNotFoundError(f"Input logs not found: {input_file!r}")

    if ndxs is None:
        spec = simpledialog.askstring("Atom indices", "Enter atom indices (e.g., 10-50 or 10,11,12):")
        ndxs = parse_ndxs_spec(spec)

    if not ndxs:
        raise ValueError("No atom indices selected (ndxs is empty).")

    if output_dir is None:
        output_dir = filedialog.askdirectory(title="Choose output folder")

    if not output_name:
        output_name = "Extracted_logs.csv"

    out_path = os.path.join(output_dir, output_name)

    sections = read_logs_sections(input_file)
    filtered = filter_logs_sections(sections, ndxs)

    _touch_build_information(filtered.get("build information"), out_path)

    write_logs_sections(out_path, filtered, renumber=True)

    print(f"Wrote extracted logs to: {out_path}")
    return out_path


if __name__ == "__main__":
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes("-topmost", 1)

    extract_logs()
