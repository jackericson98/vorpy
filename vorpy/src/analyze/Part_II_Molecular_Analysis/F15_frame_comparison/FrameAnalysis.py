"""VorPy frame-wise volume, surface area, and number-of-contacts analysis.

Compatibility:
- Group volume/area are read from the Group Information section by default.
- The original cross-network per-atom volume-comparison algorithm is retained
  as `compare_by_atom` and `volume_mode='legacy_atom_sum'`, but is NOT the
  default: atomic records can include environmental/solvent context atoms.
- Accepts f1, f0009 and frame_0009 directories, as well as time-stamped logs
  such as aw_logs(20261027-180302).csv.
- Counts each unique surface with positive Contact Area as one contact.
  This is a count of contact patches, not necessarily distinct atom pairs.

From a command line:
    python FrameAnalysis.py --folder /path/to/frame_parent --no-show
    python FrameAnalysis.py --folder /path/to/frame_parent --no-show
    python FrameAnalysis.py --log /path/to/aw_logs.csv

Running without arguments opens the original Tkinter folder picker.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

# Allow a script nested under a VorPy source checkout to find its project root,
# preserving the original FrameAnalysis.py behavior without hard-coded depth.
for _parent in Path(__file__).resolve().parents:
    if (_parent / "vorpy" / "src").is_dir():
        if str(_parent) not in sys.path:
            sys.path.insert(0, str(_parent))
        break

# Compatible with an installed/local VorPy checkout.
from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


NETWORKS = ("aw", "pow", "prm")
NETWORK_LABELS = {"aw": "AW", "pow": "Power", "prm": "Primitive"}

# Plot appearance: edit these settings to customize all six figures.
NETWORK_COLORS = {"aw": "blue", "pow": "red", "prm": "green"}
TICK_LABEL_SIZE = 14       # Font size of tick numbers on both axes
AXIS_LABEL_SIZE = 16       # Font size of x/y axis labels
SHOW_TITLES = True         # Set False to omit plot titles
TITLE_FONT_SIZE = 18       # Used only when SHOW_TITLES is True

PLOT_LIMITS = {
    "volume_aw_diff": {
        "xlim": (-0.5, 10.5),
        "ylim": (-7.0, 9.0),
    },
    "volume_frame_dev": {
        "xlim": (-0.5, 10.5),
        "ylim": (-11.0, 2.0),
    },
    "surface_aw_diff": {
        "xlim": (-0.5, 10.5),
        "ylim": (-3.0, 2.0),
    },
    "surface_frame_dev": {
        "xlim": (-0.5, 10.5),
        "ylim": (-2.5, 4.0),
    },
    "contacts_aw_diff": {
        "xlim": (-0.5, 10.5),
        "ylim": (-22.0, 12.0),
    },
    "contacts_frame_dev": {
        "xlim": (-0.5, 10.5),
        "ylim": (-8.5, 13.0),
    },
}

FRAME_PATTERN = re.compile(r"(?:f|frame)[_-]?(\d+)$", flags=re.IGNORECASE)
METRICS = {
    "vol": ("Volume", r"Volume ($\mathrm{\AA}^3$)", "volume"),
    "sa": ("Surface area", r"Surface area ($\mathrm{\AA}^2$)", "surface_area"),
    "contacts": ("Number of contacts", "Number of contacts", "number_of_contacts"),
    "contact_count": ("Number of contacts", "Number of contacts", "number_of_contacts"),
}


def _warn(message):
    warnings.warn(message, RuntimeWarning, stacklevel=2)


def _frame_id(name):
    match = FRAME_PATTERN.fullmatch(name)
    return int(match.group(1)) if match else None


def _frame_folders(folder):
    folder = Path(folder).expanduser().resolve()
    if not folder.is_dir():
        raise NotADirectoryError(f"Frame folder does not exist: {folder}")

    frames = []
    for child in folder.iterdir():
        if child.is_dir() and (number := _frame_id(child.name)) is not None:
            frames.append((number, child))

    if frames:
        return sorted(frames)

    # Support pointing the selector directly at a single f9/frame_0009 dir.
    if any((folder / network).is_dir() for network in NETWORKS) or any(
        "_logs" in p.name.lower() and p.suffix.lower() == ".csv"
        for p in folder.iterdir() if p.is_file()
    ):
        return [(_frame_id(folder.name) or 1, folder)]

    raise FileNotFoundError(
        f"No frame directories or network log folders found under {folder}. "
        "Expected f1/aw/aw_logs.csv, frame_0009/aw/aw_logs(...).csv, etc."
    )


def _find_log(frame_folder, network):
    """Prefer historical exact filenames; fall back to timestamped variants."""
    frame_folder = Path(frame_folder)
    for base in (frame_folder / network, frame_folder):
        if not base.is_dir():
            continue
        exact = base / f"{network}_logs.csv"
        if exact.is_file():
            return exact

        candidates = [
            p for p in base.iterdir()
            if p.is_file() and p.suffix.lower() == ".csv"
            and p.stem.lower().startswith(f"{network}_logs")
        ]
        if candidates:
            # Multiple time-stamped runs: prefer newest by modification time.
            return max(candidates, key=lambda p: (p.stat().st_mtime, p.name))
    return None


def _has_contact_area_column(path):
    """Don't interpret the reader's optional-field default as real zero area."""
    with open(path, newline="", encoding="utf-8-sig") as stream:
        reader = csv.reader(stream)
        for row in reader:
            if row and row[0].strip().lower() == "surfaces":
                for header in reader:
                    if header and any(col.strip() for col in header):
                        norm = {x.strip().lower().replace("_", " ") for x in header}
                        return "contact area" in norm
                return False
    return False


def _number(value):
    try:
        value = float(value)
        return value if math.isfinite(value) else math.nan
    except (TypeError, ValueError):
        return math.nan


def _surface_contacts(log, path):
    """Sum contact patch areas once per unique surface (not per atom)."""
    surfaces = log.get("surfs", pd.DataFrame())
    if not _has_contact_area_column(path):
        _warn(f"{path.name}: no Surfaces / Contact Area header; contacts unavailable")
        return math.nan, math.nan
    if surfaces.empty or "Contact Area" not in surfaces.columns:
        _warn(f"{path.name}: no parsed surface contact areas; contacts unavailable")
        return math.nan, math.nan

    if "Index" in surfaces.columns and surfaces["Index"].duplicated().any():
        n = int(surfaces["Index"].duplicated().sum())
        _warn(f"{path.name}: {n} duplicate surface IDs; counting only the first")
        surfaces = surfaces.drop_duplicates(subset="Index")

    areas = pd.to_numeric(surfaces["Contact Area"], errors="coerce")
    if not areas.notna().any():
        return math.nan, math.nan
    if (areas < 0).any():
        _warn(f"{path.name}: negative contact areas were found; check the log")
    return float(areas.sum()), int((areas > 0).sum())


def compare_by_atom(aw_balls, pw_balls, pm_balls, pwd_cut=0.3, pmd_cut=5, eps=1e-12):
    """Retain the OLD cross-network matched-atom algorithm for explicit use.

    WARNING: this intentionally includes matching solvent/context atoms and
    thus does not generally equal the group's geometric volume.
    """
    required = {"Index", "Volume"}
    for label, frame in zip(NETWORKS, (aw_balls, pw_balls, pm_balls)):
        if frame.empty or not required.issubset(frame.columns):
            raise ValueError(f"{label}: missing atom Index or Volume")
        if frame["Index"].duplicated().any():
            raise ValueError(f"{label}: duplicate atom Index values in log")

    pw_lookup = pw_balls.set_index("Index")
    pm_lookup = pm_balls.set_index("Index")
    vols = {"aw": 0.0, "pw": 0.0, "pm": 0.0, "replaced": 0, "matched_atoms": 0}

    for _, aw_ball in aw_balls.iterrows():
        atom_id = aw_ball["Index"]
        if atom_id not in pw_lookup.index or atom_id not in pm_lookup.index:
            continue
        awv = float(aw_ball["Volume"])
        pwv = float(pw_lookup.loc[atom_id, "Volume"])
        pmv = float(pm_lookup.loc[atom_id, "Volume"])
        if not all(math.isfinite(v) for v in (awv, pwv, pmv)):
            continue
        pwd = abs(awv - pwv) / max(abs(awv), eps)
        pmd = abs(awv - pmv) / max(abs(awv), eps)
        if pwd > pwd_cut or pmd > pmd_cut:
            awv = pwv
            vols["replaced"] += 1
        vols["aw"] += awv
        vols["pw"] += pwv
        vols["pm"] += pmv
        vols["matched_atoms"] += 1
    return vols


def get_frames_by_group(folder=None, get_balls=False, volume_mode="group"):
    """Return {frame: {'vol','sa','contacts','contact_count': {aw,pow,prm}}}.

    Parameters
    ----------
    get_balls : bool
        Backwards-compatible argument: load atomic tables for inspection.
        It no longer silently changes the plotted volume to an atom sum.
    volume_mode : {'group', 'legacy_atom_sum'}
        'group' (default) uses the authoritative Group Information volume.
        'legacy_atom_sum' reproduces the original matched-atom/replacement
        comparison, including solvent, for retrospective comparisons only.
    """
    if volume_mode not in {"group", "legacy_atom_sum"}:
        raise ValueError("volume_mode must be 'group' or 'legacy_atom_sum'")

    if folder is None:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        try:
            folder = filedialog.askdirectory(title="Select VorPy frames folder")
        finally:
            root.destroy()
        if not folder:
            return {}

    frame_dict = {}
    for frame_number, frame_folder in _frame_folders(folder):
        values = {key: {network: math.nan for network in NETWORKS}
                  for key in ("vol", "sa", "contacts", "contact_count")}
        logs = {}
        for network in NETWORKS:
            path = _find_log(frame_folder, network)
            if path is None:
                _warn(f"Frame {frame_number}: missing {network} log under {frame_folder}")
                continue
            try:
                record = read_logs2(str(path), all_=False,
                                    balls=(get_balls or volume_mode == "legacy_atom_sum"),
                                    surfs=True)
            except (OSError, ValueError, KeyError) as exc:
                _warn(f"Frame {frame_number}/{network}: could not parse {path}: {exc}")
                continue

            logs[network] = record
            group = record.get("group data", {})
            values["vol"][network] = _number(group.get("Volume"))
            values["sa"][network] = _number(group.get("Surface Area"))
            _, contact_count = _surface_contacts(record, path)
            values["contacts"][network] = contact_count
            values["contact_count"][network] = contact_count
            if not math.isfinite(values["vol"][network]):
                # Old logs may omit Group Information: preserve a documented
                # atom-sum fallback rather than silently treating it as zero.
                if not get_balls and volume_mode != "legacy_atom_sum":
                    record = read_logs2(str(path), all_=False, balls=True)
                atoms = record.get("atoms", pd.DataFrame())
                if "Volume" in atoms.columns and not atoms.empty:
                    values["vol"][network] = float(
                        pd.to_numeric(atoms["Volume"], errors="coerce").sum()
                    )
                    _warn(f"Frame {frame_number}/{network}: no Group Volume; "
                          "falling back to total ATOM volumes, which may include solvent")
                else:
                    _warn(f"Frame {frame_number}/{network}: no usable volume")

        if volume_mode == "legacy_atom_sum":
            if all(network in logs for network in NETWORKS):
                comparison = compare_by_atom(
                    *(logs[network]["atoms"] for network in NETWORKS)
                )
                for network, key in (("aw", "aw"), ("pow", "pw"), ("prm", "pm")):
                    values["vol"][network] = comparison[key]
                print(f"Frame {frame_number}: legacy matched-atom volume "
                      f"({comparison['matched_atoms']} atoms, "
                      f"{comparison['replaced']} AW replacements); "
                      "includes environment atoms where present")
            else:
                _warn(f"Frame {frame_number}: can't do matched-atom comparison "
                      "without all three networks; using available group volumes")

        frame_dict[frame_number] = values
    return dict(sorted(frame_dict.items()))


def _time_ns(frame, ns_per_frame=1.0):
    return float(frame) * ns_per_frame


def _comparison_plot(frame_dict, metric, kind, ns_per_frame=1.0,
                     reference_frame=0, save_dir=None):
    """Plot AW-relative differences or each method's deviation from t=0."""
    title, _, stem = METRICS[metric]
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    frames = sorted(frame_dict)
    methods = ("pow", "prm") if kind == "aw_difference" else NETWORKS
    if kind == "initial_deviation" and reference_frame not in frame_dict:
        warnings.warn(f"Reference frame f{reference_frame} is missing for {title}; no baseline available")
    for method in methods:
        x, y = [], []
        baseline = _number(frame_dict.get(reference_frame, {}).get(metric, {}).get(method))
        for frame in frames:
            val = _number(frame_dict[frame].get(metric, {}).get(method))
            if kind == "aw_difference":
                ref = _number(frame_dict[frame].get(metric, {}).get("aw"))
                if not math.isfinite(ref) or not math.isfinite(val):
                    continue
                if ref == 0:
                    continue
                diff = 100.0 * (ref - val) / ref
            else:
                if not math.isfinite(baseline) or baseline == 0 or not math.isfinite(val):
                    continue
                diff = 100.0 * (val - baseline) / baseline
            x.append(_time_ns(frame, ns_per_frame))
            y.append(diff)
        if x:
            ax.plot(x, y, "-o", linewidth=2, markersize=5, label=NETWORK_LABELS[method],
                    color=NETWORK_COLORS[method])
    ax.axhline(0, color="gray", linewidth=0.8)
    ax.set_xlabel("Time (ns)", fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel(
        "Difference from AW (%)" if kind == "aw_difference"
        else "Deviation from initial frame (%)",
        fontsize=AXIS_LABEL_SIZE
    )
    if SHOW_TITLES:
        ax.set_title(f"{title}: " + ("difference from AW" if kind == "aw_difference" else "deviation from initial frame"),
                     fontsize=TITLE_FONT_SIZE)
    if frames:
        ax.set_xticks([_time_ns(f, ns_per_frame) for f in frames])
    ax.tick_params(axis="both", labelsize=TICK_LABEL_SIZE)
    # ax.grid(alpha=0.25)
    if ax.lines and any(line.get_label() in NETWORK_LABELS.values() for line in ax.lines):
        ax.legend(frameon=False)
    else:
        ax.text(0.5, 0.5, "No valid comparisons", transform=ax.transAxes, ha="center", va="center")

    # Apply individual plot axis limits
    metric_key = {
        "vol": "volume",
        "sa": "surface",
        "contacts": "contacts",
        "contact_count": "contacts",
    }[metric]

    kind_key = "aw_diff" if kind == "aw_difference" else "frame_dev"
    plot_key = f"{metric_key}_{kind_key}"

    limits = PLOT_LIMITS[plot_key]

    if limits["xlim"] is not None:
        ax.set_xlim(*limits["xlim"])

    if limits["ylim"] is not None:
        ax.set_ylim(*limits["ylim"])
    fig.tight_layout()
    paths = []
    if save_dir is not None:
        out = Path(save_dir)
        out.mkdir(parents=True, exist_ok=True)
        for ext in ("png", "svg"):
            path = out / f"{stem}_{kind}.{ext}"
            fig.savefig(path, dpi=200)
            paths.append(path)
    return fig, paths


def plot_frame_diffs(frame_dict, metric="vol", ns_per_frame=1.0, save_dir=None):
    return _comparison_plot(frame_dict, metric, "aw_difference", ns_per_frame, save_dir=save_dir)[0]


def plot_frame_deviations(frame_dict, metric="vol", ns_per_frame=1.0,
                          reference_frame=0, save_dir=None):
    return _comparison_plot(frame_dict, metric, "initial_deviation", ns_per_frame,
                            reference_frame, save_dir)[0]


def plot_frame_metrics(frame_dict, save_dir=None, show=True,
                       plot_contact_count=False, ns_per_frame=1.0,
                       reference_frame=0):
    metrics = ("vol", "sa", "contacts")  # Contacts now means contact count.
    paths = []
    for metric in metrics:
        for kind in ("aw_difference", "initial_deviation"):
            fig, saved = _comparison_plot(frame_dict, metric, kind, ns_per_frame,
                                          reference_frame, save_dir)
            paths.extend(saved)
            if not show:
                plt.close(fig)
    if show:
        plt.show()
    return paths


def diagnose_log(path):
    """Inspect one actual log, independent of the other two network logs."""
    path = Path(path)
    record = read_logs2(str(path), all_=False, balls=True, surfs=True)
    group = record.get("group data", {})
    atoms = record.get("atoms", pd.DataFrame())
    contact_area, contact_count = _surface_contacts(record, path)
    print(f"Log: {path.name}")
    print(f"  Group Volume:       {_number(group.get('Volume')):.6f} Å³")
    print(f"  Group Surface Area: {_number(group.get('Surface Area')):.6f} Å²")
    if not atoms.empty and "Volume" in atoms:
        print(f"  ALL atom volumes:   {atoms['Volume'].sum():.6f} Å³ "
              f"(includes any solvent/context atoms)")
        if "Residue" in atoms:
            print(f"  Atom residues:      {atoms['Residue'].value_counts().to_dict()}")
    print(f"  Surface contact area: {contact_area:.6f} Å²")
    print(f"  Contact patches with positive area: {contact_count}")
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--folder", help="Parent folder containing f0, f1, ... f10")
    parser.add_argument("--log", help="Diagnose a single AW/POW/PRM CSV log, then exit")
    parser.add_argument("--save-dir", help="Where to write SVG and PNG plots (default: <folder>/FrameAnalysis_plots)")
    parser.add_argument("--no-show", action="store_true", help="Save without opening plot windows")
    parser.add_argument("--plot-contact-count", action="store_true", help="Deprecated: contact count is now the default contact metric")
    parser.add_argument("--volume-mode", default="group", choices=("group", "legacy_atom_sum"),
                        help="Use Group Information (recommended) or old matched-atom volume comparison")
    parser.add_argument("--ns-per-frame", type=float, default=1.0, help="Time spacing between frames in ns (default: 1)")
    parser.add_argument("--reference-frame", type=int, default=0, help="Baseline frame for deviations (default: f0)")
    args = parser.parse_args(argv)

    if args.log:
        diagnose_log(args.log)
        return

    folder = args.folder
    if not folder:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        try:
            folder = filedialog.askdirectory(title="Select VorPy frames folder")
        finally:
            root.destroy()
        if not folder:
            print("No folder selected")
            return
    save_dir = args.save_dir
    if save_dir is None:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        try:
            save_dir = filedialog.askdirectory(title="Select VorPy frames folder")
        finally:
            root.destroy()
        if not save_dir:
            save_dir = str(Path(folder) / "FrameAnalysis_plots")

    data = get_frames_by_group(folder=folder, volume_mode=args.volume_mode)
    print(f"\nPer-frame measurements (volume mode: {args.volume_mode}):")
    for number, values in data.items():
        print(f"{number}: {values}")
    paths = plot_frame_metrics(data, save_dir=save_dir, show=not args.no_show,
                               plot_contact_count=args.plot_contact_count,
                               ns_per_frame=args.ns_per_frame, reference_frame=args.reference_frame)
    print("\nSaved plots:")
    for path in paths:
        print(f"  {path}")


if __name__ == "__main__":
    main()
