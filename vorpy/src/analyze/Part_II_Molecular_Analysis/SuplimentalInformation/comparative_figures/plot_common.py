from __future__ import annotations

from pathlib import Path
from tkinter import Tk, filedialog

import matplotlib.pyplot as plt
import pandas as pd


REQUIRED_STATIC_COLUMNS = {
    "system",
    "kind",
    "tier",
    "context",
    "method",
    "radius_model",
    "interface_atom_count",
    "interface_element_count",
    "candidate_elements",
    "eligible_elements",
    "selected_elements",
    "excluded_elements",
    "connected_components",
    "interface_area_A2",
    "signed_raw_edge_A_rad",
    "unsigned_raw_edge_A_rad",
    "positive_edge_A_rad",
    "negative_edge_A_rad",
    "smooth_H_A",
    "conventional_edge_H_A",
    "combined_H_A",
    "H_per_area_inv_A",
    "unsigned_edge_per_area_rad_per_A",
    "cancellation_fraction",
    "unresolved_quadratures",
}


METHOD_LABELS = {
    "full_power": "Full Power",
    "full_power_aw_radii": "Full Power\n(AW radii)",
    "alpha0_power_before_M": r"$\alpha=0$ Power",
    "cazals_alpha0_power": r"$\alpha=0$ + M=5",
    "aw": "AW",
}


SYSTEM_ORDER = ["2KAI", "1AK4", "1BRC"]


def select_csv(
    title: str = "Select static_method_comparison.csv",
) -> Path:
    """Open a Tkinter file picker and return one selected CSV path."""
    root = Tk()
    root.withdraw()
    root.attributes("-topmost", True)

    filename = filedialog.askopenfilename(
        title=title,
        filetypes=[
            ("CSV files", "*.csv"),
            ("All files", "*.*"),
        ],
    )

    root.destroy()

    if not filename:
        raise SystemExit("No file selected.")

    return Path(filename)


def select_output_directory(
    initial_dir: Path | None = None,
    title: str = "Select output directory for figures",
) -> Path:
    """Open a Tkinter directory picker."""
    root = Tk()
    root.withdraw()
    root.attributes("-topmost", True)

    directory = filedialog.askdirectory(
        title=title,
        initialdir=str(initial_dir) if initial_dir else None,
    )

    root.destroy()

    if not directory:
        raise SystemExit("No output directory selected.")

    path = Path(directory)
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_static_comparison(path: Path | None = None) -> tuple[pd.DataFrame, Path]:
    """
    Load and validate static_method_comparison.csv.

    Returns
    -------
    df
        Validated DataFrame.
    path
        Path actually loaded.
    """
    if path is None:
        path = select_csv()

    df = pd.read_csv(path)

    missing = REQUIRED_STATIC_COLUMNS - set(df.columns)
    if missing:
        missing_text = "\n  - ".join(sorted(missing))
        raise ValueError(
            "Selected CSV does not contain the required comparative-study "
            f"columns:\n  - {missing_text}"
        )

    if df.empty:
        raise ValueError("Selected CSV contains no rows.")

    return df, path


def static_protein_only(df: pd.DataFrame) -> pd.DataFrame:
    """Return static H-tier protein-only observations."""
    result = df.copy()

    if "kind" in result.columns:
        result = result[result["kind"] == "static"]

    if "tier" in result.columns:
        result = result[result["tier"] == "H"]

    if "context" in result.columns:
        result = result[result["context"] == "protein_only"]

    return result.copy()


def ordered_systems(df: pd.DataFrame) -> list[str]:
    """Use preferred publication ordering, followed by any extra systems."""
    available = list(dict.fromkeys(df["system"].astype(str)))

    ordered = [system for system in SYSTEM_ORDER if system in available]
    ordered.extend(system for system in available if system not in ordered)

    return ordered


def require_methods(
    df: pd.DataFrame,
    methods: list[str],
    systems: list[str] | None = None,
) -> None:
    """Fail loudly when a required system/method combination is missing."""
    if systems is None:
        systems = ordered_systems(df)

    missing: list[str] = []

    for system in systems:
        system_methods = set(
            df.loc[df["system"] == system, "method"].astype(str)
        )
        for method in methods:
            if method not in system_methods:
                missing.append(f"{system}: {method}")

    if missing:
        raise ValueError(
            "Required comparative-study rows are missing:\n  - "
            + "\n  - ".join(missing)
        )


def method_label(method: str) -> str:
    return METHOD_LABELS.get(method, method)


def configure_matplotlib() -> None:
    """Shared publication-oriented Matplotlib defaults."""
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 13,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
            "figure.titlesize": 15,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.bbox": "tight",
        }
    )


def save_figure(
    fig,
    output_dir: Path,
    stem: str,
    dpi: int = 600,
) -> list[Path]:
    """
    Save one figure as PNG, PDF, and SVG.

    Returns paths written.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    outputs = [
        output_dir / f"{stem}.png",
        output_dir / f"{stem}.pdf",
        output_dir / f"{stem}.svg",
    ]

    fig.savefig(outputs[0], dpi=dpi)
    fig.savefig(outputs[1])
    fig.savefig(outputs[2])

    return outputs