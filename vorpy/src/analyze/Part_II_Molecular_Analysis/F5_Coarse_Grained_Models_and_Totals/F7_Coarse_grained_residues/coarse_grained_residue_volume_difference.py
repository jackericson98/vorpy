import os
import sys
import tkinter as tk
from tkinter import filedialog

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2



"""
Figure 7A: Residue volume "fingerprint"
- X axis: 25 residues
- For each residue: 5 violins (CG schemes)
- Y axis: coarse-grained residue volume (Å^3)
"""



AA_ORDER = [
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
]


NA_ORDER = ["DA", "DC", "DG", "DT", "DU"]
RES_ORDER = AA_ORDER + NA_ORDER


CG_FOLDER_TO_LABEL = {
    "1_Atom": "Atom",
    "2_Encap": "ENC",
    "3_Encap_SR": "ENC-SR",
    "4_AD": "AD",
    "5_AD_SR": "AD-SR",
    "6_AD_MW": "AD-MW",
    "7_AD_MW_SR": "AD-MW-SR",
}


DEFAULT_CG_LABEL_ORDER = ["Atom", "AD", "AD-SR", "AD-MW", "AD-MW-SR", "ENC-SR", "ENC"]

CG_LABEL_TO_FOLDER = {v: k for k, v in CG_FOLDER_TO_LABEL.items()}


VOR_SCHEMES = ("aw", "pow", "prm")


def normalize_residue_label(res: str) -> str:
    r = str(res).strip().upper()

    # Map RNA single-letter bases to deoxy-style
    if r in {"A", "C", "G", "U", "T"}:
        if r == "A":
            return "DA"
        if r == "C":
            return "DC"
        if r == "G":
            return "DG"
        if r == "T":
            return "DT"
        if r == "U":
            return "DU"

    # Already deoxy-style
    if r in {"DA", "DC", "DG", "DT", "DU"}:
        return r

    # Otherwise (amino acids etc.)
    return r


def compute_aa_rank_groups(df: pd.DataFrame,
                           vor_scheme: str,
                           atom_label: str = "Atom",
                           aa_order: list[str] | None = None) -> dict[str, list[str]]:
    """
    Returns AA groups based on Atom median volumes for the chosen vor_scheme:
      - largest_7: top 7
      - middle_6: next 6
      - smallest_7: bottom 7
    """
    if aa_order is None:
        aa_order = AA_ORDER

    sub = df[
        (df["vor_scheme"] == vor_scheme) &
        (df["cg_scheme"] == atom_label) &
        (df["Residue"].isin(aa_order))
    ].copy()

    if len(sub) == 0:
        raise RuntimeError("No Atom data available to compute amino-acid rank groups.")

    med = (
        sub.groupby("Residue")["residue_volume"]
        .median()
        .reindex(aa_order)
        .dropna()
        .sort_values(ascending=False)
    )

    aa_sorted = med.index.tolist()
    if len(aa_sorted) != 20:
        print(f"Warning: expected 20 amino acids, found {len(aa_sorted)} in Atom medians.")

    largest_7 = aa_sorted[:7]
    middle_6 = aa_sorted[7:13]
    smallest_7 = aa_sorted[13:20]

    return {
        "largest_7": largest_7,
        "middle_6": middle_6,
        "smallest_7": smallest_7,
    }


def choose_root_folder() -> str:
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes("-topmost", 1)

    folder = filedialog.askdirectory(title="Select root folder containing model folders")
    if not folder:
        raise RuntimeError("No folder selected.")
    return folder


def iter_logs_files(root_folder: str,
                    model_folders: set[str],
                    cg_labels_to_use: list[str],
                    vor_schemes: tuple[str, ...] = ("aw", "pow", "prm")) -> list[dict]:
    """
    Traverse:
        root / model / cg_scheme / {aw,pow,prm} / *_logs.csv

    cg_labels_to_use are the human-readable scheme labels ("Atom", "AD", ...),
    mapped to folder names via CG_LABEL_TO_FOLDER.
    """
    cg_folders = []
    for lab in cg_labels_to_use:
        folder = CG_LABEL_TO_FOLDER.get(lab, None)
        print(folder, lab)
        if folder is not None:
            cg_folders.append((folder, lab))

    found = []

    for model in sorted(model_folders):
        model_dir = os.path.join(root_folder, model)
        if not os.path.isdir(model_dir):
            continue

        for cg_folder, cg_label in cg_folders:
            cg_dir = os.path.join(model_dir, cg_folder)
            if not os.path.isdir(cg_dir):
                continue

            for vor in vor_schemes:
                vor_dir = os.path.join(cg_dir, vor)
                if not os.path.isdir(vor_dir):
                    continue

                canonical = os.path.join(vor_dir, f"{vor}_logs.csv")
                if os.path.isfile(canonical):
                    found.append(
                        {
                            "model": model,
                            "cg_folder": cg_folder,
                            "cg_label": cg_label,
                            "vor_scheme": vor,
                            "logs_path": canonical,
                        }
                    )
                    continue

                for fname in os.listdir(vor_dir):
                    if fname.endswith("_logs.csv") and fname.lower().startswith(vor):
                        found.append(
                            {
                                "model": model,
                                "cg_folder": cg_folder,
                                "cg_label": cg_label,
                                "vor_scheme": vor,
                                "logs_path": os.path.join(vor_dir, fname),
                            }
                        )

    return found


def build_residue_volume_df(root_folder: str,
                            vor_scheme: str,
                            cg_labels_to_use: list[str],
                            model_folders: set[str] | None = None,
                            residue_allowlist: list[str] | None = None,
                            filter_solvent: bool = True) -> pd.DataFrame:
    """
    Build a long dataframe where each row is ONE residue instance (chain, resseq)
    and the value is the sum of ball volumes for that residue instance.

    For split-residue schemes, this sums multiple beads back to the residue instance.
    """
    if residue_allowlist is None:
        residue_allowlist = RES_ORDER

    files = iter_logs_files(
        root_folder=root_folder,
        model_folders=model_folders,
        cg_labels_to_use=cg_labels_to_use,
    )
    files = [f for f in files if f["vor_scheme"] == vor_scheme]
    files = [f for f in files if f["vor_scheme"] == vor_scheme]

    rows = []
    n_files = len(files)

    if n_files == 0:
        raise RuntimeError(f"No logs files found for vor_scheme='{vor_scheme}' under: {root_folder}")

    for i, meta in enumerate(files, start=1):
        logs_path = meta["logs_path"]
        print(f"[{i}/{n_files}] Reading: {logs_path}")

        logs_obj = read_logs2(
            logs_path,
            return_dict=False,
            no_sol=False,  # avoid KeyError in read_logs2 solvent filter
            all_=False,
            balls=True,
            surfs=False,
            edges=False,
            verts=False,
        )
        atoms_df = logs_obj.get("atoms", None)

        if atoms_df is None or len(atoms_df) == 0:
            continue

        atoms_df["Residue"] = atoms_df["Residue"].apply(normalize_residue_label)

        if filter_solvent and atoms_df is not None and "Name" in atoms_df.columns:
            atom_name = atoms_df["Name"].astype(str).str.strip().str.lower()
            drop = {"hw1", "hw2", "ow", "h02", "h01", "na", "cl", "mg", "k"}
            atoms_df = atoms_df.loc[~atom_name.isin(drop)].copy()

        needed = {"Residue", "Chain", "Residue Sequence", "Volume"}
        missing = needed.difference(set(atoms_df.columns))
        if missing:
            print(f"  !! Missing columns {sorted(missing)} in {logs_path}")
            continue

        # drop solvent/ions if desired
        if "Name" in atoms_df.columns:
            atom_name = atoms_df["Name"].astype(str).str.strip().str.lower()
            drop = {"hw1", "hw2", "ow", "h02", "h01", "na", "cl", "mg", "k"}
            atoms_df = atoms_df.loc[~atom_name.isin(drop)].copy()

        g = (
            atoms_df
            .groupby(["Residue", "Chain", "Residue Sequence"], as_index=False)["Volume"]
            .sum()
            .rename(columns={"Volume": "residue_volume"})
        )
        print(atoms_df.columns.tolist())

        if len(g) == 0:
            continue

        g["model"] = meta["model"]
        g["cg_scheme"] = meta["cg_label"]
        g["vor_scheme"] = meta["vor_scheme"]
        g["source_file"] = os.path.basename(logs_path)

        rows.append(g)

    if not rows:
        raise RuntimeError("No residue volumes were extracted (all reads failed or filtered out).")

    df = pd.concat(rows, ignore_index=True)

    # Optional: enforce minimum sample counts later at plotting time; keep everything here.
    return df


def _quantiles(arr: np.ndarray) -> tuple[float, float, float]:
    q1 = float(np.quantile(arr, 0.25))
    q2 = float(np.quantile(arr, 0.50))
    q3 = float(np.quantile(arr, 0.75))
    return q1, q2, q3


def plot_figure_7a_grid(df: pd.DataFrame,
                        vor_scheme: str,
                        cg_order: list[str],
                        save: str | None = None,
                        show: bool = True,
                        n_min: int = 5,
                        fix_y_range: bool = True) -> None:

    groups = compute_aa_rank_groups(df=df, vor_scheme=vor_scheme, atom_label="Atom", aa_order=AA_ORDER)

    # -------- Optional global y-range --------
    if fix_y_range:
        all_vals = df["residue_volume"].to_numpy(dtype=float)
        all_vals = all_vals[np.isfinite(all_vals)]

        y_min = float(np.min(all_vals))
        y_max = float(np.max(all_vals))

        padding = 0.05 * (y_max - y_min) if y_max > y_min else 1.0
        y_limits = (y_min - padding, y_max + padding)
    else:
        y_limits = None
    # -----------------------------------------

    fig, axes = plt.subplots(
        nrows=2,
        ncols=2,
        figsize=(20, 10),
        sharey=False,  # we enforce equal limits manually
    )

    # Requested layout:
    # TL = nucleics
    # TR = largest 7 AAs
    # BL = middle 6 AAs
    # BR = smallest 7 AAs
    ax_tl = axes[0, 0]
    ax_tr = axes[0, 1]
    ax_bl = axes[1, 0]
    ax_br = axes[1, 1]

    TITLE_SIZE = 25
    LABEL_SIZE = 20
    TICK_SIZE = 18

    # Top-left: nucleics
    handles = plot_figure_7a(
        df=df,
        ax=ax_tl,
        title="Nucleic acids",
        residue_order=NA_ORDER,
        cg_order=cg_order,
        n_min=n_min,
        show_legend=False,
    )

    # Top-right: largest 7 amino acids
    plot_figure_7a(
        df=df,
        ax=ax_tr,
        title="Largest 7 amino acids",
        residue_order=groups["largest_7"],
        cg_order=cg_order,
        n_min=n_min,
        show_legend=False,
        ylabel="",
    )

    # Bottom-left: middle 6 amino acids
    plot_figure_7a(
        df=df,
        ax=ax_bl,
        title="Middle 6 amino acids",
        residue_order=groups["middle_6"],
        cg_order=cg_order,
        n_min=n_min,
        show_legend=False,
    )

    # Bottom-right: smallest 7 amino acids
    plot_figure_7a(
        df=df,
        ax=ax_br,
        title="Smallest 7 amino acids",
        residue_order=groups["smallest_7"],
        cg_order=cg_order,
        n_min=n_min,
        show_legend=False,
        ylabel="",
    )

    # Apply equal y-limits + font styling to all axes
    for ax in [ax_tl, ax_tr, ax_bl, ax_br]:
        ax.set_ylim(y_limits)

        ax.title.set_fontsize(TITLE_SIZE)
        ax.yaxis.label.set_fontsize(LABEL_SIZE)

        ax.tick_params(axis="both", labelsize=TICK_SIZE, width=2.5, length=8)

        for spine in ax.spines.values():
            spine.set_linewidth(2.5)

    fig.suptitle(
        f"Figure 7A | Coarse-grained residue volumes ({vor_scheme.upper()})",
        fontsize=22
    )

    # One shared legend for the full figure
    if handles:
        fig.legend(
            handles=handles,
            labels=[h.get_label() for h in handles],
            title="Coarse-graining scheme",
            fontsize=18,
            title_fontsize=20,
            frameon=False,
            loc="upper left",
            bbox_to_anchor=(0.01, 0.98),
        )

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    if save is not None:
        fig.savefig(save, dpi=300, bbox_inches="tight")

    if show:
        plt.show()

    plt.close(fig)


def plot_figure_7a(df: pd.DataFrame,
                   ax: plt.Axes,
                   title: str,
                   residue_order: list[str],
                   cg_order: list[str],
                   n_min: int = 5,
                   y_range: tuple[float, float] | None = None,
                   show_legend: bool = True,
                   ylabel: str = "Residue volume (Å³)") -> list[Patch]:
    # Stable CG -> color mapping
    cycle_colors = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if len(cycle_colors) == 0:
        cycle_colors = ["C0", "C1", "C2", "C3", "C4", "C5", "C6"]

    cg_to_color = {}
    for i, cg in enumerate(cg_order):
        cg_to_color[cg] = cycle_colors[i % len(cycle_colors)]

    offsets = np.linspace(-0.36, 0.36, len(cg_order))
    width = min(0.14, 0.75 / max(1, len(cg_order)))
    xs = np.arange(len(residue_order), dtype=float)

    ax.set_title(title, fontsize=14, pad=10)
    ax.set_xticks(xs)
    ax.set_xticklabels(residue_order, rotation=90, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=12)

    ax.grid(axis="y", linewidth=0.8, alpha=0.35)
    ax.tick_params(axis="both", width=2, length=7)
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)

    if y_range is not None:
        ax.set_ylim(y_range)

    legend_handles: list[Patch] = []

    for k, cg in enumerate(cg_order):
        color = cg_to_color[cg]
        pos = xs + offsets[k]

        data_per_res = []
        draw_pos = []

        for i, res in enumerate(residue_order):
            vals = df[(df["cg_scheme"] == cg) & (df["Residue"] == res)]["residue_volume"].to_numpy(dtype=float)
            vals = vals[np.isfinite(vals)]

            if vals.size < n_min:
                continue

            data_per_res.append(vals)
            draw_pos.append(float(pos[i]))

        if not data_per_res:
            continue

        vp = ax.violinplot(
            data_per_res,
            positions=draw_pos,
            widths=width,
            showmeans=False,
            showmedians=False,
            showextrema=False,
        )

        for body in vp["bodies"]:
            body.set_facecolor(color)
            body.set_edgecolor(color)
            body.set_alpha(0.55)
            body.set_linewidth(0.8)

        # Median + IQR overlays
        for vals, x in zip(data_per_res, draw_pos):
            q1 = float(np.quantile(vals, 0.25))
            q2 = float(np.quantile(vals, 0.50))
            q3 = float(np.quantile(vals, 0.75))

            ax.plot([x, x], [q1, q3], linewidth=2.0, color=color)
            ax.plot([x - width * 0.35, x + width * 0.35], [q2, q2], linewidth=2.5, color=color)

        legend_handles.append(Patch(facecolor=color, edgecolor=color, alpha=0.55, label=cg))

    if show_legend and legend_handles:
        ax.legend(handles=legend_handles, title="Coarse-graining scheme", frameon=False, fontsize=10, title_fontsize=11)

    return legend_handles


def main() -> None:
    root_folder = choose_root_folder()

    # -------- toggles --------
    vor_scheme = "pow"  # "aw" | "pow" | "prm"

    cg_labels_to_plot = [
        "Atom",
        "AD",
        "AD-MW",
        # "AD-SR",
        # "AD-MW-SR",
        "ENC",
        # "ENC-SR",
    ]

    include_models = {"F_BDNA", "G_Hammerhead", "I_T4LP", "K_NCP", "L_BSA"}
    filter_solvent = True
    # -------------------------

    df = build_residue_volume_df(
        root_folder=root_folder,
        vor_scheme=vor_scheme,
        cg_labels_to_use=cg_labels_to_plot,
        model_folders=include_models,
        residue_allowlist=RES_ORDER,
        filter_solvent=filter_solvent,
    )

    plot_figure_7a_grid(
        df=df,
        vor_scheme=vor_scheme,
        cg_order=cg_labels_to_plot,
        save=None,
        show=True,
        n_min=5,
        fix_y_range=False,
    )


if __name__ == "__main__":

    main()
