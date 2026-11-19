import os
import sys
import ast
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tkinter import Tk, filedialog


vorpy_root = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..")
)
sys.path.append(vorpy_root)


from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


def parse_balls(value):
    if isinstance(value, (list, tuple, np.ndarray)):
        return int(value[0]), int(value[1])

    if isinstance(value, str):
        cleaned = value.strip()

        try:
            parsed = ast.literal_eval(cleaned)
            return int(parsed[0]), int(parsed[1])
        except Exception:
            cleaned = cleaned.replace("[", "").replace("]", "").replace(",", " ")
            parts = cleaned.split()
            return int(parts[0]), int(parts[1])

    raise ValueError(f"Could not parse Balls value: {value}")


def normalize_pair(pair):
    a, b = pair
    return tuple(sorted([int(a), int(b)]))


def infer_element_from_atom_name(atom_name):
    """
    Infer chemical element from PDB-style atom name.
    Examples:
        CA -> C for protein alpha carbon
        C1' -> C
        OP1 -> O
        H5' -> H
        NA -> N unless explicitly sodium handling is needed
    """

    name = str(atom_name).strip().upper()

    if not name:
        return "UNK"

    # Remove common digits/symbols from atom naming
    cleaned = (
        name.replace("'", "")
        .replace("*", "")
        .replace('"', "")
        .replace("_", "")
        .strip()
    )

    # PDB atom names usually start with element-ish character after leading digits
    cleaned = cleaned.lstrip("0123456789")

    if not cleaned:
        return "UNK"

    # Common biomolecular elements
    first = cleaned[0]

    if first in {"C", "H", "N", "O", "P", "S"}:
        return first

    # Optional metal/ion handling
    two_letter = cleaned[:2]
    if two_letter in {"MG", "ZN", "FE", "CA", "NA", "CL", "MN", "CU", "K"}:
        return two_letter

    return first


def get_element_column(atoms_df):
    possible_cols = [
        "Element",
        "element",
        "Atom Element",
        "atom_element",
        "AtomType",
        "Atom Type",
        "Type",
    ]

    for col in possible_cols:
        if col in atoms_df.columns:
            return col

    if "Name" in atoms_df.columns:
        atoms_df["Inferred Element"] = atoms_df["Name"].apply(infer_element_from_atom_name)
        return "Inferred Element"

    raise ValueError(
        "Could not find an element/type column or Name column in atoms_df. "
        f"Available columns: {list(atoms_df.columns)}"
    )


def assign_energy_rank(pair_type):
    # crude but effective first-pass ranking

    low = {"C-C", "C-H", "H-H"}
    medium = {"C-O", "C-N", "H-O", "H-N"}
    high = {"O-P", "N-O", "H-P", "C-P", "O-O", "N-N"}

    if pair_type in low:
        return 1
    if pair_type in medium:
        return 2
    if pair_type in high:
        return 3

    return np.nan


def add_surface_pair_columns(surfs_df):
    df = surfs_df.copy()

    df[["Ball_1", "Ball_2"]] = df["Balls"].apply(
        lambda x: pd.Series(parse_balls(x))
    )

    df["Pair"] = df[["Ball_1", "Ball_2"]].apply(
        lambda row: normalize_pair((row["Ball_1"], row["Ball_2"])),
        axis=1
    )

    df["Abs Mean Curvature"] = df["Mean Curvature"].abs()

    return df



def add_atom_pair_types(atoms_df, surfs_df):
    df = add_surface_pair_columns(surfs_df)

    element_col = get_element_column(atoms_df)

    atom_type_map = {
        int(row["Index"]): str(row[element_col]).strip()
        for _, row in atoms_df.iterrows()
    }

    def get_pair_type(row):
        t1 = atom_type_map.get(int(row["Ball_1"]), "UNK")
        t2 = atom_type_map.get(int(row["Ball_2"]), "UNK")
        return "-".join(sorted([t1, t2]))

    df["AtomType_1"] = df["Ball_1"].map(atom_type_map).fillna("UNK")
    df["AtomType_2"] = df["Ball_2"].map(atom_type_map).fillna("UNK")
    df["PairType"] = df.apply(get_pair_type, axis=1)

    return df


def classify_pair_energy(pair_type):
    nonpolar = {"C", "H"}
    polar = {"O", "N", "S", "P"}

    parts = pair_type.split("-")

    if len(parts) != 2:
        return "Other"

    a, b = parts

    if "UNK" in {a, b}:
        return "External/Unknown"

    if a in nonpolar and b in nonpolar:
        return "Nonpolar"

    if (a in nonpolar and b in polar) or (a in polar and b in nonpolar):
        return "Mixed"

    if a in polar and b in polar:
        return "Polar"

    return "Other"

def add_literature_surface_energy_comparison(df):
    """
    Compares literature roughness-based surface energy scaling
    against our curvature-weighted surface-energy-like quantity.

    Literature rough surface model:
        gamma / gamma_0 = 1 + k * roughness_parameter

    Here:
        local roughness parameter ≈ |H| * sqrt(A)

    This is dimensionless:
        H has units 1/length
        sqrt(A) has units length
    """

    out_df = df.copy()

    out_df["Local Roughness Parameter"] = (
        out_df["Mean Curvature"].abs() * np.sqrt(out_df["Surface Area"])
    )

    # Wang et al. stochastic rough surface model
    out_df["Literature Surface Energy Norm"] = (
        1.0 + 1.04 * out_df["Local Roughness Parameter"]
    )

    # Your current geometry-based quantity
    out_df["Calculated Surface Energy Raw"] = (
        out_df["Mean Curvature"].abs()
    )

    median_val = out_df["Calculated Surface Energy Raw"].median()

    if median_val == 0:
        median_val = out_df["Calculated Surface Energy Raw"].replace(0, np.nan).median()

    out_df["Calculated Surface Energy Norm"] = (
        out_df["Calculated Surface Energy Raw"] / median_val
    )

    return out_df


def plot_literature_vs_calculated_surface_energy(df, output_dir):
    plot_df = add_literature_surface_energy_comparison(df)

    plot_df = plot_df[
        plot_df["EnergyClass"].isin(["Nonpolar", "Mixed", "Polar"])
    ].copy()

    fig, ax = plt.subplots(figsize=(9, 7))

    class_colors = {
        "Nonpolar": "tab:blue",
        "Mixed": "tab:orange",
        "Polar": "tab:red",
    }

    for cls, color in class_colors.items():
        cls_df = plot_df[plot_df["EnergyClass"] == cls]

        ax.scatter(
            cls_df["Literature Surface Energy Norm"],
            cls_df["Calculated Surface Energy Norm"],
            s=35,
            alpha=0.55,
            label=cls,
            edgecolors="black",
            linewidths=0.25
        )

    ax.set_xlabel("Literature Roughness-Based Surface Energy, γ/γ₀", fontsize=16)
    ax.set_ylabel("Calculated Curvature-Weighted Surface Energy, normalized", fontsize=16)

    ax.set_title(
        "Literature Surface-Energy Scaling vs Calculated AW Surface-Energy Proxy",
        fontsize=17
    )

    ax.legend(title="Interaction Class", fontsize=12, title_fontsize=13)

    ax.tick_params(axis="both", labelsize=13, width=1.5, length=6)

    for spine in ax.spines.values():
        spine.set_linewidth(1.5)

    plt.tight_layout()

    out_png = os.path.join(output_dir, "literature_vs_calculated_surface_energy.png")
    out_svg = os.path.join(output_dir, "literature_vs_calculated_surface_energy.svg")

    plt.savefig(out_png, dpi=300)
    plt.savefig(out_svg)
    plt.close()

    plot_df.to_csv(
        os.path.join(output_dir, "literature_vs_calculated_surface_energy_data.csv"),
        index=False
    )

    print(f"Saved literature comparison plot: {out_png}")


def parse_neighbors(value):
    if isinstance(value, (list, tuple, set, np.ndarray)):
        return [int(v) for v in value]

    if pd.isna(value):
        return []

    if isinstance(value, str):
        cleaned = value.strip()

        try:
            parsed = ast.literal_eval(cleaned)
            if isinstance(parsed, (list, tuple, set)):
                return [int(v) for v in parsed]
        except Exception:
            pass

        cleaned = cleaned.replace("[", "").replace("]", "").replace(",", " ")
        return [int(v) for v in cleaned.split() if v.strip().lstrip("-").isdigit()]

    try:
        return [int(value)]
    except Exception:
        return []


def plot_energy_rank_vs_curvature(df, output_dir):
    plot_df = df.dropna(subset=["EnergyRank"]).copy()

    fig, ax = plt.subplots(figsize=(8,6))

    ax.scatter(
        plot_df["EnergyRank"],
        plot_df["Mean Curvature"].abs(),
        alpha=0.5,
        s=30,
        edgecolors="black",
        linewidths=0.2
    )

    # add mean trend
    means = plot_df.groupby("EnergyRank")["Mean Curvature"].apply(lambda x: np.mean(np.abs(x)))
    ax.plot(means.index, means.values, color="red", linewidth=2)

    ax.set_xlabel("Interaction Energy Rank", fontsize=14)
    ax.set_ylabel("|Mean Curvature|", fontsize=14)
    ax.set_title("Curvature vs Interaction Energy Rank", fontsize=16)

    plt.tight_layout()

    rank_summary = (
        plot_df.groupby("EnergyRank")
        .agg(
            Count=("Abs Mean Curvature", "count"),
            MeanAbsCurvature=("Abs Mean Curvature", "mean"),
            MedianAbsCurvature=("Abs Mean Curvature", "median"),
            StdAbsCurvature=("Abs Mean Curvature", "std"),
            MeanSurfaceArea=("Surface Area", "mean"),
        )
        .reset_index()
    )

    rank_summary.to_csv(
        os.path.join(output_dir, "energy_rank_curvature_summary.csv"),
        index=False
    )

    print("\nEnergy rank curvature summary:")
    print(rank_summary)
    plt.savefig(os.path.join(output_dir, "energy_rank_vs_curvature.png"), dpi=300)
    plt.close()


def get_bonded_pairs(atoms_df, bond_col="Neighbors"):
    if bond_col not in atoms_df.columns:
        raise ValueError(
            f"Could not find bond column '{bond_col}'. "
            f"Available columns: {list(atoms_df.columns)}"
        )

    bonded_pairs = set()

    for _, row in atoms_df.iterrows():
        atom_index = int(row["Index"])
        neighbors = parse_neighbors(row[bond_col])

        for neighbor in neighbors:
            bonded_pairs.add(normalize_pair((atom_index, neighbor)))

    return bonded_pairs



def add_bonded_flag(atoms_df, surfs_df, bond_col="Neighbors"):
    df = add_atom_pair_types(atoms_df, surfs_df)

    bonded_pairs = get_bonded_pairs(atoms_df, bond_col=bond_col)

    df["IsBonded"] = df["Pair"].isin(bonded_pairs)

    return df



def plot_curvature_by_energy_class(df, save_path=None, title=None):
    plot_df = df[df["EnergyClass"] != "Other"].copy()

    preferred_order = ["Nonpolar", "Mixed", "Polar"]

    order = [
        cls for cls in preferred_order
        if cls in plot_df["EnergyClass"].unique()
        and plot_df.loc[plot_df["EnergyClass"] == cls, "Mean Curvature"].dropna().size > 0
    ]

    if len(order) == 0:
        print("No valid energy classes found for plotting.")
        print(df["EnergyClass"].value_counts(dropna=False))
        return

    data = [
        plot_df.loc[plot_df["EnergyClass"] == cls, "Mean Curvature"].dropna().to_numpy()
        for cls in order
    ]

    fig, ax = plt.subplots(figsize=(10, 7))

    ax.violinplot(
        data,
        showmeans=True,
        showmedians=True,
        showextrema=True
    )

    ax.axhline(0, linestyle="--", linewidth=1.5)

    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels(order, fontsize=14)

    ax.set_xlabel("Atom-Pair Interaction Class", fontsize=16)
    ax.set_ylabel("Mean Curvature", fontsize=16)

    if title is None:
        title = "AW Interface Curvature by Atom-Pair Interaction Class"

    ax.set_title(title, fontsize=18)

    ax.tick_params(axis="both", labelsize=13, width=1.5, length=6)

    for spine in ax.spines.values():
        spine.set_linewidth(1.5)

    plt.tight_layout()

    if save_path is not None:
        plt.savefig(save_path, dpi=300)
        plt.close()
    else:
        plt.show()



def plot_pairtype_curvature(df, save_path=None, title=None, min_count=25):
    counts = df["PairType"].value_counts()
    keep_pairtypes = counts[counts >= min_count].index.tolist()

    plot_df = df[df["PairType"].isin(keep_pairtypes)].copy()

    order = (
        plot_df.groupby("PairType")["Mean Curvature"]
        .median()
        .sort_values()
        .index
        .tolist()
    )

    data = [
        plot_df.loc[plot_df["PairType"] == pair_type, "Mean Curvature"].dropna()
        for pair_type in order
    ]

    fig, ax = plt.subplots(figsize=(max(10, 0.65 * len(order)), 7))

    ax.violinplot(
        data,
        showmeans=True,
        showmedians=True,
        showextrema=True
    )

    ax.axhline(0, linestyle="--", linewidth=1.5)

    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels(order, rotation=45, ha="right", fontsize=12)

    ax.set_xlabel("Atom Pair Type", fontsize=16)
    ax.set_ylabel("Mean Curvature", fontsize=16)

    if title is None:
        title = "AW Interface Curvature by Atom Pair Type"

    ax.set_title(title, fontsize=18)

    ax.tick_params(axis="both", labelsize=13, width=1.5, length=6)

    for spine in ax.spines.values():
        spine.set_linewidth(1.5)

    plt.tight_layout()

    if save_path is not None:
        plt.savefig(save_path, dpi=300)
        plt.close()
    else:
        plt.show()



def write_summary(df, save_path):
    summary = (
        df.groupby(["EnergyClass", "PairType"])
        .agg(
            Count=("Mean Curvature", "count"),
            MeanCurvature=("Mean Curvature", "mean"),
            MedianCurvature=("Mean Curvature", "median"),
            StdCurvature=("Mean Curvature", "std"),
            MeanAbsCurvature=("Abs Mean Curvature", "mean"),
            MeanSurfaceArea=("Surface Area", "mean"),
        )
        .reset_index()
        .sort_values(["EnergyClass", "PairType"])
    )

    summary.to_csv(save_path, index=False)

    return summary


def plot_pairtype_curvature_bar(df, output_dir, min_count=50):
    pair_summary = (
        df.dropna(subset=["Mean Curvature"])
        .groupby("PairType")
        .agg(
            Count=("Mean Curvature", "count"),
            MeanAbsCurvature=("Mean Curvature", lambda x: np.mean(np.abs(x))),
            MedianAbsCurvature=("Mean Curvature", lambda x: np.median(np.abs(x))),
            StdAbsCurvature=("Mean Curvature", lambda x: np.std(np.abs(x))),
            MeanSurfaceArea=("Surface Area", "mean"),
        )
        .query("Count >= @min_count")
        .sort_values("MedianAbsCurvature")
    )

    out_csv = os.path.join(output_dir, "pairtype_curvature_summary.csv")
    pair_summary.to_csv(out_csv)

    fig, ax = plt.subplots(figsize=(14, 6))

    x = np.arange(len(pair_summary))

    ax.bar(
        x,
        pair_summary["MedianAbsCurvature"]
    )

    ax.set_xticks(x)
    ax.set_xticklabels(pair_summary.index, rotation=45, ha="right")

    ax.set_xlabel("Atom-Pair Type", fontsize=14)
    ax.set_ylabel("Median |Mean Curvature|", fontsize=14)
    ax.set_title("Median AW Curvature by Atom-Pair Type", fontsize=16)

    ax.tick_params(axis="both", labelsize=12, width=1.4, length=5)

    for spine in ax.spines.values():
        spine.set_linewidth(1.4)

    plt.tight_layout()

    out_png = os.path.join(output_dir, "pairtype_curvature_bar.png")
    out_svg = os.path.join(output_dir, "pairtype_curvature_bar.svg")

    plt.savefig(out_png, dpi=300)
    plt.savefig(out_svg)
    plt.close()

    print(f"Saved pair-type curvature bar plot: {out_png}")
    print(f"Saved pair-type curvature summary: {out_csv}")


def run_surface_energy_curvature_analysis():
    Tk().withdraw()

    folder = filedialog.askdirectory(
        title="Choose model folder containing aw/aw_logs.csv"
    )

    if not folder:
        print("No folder selected.")
        return

    aw_logs_path = os.path.join(folder, "aw", "aw_logs.csv")
    output_dir = os.path.join(folder, "aw", "surface_energy_curvature")

    os.makedirs(output_dir, exist_ok=True)

    logs = read_logs2(aw_logs_path)

    atoms_df = logs["atoms"]
    surfs_df = logs["surfs"]

    full_df = add_bonded_flag(
        atoms_df=atoms_df,
        surfs_df=surfs_df,
        bond_col="Neighbors"
    )

    full_df["EnergyClass"] = full_df["PairType"].apply(classify_pair_energy)
    full_df["EnergyRank"] = full_df["PairType"].apply(assign_energy_rank)

    full_df["AbsCurvature"] = full_df["Mean Curvature"].abs()

    full_df["AreaWeightedCurvature"] = (
            full_df["Surface Area"] * full_df["AbsCurvature"]
    )

    bonded_df = full_df[full_df["IsBonded"]].copy()
    nonbonded_df = full_df[~full_df["IsBonded"]].copy()

    full_csv = os.path.join(output_dir, "all_surfaces_curvature_energy_classes.csv")
    bonded_csv = os.path.join(output_dir, "bonded_surfaces_curvature_energy_classes.csv")
    nonbonded_csv = os.path.join(output_dir, "nonbonded_surfaces_curvature_energy_classes.csv")
    summary_csv = os.path.join(output_dir, "curvature_by_energy_class_summary.csv")

    full_df.to_csv(full_csv, index=False)
    bonded_df.to_csv(bonded_csv, index=False)
    nonbonded_df.to_csv(nonbonded_csv, index=False)

    write_summary(full_df, summary_csv)

    energy_surface_summary = (
        full_df.dropna(subset=["EnergyRank"])
        .groupby("EnergyRank")
        .agg(
            TotalArea=("Surface Area", "sum"),
            TotalCurvature=("AreaWeightedCurvature", "sum"),
            MeanCurvature=("AbsCurvature", "mean"),
            AreaWeightedMeanCurvature=("AreaWeightedCurvature", "sum"),
        )
    )

    energy_surface_summary["SurfaceWeightedCurvature"] = (
            energy_surface_summary["TotalCurvature"] /
            energy_surface_summary["TotalArea"]
    )

    plot_pairtype_curvature_bar(
        df=full_df,
        output_dir=output_dir,
        min_count=50
    )

    print("\nAll surfaces energy classes:")
    print(full_df["EnergyClass"].value_counts(dropna=False))

    print("\nBonded surfaces energy classes:")
    print(bonded_df["EnergyClass"].value_counts(dropna=False))

    print("\nNonbonded surfaces energy classes:")
    print(nonbonded_df["EnergyClass"].value_counts(dropna=False))

    print(full_df["PairType"].value_counts().head(40))
    print(full_df.loc[full_df["EnergyClass"] == "Other", "PairType"].value_counts().head(40))

    plot_curvature_by_energy_class(
        full_df,
        save_path=os.path.join(output_dir, "all_surfaces_curvature_by_energy_class.png"),
        title="All AW Interfaces: Curvature by Atom-Pair Interaction Class"
    )

    plot_curvature_by_energy_class(
        bonded_df,
        save_path=os.path.join(output_dir, "bonded_surfaces_curvature_by_energy_class.png"),
        title="Bonded AW Interfaces: Curvature by Atom-Pair Interaction Class"
    )

    plot_curvature_by_energy_class(
        nonbonded_df,
        save_path=os.path.join(output_dir, "nonbonded_surfaces_curvature_by_energy_class.png"),
        title="Nonbonded AW Interfaces: Curvature by Atom-Pair Interaction Class"
    )

    plot_pairtype_curvature(
        full_df,
        save_path=os.path.join(output_dir, "all_surfaces_curvature_by_pairtype.png"),
        title="All AW Interfaces: Curvature by Atom Pair Type",
        min_count=25
    )

    plot_literature_vs_calculated_surface_energy(
        df=full_df,
        output_dir=output_dir
    )

    plot_energy_rank_vs_curvature(
        df=full_df,
        output_dir=output_dir
    )

    print("\nSaved outputs:")
    print(full_csv)
    print(bonded_csv)
    print(nonbonded_csv)
    print(summary_csv)
    print("\nCounts by energy class:")
    print(full_df["EnergyClass"].value_counts())
    print("\nBonded surfaces:", len(bonded_df))
    print("Nonbonded surfaces:", len(nonbonded_df))
    print("\nCounts by energy rank:")
    print(full_df["EnergyRank"].value_counts(dropna=False).sort_index())


if __name__ == "__main__":
    run_surface_energy_curvature_analysis()
