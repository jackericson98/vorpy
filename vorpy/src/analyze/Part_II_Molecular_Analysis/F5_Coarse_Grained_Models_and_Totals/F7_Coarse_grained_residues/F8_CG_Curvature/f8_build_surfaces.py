import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from tkinter import Tk
from tkinter import filedialog

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


# ----------------------------- residue classification -----------------------------

PROTEIN_RES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

DNA_RES = {
    "DA", "DT", "DG", "DC", "DI",
    "ADE", "THY", "GUA", "CYT",
}

RNA_RES = {
    "A", "U", "G", "C", "I",
    "RA", "RU", "RG", "RC", "RI",
}


def classify_residue(res_name: str) -> str:
    res = str(res_name).strip().upper()

    if res in PROTEIN_RES:
        return "protein"

    if res in DNA_RES:
        return "dna"

    if res in RNA_RES:
        return "dna"

    return "other"


# ----------------------------- directory walking -----------------------------

def iter_log_files(main_systems_dir: str):
    """
    Directory convention:
      Main_Systems/
        K_NCP/
          3_Encap_SR/
            aw/aw_logs.csv
            pow/pow_logs.csv
            prm/prm_logs.csv

    Yields: (model, cg_scheme, partition, log_path)
    """
    root = Path(main_systems_dir)

    for model_dir in sorted(root.iterdir()):
        if not model_dir.is_dir():
            continue

        model = model_dir.name.split("_", 1)[-1]

        for cg_dir in sorted(model_dir.iterdir()):
            if not cg_dir.is_dir():
                continue

            cg_scheme = cg_dir.name.split("_", 1)[-1]

            for partition_dir in sorted(cg_dir.iterdir()):
                if not partition_dir.is_dir():
                    continue

                partition = partition_dir.name

                log_path = partition_dir / f"{partition}_logs.csv"
                if log_path.exists():
                    yield model, cg_scheme, partition, str(log_path)


# ----------------------------- atom helpers -----------------------------

def build_atom_tables(atoms_df: pd.DataFrame,
                      require_complete_cell: bool = True,
                      require_positive_area: bool = True) -> tuple[pd.DataFrame, dict[int, tuple]]:
    """
    Adds:
      res_id       = (Chain, Residue Sequence, Residue)
      res_class    = protein/dna/other
      is_outer_atom

    Returns:
      atoms_df (augmented),
      atom_to_res_id: atom_index -> res_id
    """
    df = atoms_df.copy()

    df["res_id"] = list(zip(df["Chain"], df["Residue Sequence"], df["Residue"]))
    df["res_class"] = df["Residue"].apply(classify_residue)

    outer = ~df["Inner Ball?"].astype(bool)

    if require_complete_cell and "Complete Cell?" in df.columns:
        outer = outer & df["Complete Cell?"].astype(bool)

    if require_positive_area and "Surface Area" in df.columns:
        outer = outer & (df["Surface Area"].astype(float) > 0.0)

    df["is_outer_atom"] = outer

    atom_to_res_id = dict(zip(df["Index"].astype(int).tolist(), df["res_id"].tolist()))

    return df, atom_to_res_id


def _empty_surface_df() -> pd.DataFrame:
    return pd.DataFrame(columns=[
        "surface_kind",
        "res_i_id", "res_i_class",
        "res_j_id", "res_j_class",
        "pair_class",
        "area", "H_mean", "K_mean",
        "n_terms",
    ])


def _classify_ball(ball_idx: int,
                   atom_to_res_id: dict[int, tuple],
                   res_class_map: dict[tuple, str],
                   outer_atom_map: dict[int, bool]) -> tuple[tuple | str | None, str | None, bool]:
    """
    Returns:
      entity_id, entity_class, is_valid_outer_side

    Conventions:
      - if ball_idx exists in atom_to_res_id -> atomic/residue side
      - if ball_idx does not exist in atom_to_res_id -> solvent side

    This assumes solvent pseudo-balls are not present in atoms_df.
    If your solvent encoding differs, this is the one helper to edit.
    """
    if ball_idx in atom_to_res_id:
        res_id = atom_to_res_id[ball_idx]
        cls = res_class_map.get(res_id, "other")
        is_valid_outer = bool(outer_atom_map.get(ball_idx, False))

        return res_id, cls, is_valid_outer

    return "SOLV", "solv", True


def _canonicalize_surface_pair(entity_a,
                               class_a: str,
                               entity_b,
                               class_b: str) -> tuple[object, str, object, str]:
    """
    Canonical ordering so equivalent surfaces aggregate together.

    Rules:
      - residue/residue pairs sorted lexicographically by entity id
      - residue/solvent keeps residue first, solvent second
      - solvent/solvent stays solvent/solvent
    """
    if class_a == "solv" and class_b != "solv":
        return entity_b, class_b, entity_a, class_a

    if class_b == "solv" and class_a != "solv":
        return entity_a, class_a, entity_b, class_b

    if entity_a <= entity_b:
        return entity_a, class_a, entity_b, class_b

    return entity_b, class_b, entity_a, class_a


def _pair_class(ci: str, cj: str) -> str:
    pair = {ci, cj}

    if ci == "protein" and cj == "protein":
        return "P-P"

    if ci == "dna" and cj == "dna":
        return "D-D"

    if pair == {"protein", "dna"}:
        return "P-D"

    if pair == {"protein", "solv"}:
        return "P-S"

    if pair == {"dna", "solv"}:
        return "D-S"

    if ci == "solv" and cj == "solv":
        return "S-S"

    return "O-O"


# ----------------------------- core surface builder -----------------------------

def build_surface_table(surfs_df: pd.DataFrame,
                        atoms_df: pd.DataFrame,
                        atom_to_res_id: dict[int, tuple],
                        require_contact_area: bool = False,
                        include_ss: bool = False) -> pd.DataFrame:
    """
    Build a unified SURFACE-LEVEL interface table from surfs_df.

    Every row in the output is an aggregation of one or more surfaces with the same
    canonical entity pair:
      - P-P
      - P-D
      - D-D
      - P-S
      - D-S
      - optionally S-S

    Filtering:
      - positive surface area
      - optional positive contact area
      - atomic sides must be outer atoms
      - solvent sides are always allowed
    """
    if surfs_df is None or len(surfs_df) == 0:
        return _empty_surface_df()

    outer_atom_map = dict(zip(
        atoms_df["Index"].astype(int).tolist(),
        atoms_df["is_outer_atom"].astype(bool).tolist()
    ))

    res_class_map = dict(zip(
        atoms_df["res_id"].tolist(),
        atoms_df["res_class"].tolist()
    ))

    rows = []

    for _, r in surfs_df.iterrows():
        balls = r["Balls"]
        if balls is None or len(balls) != 2:
            continue

        a = int(balls[0])
        b = int(balls[1])

        area = float(r["Surface Area"])
        if area <= 0.0:
            continue

        if require_contact_area and "Contact Area" in r and float(r["Contact Area"]) <= 0.0:
            continue

        entity_a, class_a, valid_a = _classify_ball(
            a,
            atom_to_res_id=atom_to_res_id,
            res_class_map=res_class_map,
            outer_atom_map=outer_atom_map,
        )

        entity_b, class_b, valid_b = _classify_ball(
            b,
            atom_to_res_id=atom_to_res_id,
            res_class_map=res_class_map,
            outer_atom_map=outer_atom_map,
        )

        if not valid_a or not valid_b:
            continue

        pair_class = _pair_class(class_a, class_b)

        if pair_class == "O-O":
            continue

        if pair_class == "S-S" and not include_ss:
            continue

        res_i, res_i_class, res_j, res_j_class = _canonicalize_surface_pair(
            entity_a, class_a, entity_b, class_b
        )

        H = float(r["Mean Curvature"])
        K = float(r["Gauss Curvature"]) if "Gauss Curvature" in r else np.nan

        rows.append((
            "surface",
            res_i, res_i_class,
            res_j, res_j_class,
            pair_class,
            area, H, K
        ))

    if len(rows) == 0:
        return _empty_surface_df()

    tmp = pd.DataFrame(rows, columns=[
        "surface_kind",
        "res_i_id", "res_i_class",
        "res_j_id", "res_j_class",
        "pair_class",
        "area", "H", "K",
    ])

    tmp["aH"] = tmp["area"] * tmp["H"]
    tmp["aK"] = tmp["area"] * tmp["K"]

    out = tmp.groupby(
        ["surface_kind", "res_i_id", "res_i_class", "res_j_id", "res_j_class", "pair_class"],
        dropna=False
    ).agg(
        area=("area", "sum"),
        aH=("aH", "sum"),
        aK=("aK", "sum"),
        n_terms=("H", "count"),
    ).reset_index()

    out["H_mean"] = out["aH"] / out["area"]
    out["K_mean"] = out["aK"] / out["area"]

    return out[[
        "surface_kind",
        "res_i_id", "res_i_class",
        "res_j_id", "res_j_class",
        "pair_class",
        "area", "H_mean", "K_mean",
        "n_terms",
    ]]


def build_surf_df_for_log(log_path: str,
                          model: str,
                          cg_scheme: str,
                          partition: str,
                          require_complete_cell: bool = True,
                          require_positive_area: bool = True,
                          require_contact_area: bool = False,
                          include_ss: bool = False) -> pd.DataFrame:
    logs_obj = read_logs2(
        log_path,
        return_dict=False,
        all_=True,
        balls=True,
        surfs=True,
        edges=False,
        verts=False,
    )

    atoms_df, atom_to_res_id = build_atom_tables(
        logs_obj["atoms"],
        require_complete_cell=require_complete_cell,
        require_positive_area=require_positive_area,
    )

    surf_df = build_surface_table(
        logs_obj["surfs"],
        atoms_df,
        atom_to_res_id,
        require_contact_area=require_contact_area,
        include_ss=include_ss,
    )

    surf_df["model"] = model
    surf_df["cg_scheme"] = cg_scheme
    surf_df["partition"] = partition
    surf_df["log_path"] = log_path

    return surf_df


# ----------------------------- public entrypoint -----------------------------

def build_all_surfaces(main_systems_dir: str,
                       include_models: set[str] | None = None,
                       include_cg_schemes: set[str] | None = None,
                       include_partitions: set[str] | None = None,
                       out_dir: str = "figure8_outputs",
                       include_ss: bool = False) -> str:
    """
    Walks the directory tree, builds a unified surf_df, and writes surf_df.csv.
    Returns the output CSV path.
    """
    rows = []

    for model, cg_scheme, partition, log_path in iter_log_files(main_systems_dir):
        if include_models is not None and model not in include_models:
            continue

        if include_cg_schemes is not None and cg_scheme not in include_cg_schemes:
            continue

        if include_partitions is not None and partition not in include_partitions:
            continue

        print(f"[F8] model={model} cg={cg_scheme} part={partition} | {log_path}")

        surf_df = build_surf_df_for_log(
            log_path=log_path,
            model=model,
            cg_scheme=cg_scheme,
            partition=partition,
            require_complete_cell=True,
            require_positive_area=True,
            require_contact_area=False,
            include_ss=include_ss,
        )

        qa = surf_df.groupby(["pair_class"], dropna=False).agg(
            n=("pair_class", "count"),
            area=("area", "sum"),
        )

        print(qa)

        rows.append(surf_df)

    if len(rows) == 0:
        raise RuntimeError("No logs found / matched filters.")

    out = pd.concat(rows, ignore_index=True, sort=False)

    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    csv_path = out_path / "surf_df.csv"
    out.to_csv(csv_path, index=False)

    print(f"[F8] wrote {csv_path} | rows={len(out)}")

    return str(csv_path)


# ----------------------------- script runner -----------------------------

def main() -> None:
    root = Tk()
    root.withdraw()

    main_systems_dir = filedialog.askdirectory(title="Select Main_Systems folder")
    if not main_systems_dir:
        print("[F8] No folder selected. Exiting.")
        return

    csv_path = build_all_surfaces(
        main_systems_dir=main_systems_dir,
        include_models={"NCP"},
        include_cg_schemes=None,
        include_partitions={"aw", "pow", "prm"},
        out_dir="figure8_outputs",
        include_ss=False,
    )

    print(f"[F8] Done. Output: {csv_path}")


if __name__ == "__main__":

    main()
