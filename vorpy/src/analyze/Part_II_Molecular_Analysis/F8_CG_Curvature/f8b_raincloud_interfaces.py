import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from tkinter import Tk
from tkinter import filedialog

from matplotlib.colors import to_rgb, to_hex

try:
    from scipy.stats import gaussian_kde
except Exception:
    gaussian_kde = None


# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
sys.path.append(vorpy_root)

from vorpy.src.analyze.Part_II_Molecular_Analysis.F8_CG_Curvature.f8_build_surfaces import build_all_surfaces

SHOW_MODES = True
SHOW_REGIMES = False

PAIR_COLORS = {
    "P-P": "#1f77b4",
    "P-D": "#ff7f0e",
    "D-D": "#2ca02c",
    "P-S": "#9467bd",
    "D-S": "#8c564b",
}


def add_regime_shading(ax, y_max):
    """
    Adds curvature regime shading to the background.
    """
    LOW_MAX = 0.15
    HIGH_MIN = 0.25

    ax.axhspan(0.0, LOW_MAX, alpha=0.08, color='blue', zorder=0)
    ax.axhspan(LOW_MAX, HIGH_MIN, alpha=0.05, color='gray', zorder=0)
    ax.axhspan(HIGH_MIN, y_max, alpha=0.08, color='red', zorder=0)


def compute_dual_modes(vals):
    """
    Compute low/high curvature modes using threshold split.
    Returns (low_mode, high_mode, low_frac, high_frac)
    """
    vals = vals[np.isfinite(vals)]
    if len(vals) < 10:
        return None, None, None, None

    low = vals[vals < 0.15]
    high = vals[vals > 0.25]

    if len(low) < 5 or len(high) < 5:
        return None, None, None, None

    low_mode = np.median(low)
    high_mode = np.median(high)

    low_frac = len(low) / len(vals)
    high_frac = len(high) / len(vals)

    return low_mode, high_mode, low_frac, high_frac


def plot_dual_mode_markers(ax, x_pos, vals, color):
    """
    Plot dual-mode markers (circle = low, diamond = high)
    """
    low_mode, high_mode, low_frac, high_frac = compute_dual_modes(vals)

    if low_mode is None:
        return

    # scale marker size by fraction
    base_size = 80

    ax.scatter(
        x_pos,
        low_mode,
        s=base_size * (0.5 + low_frac),
        color=color,
        edgecolor='black',
        linewidth=1.2,
        zorder=6,
        marker='o'
    )

    ax.scatter(
        x_pos,
        high_mode,
        s=base_size * (0.5 + high_frac),
        color=color,
        edgecolor='black',
        linewidth=1.2,
        zorder=6,
        marker='D'
    )


def adjust_color(color: str, factor: float) -> str:
    """
    factor > 1.0  -> lighten
    factor < 1.0  -> darken
    """
    r, g, b = to_rgb(color)

    if factor >= 1.0:
        r = 1 - (1 - r) / factor
        g = 1 - (1 - g) / factor
        b = 1 - (1 - b) / factor
    else:
        r = r * factor
        g = g * factor
        b = b * factor

    return to_hex((r, g, b))


PAIR_COLOR_STYLE = {
    cls: {
        "violin": base,
        "box": adjust_color(base, 0.75),
        "scatter": adjust_color(base, 1.45),
    }
    for cls, base in PAIR_COLORS.items()
}


def _trim_to_y_range(vals: np.ndarray,
                     areas: np.ndarray | None,
                     y_range: tuple[float, float] | None):
    """
    Trim values (and matching areas) to the visible y-range.
    """
    vals = np.asarray(vals, dtype=float)
    ok = np.isfinite(vals)

    if y_range is not None:
        y0, y1 = float(y_range[0]), float(y_range[1])
        ok = ok & (vals >= y0) & (vals <= y1)

    vals_out = vals[ok]

    if areas is None:
        return vals_out, None

    areas = np.asarray(areas, dtype=float)
    areas_out = areas[ok]

    return vals_out, areas_out


def _square_grid_shape(n_panels: int) -> tuple[int, int]:
    """
    Choose a near-square grid.

    Examples:
      1 -> (1, 1)
      2 -> (1, 2)
      3 -> (2, 2)
      4 -> (2, 2)
      5 -> (2, 3)
      6 -> (2, 3)
      7 -> (3, 3)
      8 -> (3, 3)
      9 -> (3, 3)
    """
    if n_panels <= 0:
        return 1, 1

    ncols = int(np.ceil(np.sqrt(n_panels)))
    nrows = int(np.ceil(n_panels / ncols))

    return nrows, ncols


def _ensure_surf_df(main_systems_dir: str,
                    out_dir: str,
                    include_models: set[str] | None,
                    include_partitions: set[str] | None,
                    include_cg_schemes: set[str] | None) -> str:
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    csv_path = out_path / "surf_df.csv"
    if csv_path.exists():
        return str(csv_path)

    return build_all_surfaces(
        main_systems_dir=main_systems_dir,
        include_models=include_models,
        include_partitions=include_partitions,
        include_cg_schemes=include_cg_schemes,
        out_dir=out_dir,
        include_ss=False,
    )


def _cg_order_key(name: str) -> tuple:
    order = {
        "Atom": 1,
        "Encap": 2,
        "Encap_SR": 3,
        "AD": 4,
        "AD_SR": 5,
        "AD_MW": 6,
        "AD_MW_SR": 7,
    }

    if name in order:
        return 0, order[name]

    return 1, str(name)


def _resample_weighted(values: np.ndarray,
                       weights: np.ndarray,
                       n_sample: int,
                       rng: np.random.Generator) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)

    ok = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    values = values[ok]
    weights = weights[ok]

    if values.size == 0:
        return values

    if values.size <= 3:
        return values

    p = weights / np.sum(weights)
    idx = rng.choice(
        np.arange(values.size),
        size=min(n_sample, max(values.size, 10)),
        replace=True,
        p=p
    )

    return values[idx]


def _half_violin(ax,
                 data: np.ndarray,
                 x0: float,
                 side: str = "left",
                 width: float = 0.35,
                 bw_adjust: float = 1.0,
                 n_grid: int = 200,
                 alpha: float = 0.35,
                 color: str = "#333333",
                 kde_jitter_frac: float = 1e-3,
                 fallback_bins: int = 30) -> None:

    data = np.asarray(data, dtype=float)
    data = data[np.isfinite(data)]

    if data.size < 5:
        return

    y_min = float(np.min(data))
    y_max = float(np.max(data))
    span = y_max - y_min

    if span <= 1e-12:
        ys = np.linspace(y_min - 1e-6, y_max + 1e-6, n_grid)
        dens = np.zeros_like(ys)
        dens[len(dens) // 2] = 1.0
        dens = dens / np.max(dens) * width

        if side.lower() == "left":
            xs = x0 - dens
            poly = ax.fill_betweenx(
                ys,
                xs,
                x0,
                alpha=alpha,
                color=color
            )
            poly.set_clip_path(ax.patch)
        else:
            xs = x0 + dens
            poly = ax.fill_betweenx(
                ys,
                xs,
                x0,
                alpha=alpha,
                color=color
            )
            poly.set_clip_path(ax.patch)

        return

    pad = 0.05 * (span + 1e-12)
    ys = np.linspace(y_min - pad, y_max + pad, n_grid)

    def _draw_from_density(dens_vals: np.ndarray):
        if np.max(dens_vals) <= 0:
            return

        dens_vals = dens_vals / np.max(dens_vals) * width

        if side.lower() == "left":
            xs = x0 - dens_vals
            poly = ax.fill_betweenx(ys, xs, x0, alpha=alpha, color=color)
            poly.set_clip_path(ax.patch)
        else:
            xs = x0 + dens_vals
            poly = ax.fill_betweenx(ys, x0, xs, alpha=alpha, color=color)
            poly.set_clip_path(ax.patch)

    if gaussian_kde is not None:
        jitter_scale = float(kde_jitter_frac) * span
        data_kde = data.copy()

        if jitter_scale > 0:
            rng_local = np.random.default_rng(0)
            data_kde = data_kde + rng_local.normal(loc=0.0, scale=jitter_scale, size=data_kde.size)

        try:
            kde = gaussian_kde(data_kde)
            kde.set_bandwidth(kde.factor * float(bw_adjust))
            dens = kde(ys)
            _draw_from_density(dens)
            return
        except Exception:
            pass

    hist, bin_edges = np.histogram(data, bins=fallback_bins, density=True)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    dens = np.interp(ys, bin_centers, hist, left=0.0, right=0.0)
    _draw_from_density(dens)



def _boxplot(ax,
             data: np.ndarray,
             x0: float,
             box_width: float = 0.18,
             color: str = "#333333") -> None:
    if data.size == 0:
        return

    ax.boxplot(
        [data],
        positions=[x0],
        widths=box_width,
        vert=True,
        showfliers=False,
        patch_artist=True,
        boxprops={"facecolor": color, "linewidth": 2.0, "alpha": 0.6},
        medianprops={"linewidth": 2.0, "color": "black"},
        whiskerprops={"linewidth": 2.0},
        capprops={"linewidth": 2.0},
    )



def _jitter_scatter(ax,
                    data: np.ndarray,
                    x0: float,
                    side: str = "right",
                    jitter: float = 0.12,
                    size: float = 10.0,
                    alpha: float = 0.45,
                    rng: np.random.Generator | None = None,
                    max_points: int = 2500,
                    color: str = "#333333") -> None:
    if data.size == 0:
        return

    rng = rng or np.random.default_rng(0)

    if data.size > max_points:
        idx = rng.choice(np.arange(data.size), size=max_points, replace=False)
        data = data[idx]

    j = rng.uniform(0.0, jitter, size=data.size)
    if side.lower() == "left":
        xs = x0 - j
    else:
        xs = x0 + j

    ax.scatter(xs, data, s=size, alpha=alpha, edgecolors="none", color=color)



def _draw_raincloud_group(ax,
                          x0: float,
                          vals: np.ndarray,
                          areas: np.ndarray | None,
                          weight_by_area: bool,
                          weighted_n_sample: int,
                          rng: np.random.Generator,
                          violin_width: float = 0.35,
                          box_width: float = 0.18,
                          jitter: float = 0.12,
                          point_size: float = 10.0,
                          color_violin: str = "#333333",
                          color_box: str = "#222222",
                          color_scatter: str = "#555555",
                          y_range: tuple[float, float] | None = None,
                          gap: float = 0.04) -> np.ndarray:
    vals_trim, areas_trim = _trim_to_y_range(vals, areas, y_range)

    if vals_trim.size == 0:
        return vals_trim

    vals_raw = vals_trim[np.isfinite(vals_trim)]

    if weight_by_area and areas_trim is not None:
        vals_plot = _resample_weighted(vals_trim, areas_trim, n_sample=weighted_n_sample, rng=rng)
        if np.unique(vals_plot).size < 10:
            vals_plot = vals_raw
    else:
        vals_plot = vals_raw

    _half_violin(
        ax,
        vals_plot,
        x0 - gap,
        side="left",
        width=violin_width,
        bw_adjust=1.0,
        alpha=0.35,
        color=color_violin
    )

    # _boxplot(
    #     ax,
    #     vals_raw,
    #     x0,
    #     box_width=box_width,
    #     color=color_box
    # )

    _jitter_scatter(
        ax,
        vals_plot,
        x0 + gap,
        side="right",
        jitter=jitter,
        size=point_size,
        alpha=0.45,
        rng=rng,
        color=color_scatter
    )

    return vals_plot


def _resolve_cg_schemes(all_cg_schemes: list[str],
                        include_cg_schemes: list[str] | None) -> list[str]:
    """
    Returns CG schemes in requested order if provided, otherwise in standard order.
    Ignores requested schemes that are not present in the data.
    """
    all_cg_schemes = [str(x) for x in all_cg_schemes]

    if include_cg_schemes is None:
        return sorted(all_cg_schemes, key=_cg_order_key)

    requested = [str(x) for x in include_cg_schemes]
    present = set(all_cg_schemes)

    return [cg for cg in requested if cg in present]


def plot_f8b_raincloud_interfaces(surf_df: pd.DataFrame,
                                  out_dir: str,
                                  model: str = "NCP",
                                  include_cg_schemes: list[str] | None = None,
                                  pair_classes: list[str] | None = None,
                                  curvature_col: str = "H_mean",
                                  weight_by_area: bool = True,
                                  weighted_n_sample: int = 6000,
                                  fixed_y: bool = False,
                                  y_range: tuple[float, float] | None = None,
                                  title: str = "Figure 8B | AW interface curvature by CG scheme",
                                  rng_seed: int = 7,
                                  save_name: str = "F8B_raincloud_interfaces_AW.png",
                                  show: bool = True,
                                  layout: str = "grid") -> str:
    out_path = Path(out_dir)
    (out_path / "plots").mkdir(parents=True, exist_ok=True)
    save_path = str(out_path / "plots" / save_name)

    if pair_classes is None:
        pair_classes = ["P-P", "P-D", "D-D", "P-S", "D-S"]

    df = surf_df.copy()

    df = df[df["model"].astype(str) == str(model)]
    df = df[df["partition"].astype(str) == "aw"]
    df = df[df["surface_kind"].astype(str) == "surface"]
    df = df[df["pair_class"].isin(pair_classes)]
    df = df[np.isfinite(df[curvature_col].astype(float))]
    df = df[np.isfinite(df["area"].astype(float))]
    df = df[df["area"].astype(float) > 0.0]

    if include_cg_schemes is not None:
        df = df[df["cg_scheme"].isin(include_cg_schemes)]

    if df.empty:
        raise RuntimeError("No rows after filtering. Check model/cg_scheme filters and AW data.")

    cg_schemes = _resolve_cg_schemes(
        all_cg_schemes=df["cg_scheme"].astype(str).unique().tolist(),
        include_cg_schemes=include_cg_schemes,
    )

    if len(cg_schemes) == 0:
        raise RuntimeError("No CG schemes remain after applying include_cg_schemes filter.")

    rng = np.random.default_rng(rng_seed)
    all_vals = []

    layout = str(layout).strip().lower()

    if layout == "grid":
        n = len(cg_schemes)
        nrows, ncols = _square_grid_shape(n)

        fig_w = 5.2 * ncols
        fig_h = 4.6 * nrows

        fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h), sharey=True)
        axes = np.array(axes).reshape(-1)



        for ax_i, ax in enumerate(axes):
            if ax_i >= n:
                ax.axis("off")
                continue

            if SHOW_REGIMES:
                y_max = ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 0.7
                add_regime_shading(ax, y_max)

            cg = cg_schemes[ax_i]
            sub = df[df["cg_scheme"].astype(str) == cg].copy()

            ax.set_title(cg, fontsize=16)

            x_positions = np.arange(len(pair_classes), dtype=float)

            for x0, cls in zip(x_positions, pair_classes):
                d = sub[sub["pair_class"] == cls]
                vals = d[curvature_col].astype(float).to_numpy()
                areas = d["area"].astype(float).to_numpy()

                style = PAIR_COLOR_STYLE.get(cls, {
                    "violin": "#777777",
                    "box": "#555555",
                    "scatter": "#999999",
                })
                print(f"\n=== {cg} | {cls} ===")

                vals_raw = vals[np.isfinite(vals)]
                areas_raw = areas[np.isfinite(vals)]

                if len(vals_raw) == 0:
                    print("No data")
                else:
                    print(f"N surfaces: {len(vals_raw)}")
                    print(f"Mean curvature: {np.mean(vals_raw):.4f}")
                    print(f"Median curvature: {np.median(vals_raw):.4f}")

                    print("Percentiles:")
                    for p in [5, 25, 50, 75, 95]:
                        print(f"  p{p}: {np.percentile(vals_raw, p):.4f}")

                    print(f"Min / Max: {np.min(vals_raw):.4f} / {np.max(vals_raw):.4f}")

                    # Check bimodality hint
                    hist, bins = np.histogram(vals_raw, bins=20)
                    print("Histogram counts:", hist.tolist())

                vals_plot = _draw_raincloud_group(
                    ax=ax,
                    x0=x0,
                    vals=vals,
                    areas=areas,
                    weight_by_area=weight_by_area,
                    weighted_n_sample=weighted_n_sample,
                    rng=rng,
                    violin_width=0.32,
                    box_width=0.16,
                    jitter=0.10,
                    point_size=9.0,
                    color_violin=style["violin"],
                    color_box=style["box"],
                    color_scatter=style["scatter"],
                    y_range=y_range if fixed_y else None,
                )
                all_vals.append(vals_plot)

                if SHOW_MODES:
                    vals_mode, _ = _trim_to_y_range(
                        vals,
                        None,
                        y_range if fixed_y else None
                    )

                    plot_dual_mode_markers(
                        ax=ax,
                        x_pos=x0,
                        vals=vals_mode,
                        color=style["box"]
                    )

            ax.set_xticks(x_positions)
            ax.set_xticklabels(pair_classes, fontsize=11, rotation=20)

            ax.tick_params(axis="both", which="major", labelsize=12, width=2, length=7)
            for spine in ax.spines.values():
                spine.set_linewidth(2)

        if fixed_y and y_range is not None:
            y0 = float(y_range[0])
            y1 = float(y_range[1])
        else:
            flat = np.concatenate([v for v in all_vals if v.size > 0]) if len(all_vals) else np.array([])
            if flat.size > 0:
                lo = float(np.nanpercentile(flat, 1.0))
                hi = float(np.nanpercentile(flat, 99.0))
                pad = 0.08 * (hi - lo + 1e-12)
                y0 = lo - pad
                y1 = hi + pad
            else:
                y0, y1 = -0.05, 1.05

        for ax in axes[:n]:
            ax.set_ylim(y0, y1)
            for coll in ax.collections:
                coll.set_clip_path(ax.patch)

        fig.suptitle(title, fontsize=18)
        fig.text(0.5, 0.04, "Interface class", ha="center", fontsize=16)
        fig.text(0.04, 0.5, curvature_col, va="center", rotation="vertical", fontsize=16)

        fig.tight_layout(rect=[0.06, 0.06, 1, 0.94])

    elif layout == "single":
        n = len(cg_schemes)

        group_gap = 2.2
        within_gap = 0.38

        fig_w = max(16.0, 2.5 * n)
        fig_h = 7.0

        fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h))


        x_ticks = []
        x_ticklabels = []

        x_base = 0.0

        for cg in cg_schemes:
            sub = df[df["cg_scheme"].astype(str) == cg].copy()
            xs = [x_base + i * within_gap for i in range(len(pair_classes))]

            for x0, cls in zip(xs, pair_classes):
                d = sub[sub["pair_class"] == cls]
                vals = d[curvature_col].astype(float).to_numpy()
                areas = d["area"].astype(float).to_numpy()

                style = PAIR_COLOR_STYLE.get(cls, {
                    "violin": "#777777",
                    "box": "#555555",
                    "scatter": "#999999",
                })

                vals_plot = _draw_raincloud_group(
                    ax=ax,
                    x0=x0,
                    vals=vals,
                    areas=areas,
                    weight_by_area=weight_by_area,
                    weighted_n_sample=weighted_n_sample,
                    rng=rng,
                    violin_width=0.32,
                    box_width=0.16,
                    jitter=0.10,
                    point_size=9.0,
                    color_violin=style["violin"],
                    color_box=style["box"],
                    color_scatter=style["scatter"],
                    y_range=y_range if fixed_y else None,
                )
                all_vals.append(vals_plot)

                if SHOW_MODES:
                    vals_mode, _ = _trim_to_y_range(
                        vals,
                        None,
                        y_range if fixed_y else None
                    )

                    plot_dual_mode_markers(
                        ax=ax,
                        x_pos=x0,
                        vals=vals_mode,
                        color=style["box"]
                    )

            group_center = float(np.mean(xs))
            x_ticks.append(group_center)
            x_ticklabels.append(cg)

            x_base = x_base + group_gap

        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticklabels, rotation=25, ha="right", fontsize=12)

        ax.set_xlabel("Coarse-graining scheme", fontsize=16)
        ax.set_ylabel(curvature_col, fontsize=16)

        ax.tick_params(axis="both", which="major", labelsize=12, width=2, length=7)
        for spine in ax.spines.values():
            spine.set_linewidth(2)

        ax.text(
            0.01, 0.99,
            "Within each CG group: left→right = P-P, P-D, D-D, P-S, D-S",
            transform=ax.transAxes,
            ha="left", va="top", fontsize=12
        )

        if fixed_y and y_range is not None:
            y0 = float(y_range[0])
            y1 = float(y_range[1])
        else:
            flat = np.concatenate([v for v in all_vals if v.size > 0]) if len(all_vals) else np.array([])
            if flat.size > 0:
                lo = float(np.nanpercentile(flat, 1.0))
                hi = float(np.nanpercentile(flat, 99.0))
                pad = 0.08 * (hi - lo + 1e-12)
                y0 = lo - pad
                y1 = hi + pad
            else:
                y0, y1 = -0.05, 1.05

        ax.set_ylim(y0, y1)
        for coll in ax.collections:
            coll.set_clip_path(ax.patch)

        fig.suptitle(title, fontsize=18)
        fig.tight_layout(rect=[0, 0, 1, 0.95])

    else:
        raise ValueError('layout must be "grid" or "single"')

    fig.savefig(save_path, dpi=300)
    if show:
        plt.show()
    plt.close(fig)

    print(f"[F8B] saved: {save_path}")
    return save_path


def main() -> None:
    root = Tk()
    root.withdraw()

    main_systems_dir = filedialog.askdirectory(title="Select Main_Systems folder")
    if not main_systems_dir:
        print("[F8B] No folder selected. Exiting.")
        return

    out_dir = "figure8_outputs"

    surf_csv = _ensure_surf_df(
        main_systems_dir=main_systems_dir,
        out_dir=out_dir,
        include_models={"NCP"},
        include_partitions={"aw"},
        include_cg_schemes=None,
    )

    surf_df = pd.read_csv(surf_csv)

    include_cg_schemes = ["Atom", "Encap_SR", "AD", "AD_MW", "AD_SR", "AD_MW_SR"]

    fixed_y = True
    y_range = (-0.05, 1.05)

    # layout = "single"
    layout = "grid"


    plot_f8b_raincloud_interfaces(
        surf_df=surf_df,
        out_dir=out_dir,
        model="NCP",
        include_cg_schemes=include_cg_schemes,
        pair_classes=["P-P", "P-D", "D-D", "P-S", "D-S"],
        curvature_col="H_mean",
        weight_by_area=True,
        weighted_n_sample=6000,
        fixed_y=fixed_y,
        y_range=y_range,
        title="Figure 8B | AW interface curvature by CG scheme",
        save_name=f"F8B_raincloud_interfaces_AW_{layout}.png",
        show=True,
        layout=layout,
    )


if __name__ == "__main__":


    main()