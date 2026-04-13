import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize


# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
sys.path.append(vorpy_root)

from vorpy.src.calculations.surf import calc_surf_func
from vorpy.src.calculations.curvature import mean_curvature
from vorpy.src.analyze.tools.plot_templates.line import line_plot


def curvature_distance_ratio_sweep(
    r1: float = 1.0,
    surface_distance_min: float = -1.0,
    surface_distance_max: float = 2.0,
    n_distance_lines: int = 10,
    ratio_min: float = 1.0,
    ratio_max: float = 10.0,
    n_ratio_points: int = 40,
    skip_equal_radii: bool = True,
    fill_value: float = np.nan,
    verbose: bool = True,
) -> dict:
    """
    Sweep mean curvature at the midpoint between the two sphere surfaces as a function of:
      - surface_distance = (center_distance - r1 - r2)  [negative => overlap]
      - ratio = r2 / r1

    Returns a dict with:
      {
        "surface_distances": np.ndarray shape (n_distance_lines,),
        "ratios": np.ndarray shape (n_ratio_points,),
        "curvatures": np.ndarray shape (n_distance_lines, n_ratio_points),
      }
    """

    if n_distance_lines < 1:
        raise ValueError("n_distance_lines must be >= 1")

    if n_ratio_points < 2:
        raise ValueError("n_ratio_points must be >= 2")

    if ratio_min <= 0 or ratio_max <= 0:
        raise ValueError("ratio_min and ratio_max must be > 0")

    if ratio_max < ratio_min:
        raise ValueError("ratio_max must be >= ratio_min")

    surface_distances = np.linspace(surface_distance_min, surface_distance_max, n_distance_lines)
    ratios = np.linspace(ratio_min, ratio_max, n_ratio_points)

    curvatures = np.full((n_distance_lines, n_ratio_points), fill_value, dtype=float)

    total = n_distance_lines * n_ratio_points
    count = 0

    for j, sdist in enumerate(surface_distances):
        # Midpoint between surfaces (along x) measured from small sphere center
        eval_x = r1 + 0.5 * sdist
        eval_point = np.array([eval_x, 0.0, 0.0], dtype=float)

        for i, ratio in enumerate(ratios):
            r2 = r1 * ratio

            if skip_equal_radii and np.isclose(r2, r1):
                curvatures[j, i] = fill_value
                count += 1
                if verbose:
                    print(f"\rProgress: {100.0 * count / total:6.2f} %", end="")
                continue

            # Center of the large sphere so that surface_distance = (d - r1 - r2)
            large_center_x = r1 + r2 + sdist

            func = calc_surf_func([0, 0, 0], r1, [large_center_x, 0, 0], r2)

            try:
                curvatures[j, i] = mean_curvature(func, eval_point)
            except Exception:
                # Keep fill_value if the surface function / curvature eval blows up
                curvatures[j, i] = fill_value

            count += 1
            if verbose:
                print(f"\rProgress: {100.0 * count / total:6.2f} %", end="")

    if verbose:
        print()

    return {
        "surface_distances": surface_distances,
        "ratios": ratios,
        "curvatures": curvatures,
    }


def plot_curvature_distance_ratio_lines(
    sweep: dict,
    title: str = "Curvature by Distance and Radii",
    cmap=plt.cm.rainbow,
    Show: bool = True,
    save: str | None = None,
    figsize: tuple[float, float] = (8, 6),
    linewidth: float = 3.0,
    tick_val_size: int = 30,
    x_label_size: int = 30,
    y_label_size: int = 30,
    legend_label_size: int = 30,
    legend_title_size: int = 30,
) -> None:
    """
    Plot each surface-distance line as curvature vs radius ratio.
    Uses your existing line_plot() template and adds a colorbar keyed to surface distance.
    """

    surface_distances = sweep["surface_distances"]
    ratios = sweep["ratios"]
    curvatures = sweep["curvatures"]

    norm = Normalize(vmin=float(np.min(surface_distances)), vmax=float(np.max(surface_distances)))
    sm = ScalarMappable(norm=norm, cmap=cmap)

    line_colors = [cmap(norm(float(d))) for d in surface_distances]

    xs = [ratios.tolist() for _ in surface_distances]
    ys = [curvatures[k, :].tolist() for k in range(len(surface_distances))]

    line_plot(
        xs,
        ys,
        colors=line_colors,
        Show=Show,
        x_label="Radii Ratio",
        y_label="Mean Curvature",
        labels=[float(d) for d in surface_distances],
        legend_title="Surface Distance",
        legend_label_size=legend_label_size,
        legend_title_size=legend_title_size,
        title=title,
        tick_val_size=tick_val_size,
        x_label_size=x_label_size,
        y_label_size=y_label_size,
        linewidth=linewidth,
        colorbar=sm,
        figsize=figsize,
        x_ticks=[1, 2, 3],
        y_ticks=[0.0, 0.5, 1.0, 1.5, 2.0],
        xlim=[0.8, 3.2],
        ylim=[-0.1, 2.1]
    )


if __name__ == "__main__":
    sweep = curvature_distance_ratio_sweep(
        r1=1.0,
        surface_distance_min=-1.5,
        surface_distance_max=1.5,
        n_distance_lines=15,
        ratio_min=1.0,
        ratio_max=3.0,
        n_ratio_points=40,
        skip_equal_radii=True,
        verbose=True,
    )

    plot_curvature_distance_ratio_lines(
        sweep,
        title="Curvature by Distance and Radii",
        Show=True,
    )
