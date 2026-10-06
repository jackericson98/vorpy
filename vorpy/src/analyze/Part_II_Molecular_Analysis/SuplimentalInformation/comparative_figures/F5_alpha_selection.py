from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from .plot_common import (
    configure_matplotlib,
    load_static_comparison,
    method_label,
    ordered_systems,
    require_methods,
    save_figure,
    select_output_directory,
    static_protein_only,
)


METHODS = [
    "full_power",
    "alpha0_power_before_M",
    "cazals_alpha0_power",
]


def _extract_progression(df, system: str, column: str) -> list[float]:
    values = []

    for method in METHODS:
        rows = df[
            (df["system"] == system)
            & (df["method"] == method)
        ]

        if len(rows) != 1:
            raise ValueError(
                f"Expected exactly one row for {system} / {method}; "
                f"found {len(rows)}."
            )

        values.append(float(rows.iloc[0][column]))

    return values


def _plot_absolute(df, systems, output_dir):
    x = np.arange(len(METHODS))
    xlabels = [method_label(method) for method in METHODS]

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(14.0, 4.5),
    )

    specifications = [
        (
            "selected_elements",
            "Selected interface elements",
            "Interface elements",
        ),
        (
            "interface_area_A2",
            r"Interface area",
            r"Area ($\AA^2$)",
        ),
        (
            "combined_H_A",
            r"Integrated mean curvature",
            r"$C_H$ ($\AA$)",
        ),
    ]

    for ax, (column, title, ylabel) in zip(axes, specifications):
        for system in systems:
            y = _extract_progression(df, system, column)

            ax.plot(
                x,
                y,
                marker="o",
                linewidth=2,
                markersize=6,
                label=system,
            )

        ax.set_xticks(x)
        ax.set_xticklabels(xlabels)
        ax.set_title(title)
        ax.set_ylabel(ylabel)

        if column == "combined_H_A":
            ax.axhline(0.0, linewidth=0.8)

        ax.grid(axis="y", alpha=0.2)

    axes[0].legend(frameon=False)

    fig.suptitle(
        "Effect of Cazals selection on the Power interface"
    )

    fig.tight_layout()

    paths = save_figure(
        fig,
        output_dir,
        "F5a_alpha_selection_absolute",
    )

    plt.close(fig)
    return paths


def _plot_retention(df, systems, output_dir):
    x = np.arange(len(METHODS))
    xlabels = [method_label(method) for method in METHODS]

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(9.5, 4.5),
    )

    specifications = [
        (
            "selected_elements",
            "Interface-element retention",
        ),
        (
            "interface_area_A2",
            "Interface-area retention",
        ),
    ]

    for ax, (column, title) in zip(axes, specifications):
        for system in systems:
            raw = np.asarray(
                _extract_progression(df, system, column),
                dtype=float,
            )

            baseline = raw[0]

            if not np.isfinite(baseline) or baseline <= 0:
                raise ValueError(
                    f"Cannot normalize {column} for {system}: "
                    f"Full Power baseline = {baseline}"
                )

            retained = 100.0 * raw / baseline

            ax.plot(
                x,
                retained,
                marker="o",
                linewidth=2,
                markersize=6,
                label=system,
            )

            # Give exact final alpha0+M5 retention beside the last point.
            ax.annotate(
                f"{retained[-1]:.1f}%",
                (x[-1], retained[-1]),
                xytext=(5, 0),
                textcoords="offset points",
                va="center",
                fontsize=9,
            )

        ax.set_xticks(x)
        ax.set_xticklabels(xlabels)
        ax.set_title(title)
        ax.set_ylabel("Retained from Full Power (%)")
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", alpha=0.2)

    axes[0].legend(frameon=False)

    fig.suptitle(
        "Fraction of the Full Power interface retained by Cazals selection"
    )

    fig.tight_layout()

    paths = save_figure(
        fig,
        output_dir,
        "F5b_alpha_selection_retention",
    )

    plt.close(fig)
    return paths


def _print_summary(df, systems):
    print("\nALPHA-SELECTION SUMMARY")
    print("=" * 72)

    for system in systems:
        full_elements, alpha_elements, final_elements = _extract_progression(
            df,
            system,
            "selected_elements",
        )

        full_area, alpha_area, final_area = _extract_progression(
            df,
            system,
            "interface_area_A2",
        )

        print(f"\n{system}")
        print(
            "  Interface elements: "
            f"{full_elements:.0f} -> "
            f"{alpha_elements:.0f} -> "
            f"{final_elements:.0f}"
        )
        print(
            "  Element retention:   "
            f"{100.0 * final_elements / full_elements:.2f}%"
        )
        print(
            "  Interface area:       "
            f"{full_area:.3f} -> "
            f"{alpha_area:.3f} -> "
            f"{final_area:.3f} A^2"
        )
        print(
            "  Area retention:       "
            f"{100.0 * final_area / full_area:.2f}%"
        )

        if np.isclose(alpha_elements, final_elements) and np.isclose(
            alpha_area,
            final_area,
        ):
            print(
                "  M=5 effect:           "
                "none for these reported metrics"
            )


def main():
    configure_matplotlib()

    df, input_path = load_static_comparison()
    df = static_protein_only(df)

    systems = ordered_systems(df)

    require_methods(
        df,
        METHODS,
        systems=systems,
    )

    output_dir = select_output_directory(
        initial_dir=input_path.parent,
        title="Select output directory for F5 figures",
    )

    _print_summary(df, systems)

    absolute_paths = _plot_absolute(
        df,
        systems,
        output_dir,
    )

    retention_paths = _plot_retention(
        df,
        systems,
        output_dir,
    )

    print("\nFIGURES WRITTEN")
    print("=" * 72)

    for path in absolute_paths + retention_paths:
        print(path)


if __name__ == "__main__":
    main()