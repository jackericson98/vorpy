"""Publication plots generated only from persisted comparative-study CSVs."""
from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


SYSTEMS = ("2KAI", "1AK4", "1BRC")
METHODS = ("cazals_alpha0_power", "full_power", "aw", "full_power_aw_radii")
LABELS = {"cazals_alpha0_power": "Cazals alpha-0 Power", "full_power": "Full Power (Chothia radii)",
          "aw": "AW (frozen radii)", "full_power_aw_radii": "Full Power (AW radii)"}
COLORS = {"cazals_alpha0_power": "#2A8C55", "full_power": "#2864A5", "aw": "#D65332",
          "full_power_aw_radii": "#8064A2"}


def read(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def finite(rows, key, context):
    values = []
    for row in rows:
        value = row.get(key, "")
        if value == "":
            raise ValueError(f"Undefined {key} in {context}: {row}")
        number = float(value)
        if not math.isfinite(number):
            raise ValueError(f"Nonfinite {key} in {context}: {row}")
        values.append(number)
    return values


def save(fig, directory, name):
    directory.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(directory / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(directory / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(directory / f"{name}.svg", bbox_inches="tight")
    plt.close(fig)


def make_figures(output):
    output = Path(output)
    figures = output / "figures"
    static = read(output / "static_method_comparison.csv")
    trajectory = read(output / "trajectory_metrics.csv")
    summary = read(output / "trajectory_summary.csv")
    paired = read(output / "solvent_paired_effects.csv")
    if not static:
        raise ValueError("Static metrics CSV contains no rows")
    static2 = {r["method"]: r for r in static if r["system"] == "2KAI"}
    groups = [{r["method"]: r for r in static if r["system"] == system} for system in SYSTEMS]
    if any(not all(method in group for method in METHODS) for group in groups):
        raise ValueError("Static figure inputs must contain all four declared methods for every system")

    # A. Area
    fig, ax = plt.subplots(figsize=(10.4, 5.6))
    x = np.arange(len(SYSTEMS)); width = .19
    for offset, method in enumerate(METHODS):
        vals = [float(group[method]["interface_area_A2"]) for group in groups]
        ax.bar(x + (offset - 1.5) * width, vals, width, label=LABELS[method], color=COLORS[method])
    ax.set(xticks=x, xticklabels=SYSTEMS, ylabel=r"Interface area ($\mathrm{\AA}^2$)", title="Static interface area")
    ax.legend(frameon=False, ncol=2); ax.grid(axis="y", alpha=.2)
    save(fig, figures, "figure_A_static_interface_area")

    # B. Conventional integrated mean curvature
    fig, ax = plt.subplots(figsize=(10.4, 5.6))
    for offset, method in enumerate(METHODS):
        vals = [float(group[method]["combined_H_A"]) for group in groups]
        ax.bar(x + (offset - 1.5) * width, vals, width, label=LABELS[method], color=COLORS[method])
    ax.axhline(0, color="black", linewidth=.8)
    ax.set(xticks=x, xticklabels=SYSTEMS, ylabel=r"Signed conventional $C_H$ ($\mathrm{\AA}$)", title="Conventional integrated mean curvature")
    ax.legend(frameon=False, ncol=2); ax.grid(axis="y", alpha=.2)
    save(fig, figures, "figure_B_conventional_curvature")

    # C. Signed and unsigned raw edge curvature
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4), sharex=True)
    for ax, key, title in zip(axes, ("signed_raw_edge_A_rad", "unsigned_raw_edge_A_rad"),
                              ("Signed raw edge integral", "Unsigned raw edge integral")):
        for offset, method in enumerate(METHODS):
            vals = [float(group[method][key]) for group in groups]
            ax.bar(x + (offset - 1.5) * width, vals, width, label=LABELS[method], color=COLORS[method])
        ax.axhline(0, color="black", linewidth=.8); ax.set_xticks(x, SYSTEMS)
        ax.set_title(title); ax.set_ylabel(r"$\int\beta\,ds$ or $\int|\beta|\,ds$ ($\mathrm{\AA}\,\mathrm{rad}$)")
        ax.grid(axis="y", alpha=.2)
    axes[0].legend(frameon=False, fontsize=8)
    save(fig, figures, "figure_C_signed_unsigned_curvature")

    # D. AW decomposition: grouped signed components, never stacked.
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.7), sharex=True)
    for ax, key, label in zip(axes, ("smooth_H_A", "conventional_edge_H_A", "combined_H_A"),
                              (r"$\int H\,dA$", r"$\frac{1}{2}\int\beta\,ds$", r"$C_H^{AW}$")):
        vals = [float(next(r for r in static if r["system"] == system and r["method"] == "aw")[key]) for system in SYSTEMS]
        ax.bar(np.arange(3), vals, color="#D65332")
        ax.axhline(0, color="black", linewidth=.8); ax.set_xticks(np.arange(3), SYSTEMS)
        ax.set_title(label); ax.grid(axis="y", alpha=.2)
    axes[0].set_ylabel(r"Integrated mean curvature ($\mathrm{\AA}$)")
    save(fig, figures, "figure_D_aw_decomposition")

    # E. Pairwise full Power (AW radii) vs AW across static cases and frames.
    comparison = [(r, next(q for q in (static if r["kind"] == "static" else trajectory)
                            if q["filename"] == r["filename"] and q["method"] == "aw"))
                  for r in [*static, *trajectory] if r["method"] == "full_power_aw_radii"]
    if not comparison:
        raise ValueError("No saved rows for radius-matched Power/AW parity plots")
    specs = [("interface_area_A2", "Interface area (Å²)"), ("combined_H_A", r"$C_H$ (Å)"),
             ("unsigned_raw_edge_A_rad", "Unsigned edge integral (Å·rad)"),
             ("interface_atom_count", "Interface atom count")]
    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    for ax, (metric, ylabel) in zip(axes.flat, specs):
        px, py = [], []
        for power, aw in comparison:
            px.append(float(power[metric])); py.append(float(aw[metric]))
        ax.scatter(px, py, s=32, c="#4B7091", alpha=.82)
        low, high = min(px + py), max(px + py)
        pad = .04 * (high - low if high > low else max(1, abs(low)))
        ax.plot([low - pad, high + pad], [low - pad, high + pad], color="black", linestyle="--", linewidth=1)
        ax.set(xlabel="Full Power (AW radii)", ylabel="AW", title=ylabel); ax.grid(alpha=.2)
    save(fig, figures, "figure_E_power_aw_parity")

    # F. 2KAI trajectory with frozen crystal reference.
    frows = [r for r in trajectory if r["method"] in {"aw", "full_power_aw_radii"}]
    if frows:
        times = sorted({float(r["time_ps"]) for r in frows})
        keys = ("interface_area_A2", "combined_H_A", "unsigned_edge_per_area_rad_per_A")
        labels = ("Interface area (Å²)", r"Conventional $C_H$ (Å)", "Unsigned curvature / area (rad/Å)")
        fig, axes = plt.subplots(3, 1, figsize=(10.5, 10), sharex=True)
        for ax, key, ylabel in zip(axes, keys, labels):
            for method in ("full_power_aw_radii", "aw"):
                selected = sorted((r for r in frows if r["method"] == method), key=lambda r: float(r["time_ps"]))
                ax.plot([float(r["time_ps"]) for r in selected], finite(selected, key, f"trajectory {key}"),
                        marker="o", label=LABELS[method], color=COLORS[method])
                baseline = next((r for r in static if r["system"] == "2KAI" and r["method"] == method), None)
                if baseline is not None:
                    ax.axhline(float(baseline[key]), color=COLORS[method], alpha=.42, linestyle="--")
            ax.set_ylabel(ylabel); ax.grid(alpha=.2)
        axes[-1].set_xlabel("Frame time (ps)"); axes[0].legend(frameon=False)
        save(fig, figures, "figure_F_2kai_trajectory")

    # G. Paired finite-solvent effects (delta only).
    pp = [r for r in paired if r.get("interface_kind") == "protein_protein" and
          r["method"] in {"aw", "full_power_aw_radii"}]
    if pp:
        fig, axes = plt.subplots(2, 1, figsize=(10.4, 7.6), sharex=True)
        for ax, metric, ylabel in zip(axes, ("interface_area_A2", "combined_H_A"),
                                     ("Δ area (Å²)", r"Δ conventional $C_H$ (Å)")):
            subset = [r for r in pp if r["metric"] == metric]
            for method in ("full_power_aw_radii", "aw"):
                selected = sorted((r for r in subset if r["method"] == method), key=lambda r: float(r["time_ps"]))
                ax.plot([float(r["time_ps"]) for r in selected], finite(selected, "delta_solvent_HS_minus_H", f"solvent {metric}"),
                        marker="o", label=LABELS[method], color=COLORS[method])
            ax.axhline(0, color="black", linewidth=.8); ax.set_ylabel(ylabel); ax.grid(alpha=.2)
        axes[0].legend(frameon=False); axes[-1].set_xlabel("Frame time (ps)")
        save(fig, figures, "figure_G_finite_solvent_effect")

    # H. Three effect sizes in separate unit-consistent rows.
    if frows and pp:
        fig, axes = plt.subplots(3, 1, figsize=(10.4, 10), sharex=True)
        for ax, metric, ylabel in zip(axes, ("interface_area_A2", "combined_H_A", "unsigned_edge_per_area_rad_per_A"),
                                      ("Area difference (Å²)", r"$C_H$ difference (Å)", "Unsigned / area difference (rad/Å)")):
            for effect, method, style in (("geometry", "full_power_aw_radii", "-"),
                                          ("conformation", "aw", "--"), ("solvent", "aw", ":")):
                vals, tx = [], []
                for frame in sorted({r["filename"] for r in frows}):
                    indexed = {r["method"]: r for r in frows if r["filename"] == frame}
                    t = float(indexed["aw"]["time_ps"])
                    if effect == "geometry":
                        value = float(indexed["aw"][metric]) - float(indexed[method][metric])
                    elif effect == "conformation":
                        crystal = next(r for r in static if r["system"] == "2KAI" and r["method"] == "aw")
                        value = float(indexed["aw"][metric]) - float(crystal[metric])
                    else:
                        matched = [r for r in paired if r["filename"] == frame and r["method"] == "aw" and r["metric"] == metric]
                        if not matched:
                            continue
                        value = float(matched[0]["delta_solvent_HS_minus_H"])
                    tx.append(t); vals.append(value)
                ax.plot(tx, vals, linestyle=style, marker="o", label=f"{effect}: {LABELS[method]}")
            ax.axhline(0, color="black", linewidth=.8); ax.set_ylabel(ylabel); ax.grid(alpha=.2)
        axes[0].legend(frameon=False, fontsize=8, ncol=2); axes[-1].set_xlabel("Frame time (ps)")
        save(fig, figures, "figure_H_method_effect_decomposition")
    return sorted(p.name for p in figures.glob("figure_*.png"))
