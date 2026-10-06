"""Experimental AW-alpha CLI option and solved-network runner."""

from __future__ import annotations

import csv
import re
import tempfile
from pathlib import Path


def extract_aw_alpha_options(arguments):
    """Remove ``--aw-alpha-experimental [A ...]`` from legacy CLI arguments."""
    remaining, values = [], []
    args = iter(arguments)
    for arg in args:
        if arg != "--aw-alpha-experimental":
            remaining.append(arg)
            continue
        try:
            value = float(next(args))
        except (StopIteration, TypeError, ValueError):
            raise SystemExit(
                "--aw-alpha-experimental requires a finite alpha in Angstroms."
            ) from None
        if not (float("-inf") < value < float("inf")):
            raise SystemExit(
                "--aw-alpha-experimental requires a finite alpha in Angstroms."
            )
        values.append(value)
    return remaining, values


def extract_aw_alpha_pair_only(arguments):
    flag = "--aw-alpha-pairs-only-experimental"
    values = list(arguments)
    if values.count(flag) > 1:
        raise SystemExit(f"{flag} may only be specified once.")
    return [value for value in values if value != flag], flag in values


def run_aw_alpha_experimental(
    networks, destination, alpha_values=(0.0,), tolerance=1e-6, pair_overlap_only=False
):
    """Export experimental filtration data for already-solved AW networks."""
    from vorpy.src.geometry.aw_alpha import AWAlphaFiltration
    from vorpy.src.geometry.aw_alpha.export import export_aw_alpha
    from vorpy.src.geometry.aw_alpha.power_comparison import compare_alpha0_pairs

    destination = Path(destination)
    results = []
    comparison_rows, comparison_summaries = [], []
    eligible = [
        (name, net)
        for name, net in networks
        if (getattr(net, "settings", None) or {}).get("net_type", "aw") == "aw"
    ]
    if not eligible:
        raise ValueError(
            "No solved AW network is available for experimental AW alpha analysis."
        )
    for name, network in eligible:
        slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(name)).strip("._") or "network"
        filtration = AWAlphaFiltration(network, name=name, tolerance=tolerance)
        groups = getattr(network, "iface_grps", None)
        a, b = groups if groups is not None else (None, None)
        if pair_overlap_only:
            if groups is None:
                raise ValueError(
                    "Pair-overlap-only mode requires a two-group interface network."
                )
            filtration.calculate_births(max_dimension=1)
            result = {
                "unresolved": [
                    record
                    for record in filtration.records[1].values()
                    if record.geometric_birth is None
                ]
            }
        else:
            result = export_aw_alpha(
                filtration,
                destination / slug,
                alpha_values=alpha_values,
                group_a=a,
                group_b=b,
                radius_convention="radii supplied to the existing solved AW network",
            )
        results.append((name, result))
        groups = getattr(network, "iface_grps", None)
        if groups is not None and "water" not in str(name).lower():
            try:
                pairs, summary = compare_alpha0_pairs(
                    network,
                    groups[0],
                    groups[1],
                    name=name,
                    tolerance=tolerance,
                    filtration=filtration,
                )
            except (ImportError, AttributeError, RuntimeError, ValueError) as error:
                summary = {"system": str(name), "status": "failed", "error": str(error)}
                pairs = []
            else:
                summary["status"] = "complete"
            comparison_rows.extend(pairs)
            comparison_summaries.append(summary)
    if comparison_summaries:
        destination.mkdir(parents=True, exist_ok=True)
        _write_rows(destination / "alpha0_power_aw_pair_overlap.csv", comparison_rows)
        _write_rows(
            destination / "alpha0_power_aw_pair_summary.csv", comparison_summaries
        )
    return results


def _write_rows(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parser_safe_aw_alpha_structure(structure):
    """Relabel solvent records to the reserved SOL chain in a temporary PDB.

    The legacy PDB reader assumes solvent chain IDs do not collide with protein
    chains. Some crystallographic files violate that assumption. This preserves
    every record and coordinate while routing solvent into its existing solvent
    container; the original file is never changed.
    """
    from vorpy.src.boundary import SOLVENT_RESIDUES

    source = Path(structure)
    lines = source.read_text(encoding="utf-8", errors="replace").splitlines(
        keepends=True
    )
    normalized, count = [], 0
    for line in lines:
        if (
            line.startswith("HETATM")
            and len(line) >= 22
            and line[17:20].strip().upper() in SOLVENT_RESIDUES
            and line[21] != " "
        ):
            line = line[:21] + " " + line[22:]
            count += 1
        normalized.append(line)
    if not count:
        return None, 0
    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".pdb",
        prefix="vorpy_aw_alpha_",
        encoding="utf-8",
        delete=False,
    ) as handle:
        handle.writelines(normalized)
        temporary = Path(handle.name)
    return temporary, count
