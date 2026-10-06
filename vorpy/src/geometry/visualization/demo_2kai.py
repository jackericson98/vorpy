"""Build the matched-radius 2KAI Power/AW visualization checkpoint.

Example:
    python -m vorpy.src.geometry.visualization.demo_2kai \
      --aw-network output/comparative_study/cases/2KAI/H/full_aw.vpy \
      --pdb vorpy/data/2KAI.pdb \
      --phase5-dir output/comparative_study/phase5/2KAI \
      --output-dir output/2KAI/visualization
"""

from __future__ import annotations

import argparse
import ast
import csv
import json
from pathlib import Path

import numpy as np

from vorpy.src.geometry.duals import build_dual
from vorpy.src.geometry.filtrations import (
    AlphaFiltration,
    AlphaSimplex,
    build_alpha_filtration,
)
from vorpy.src.geometry.interfaces import build_alpha_interface, build_full_interface
from vorpy.src.geometry.visualization import export_dual_visualization
from vorpy.src.network.network import Network

POWER_AREA_REFERENCE = 131.679423
AW_AREA_REFERENCE = 127.679748


def _read_alpha_pairs(path, column):
    pairs = set()
    with Path(path).open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if row.get(column, "").lower() == "true":
                chain_i = row.get("atom_i_chain", "")
                chain_j = row.get("atom_j_chain", "")
                bicolor = ((chain_i in {"A", "B"} and chain_j == "I")
                           or (chain_j in {"A", "B"} and chain_i == "I"))
                if (chain_i or chain_j) and not bicolor:
                    continue
                ids = tuple(sorted(map(int, ast.literal_eval(row["generator_tuple"]))))
                pairs.add(ids)
    return pairs


def _groups_from_stable_ids(path):
    group_a, group_b = set(), set()
    with Path(path).open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            stable_id = json.loads(row["atom_id"])
            chain, generator_id = stable_id[0], int(row["aw_generator_id"])
            if chain in {"A", "B"}:
                group_a.add(generator_id)
            elif chain == "I":
                group_b.add(generator_id)
    return group_a, group_b


def _stable_metadata(path, generator_column):
    result = {}
    with Path(path).open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            chain, resnum, icode, resname, atom = json.loads(row["atom_id"])
            generator_id = int(row[generator_column])
            stable = f"{chain}|{resnum}{icode}|{resname}|{atom}"
            result[generator_id] = {
                "chain": chain, "residue_number": resnum,
                "insertion_code": icode, "residue_name": resname,
                "atom_name": atom, "stable_atom_id": stable,
            }
    return result


def _matched_power_network(aw_network):
    balls = aw_network.balls[["loc", "rad", "num", "name", "mass"]].copy(deep=True)
    settings = dict(aw_network.settings)
    settings["net_type"] = "pow"
    settings["print_metrics"] = False
    network = Network(
        locs=balls["loc"].tolist(), rads=balls["rad"].tolist(),
        names=balls["name"].tolist(), group=list(map(int, balls.index)),
        group_name="2KAI_matched_radius_power", settings=settings, balls=balls,
        sort_balls=False, system=None,
    )
    network.build()
    if list(map(int, network.balls.index)) != list(map(int, aw_network.balls.index)):
        raise RuntimeError("Power build changed network-local generator IDs relative to AW.")
    if not np.array_equal(np.asarray(network.balls["rad"], dtype=float),
                          np.asarray(aw_network.balls["rad"], dtype=float)):
        raise RuntimeError("Matched-radius Power radii differ from the AW reference network.")
    if not np.array_equal(np.asarray(network.balls["loc"].tolist(), dtype=float),
                          np.asarray(aw_network.balls["loc"].tolist(), dtype=float)):
        raise RuntimeError("Matched-radius Power coordinates differ from AW reference network.")
    return network


def _aw_pair_filtration(network, dual, pairs, alpha=0.0, tolerance=1e-6):
    """Query exact existing AW surface_birth only for persisted validated pairs."""
    from vorpy.src.geometry.aw_alpha.births import surface_birth

    records = {dimension: {} for dimension in range(4)}
    for generator_id, row in network.balls.iterrows():
        birth = -float(row["rad"])
        item = AlphaSimplex(0, (int(generator_id),), birth, birth,
                            "generator_radius_convention", True, residual=0.0)
        records[0][item.generator_tuple] = item
    for pair in sorted(pairs):
        simplex = dual.simplex(1, pair)
        if simplex is None or not simplex.supported:
            raise RuntimeError(f"Validated AW alpha pair {pair} has no supported dual edge.")
        birth, diagnostic = surface_birth(network, simplex, tolerance=tolerance)
        if birth is None or birth > alpha + tolerance:
            raise RuntimeError(f"Stored AW alpha=0 pair {pair} failed its current birth audit: {diagnostic}")
        item = AlphaSimplex(1, tuple(pair), birth, birth,
                            diagnostic.get("birth_source", "unresolved"), True,
                            residual=diagnostic.get("residual_A"),
                            notes=diagnostic.get("notes", ""), diagnostics=diagnostic)
        records[1][item.generator_tuple] = item
    return AlphaFiltration(
        network, "aw", records,
        {"filtration_kind": "experimental AW restricted-cell nerve; pair-level reference query",
         "experimental": True, "native_alpha_units": "A",
         "birth_dimensions_calculated": 1,
         "higher_dimension_limitations": "higher-dimensional births not calculated and not rendered"},
        tolerance=tolerance,
    )


def run_demo(aw_network_path, pdb_path, phase5_dir, output_dir):
    from vorpy.src.io import load_network

    aw_network_path = Path(aw_network_path).expanduser().resolve()
    pdb_path = Path(pdb_path).expanduser().resolve()
    phase5_dir = Path(phase5_dir).expanduser().resolve()
    output_dir = Path(output_dir).expanduser().resolve()
    overlap_file = phase5_dir.parent.parent / "alpha0_power_aw_pair_overlap.csv"
    if not overlap_file.exists():
        # The phase-5 root can also be supplied directly as the overlap table's parent.
        overlap_file = phase5_dir / "alpha0_power_aw_pair_overlap.csv"
    if not overlap_file.exists():
        raise FileNotFoundError(f"Validated Phase-5 alpha-pair table not found: {overlap_file}")
    aw_pairs = _read_alpha_pairs(overlap_file, "in_aw_alpha0")
    power_pairs = _read_alpha_pairs(overlap_file, "in_power_alpha0")
    if len(aw_pairs) != 27 or len(power_pairs) != 27 or aw_pairs != power_pairs:
        raise RuntimeError(f"Expected the validated 27 shared 2KAI pairs; got AW={len(aw_pairs)}, Power={len(power_pairs)}")

    aw_network = load_network(aw_network_path)
    audit_path = phase5_dir.parent.parent / "cases" / "2KAI" / "H" / "radii_audit.csv"
    group_a, group_b = _groups_from_stable_ids(audit_path)
    aw_atom_metadata = _stable_metadata(audit_path, "aw_generator_id")
    power_atom_metadata = _stable_metadata(audit_path, "power_generator_id")
    if not group_a or not group_b:
        raise RuntimeError("Could not recover the validated 2KAI A/B generator groups from radii_audit.csv")
    available = set(map(int, aw_network.balls.index))
    if any(not set(pair) <= available for pair in aw_pairs):
        raise RuntimeError("Alpha pair table contains a generator ID absent from the AW network.")

    aw_dual = build_dual(aw_network)
    aw_filtration = _aw_pair_filtration(aw_network, aw_dual, aw_pairs)
    aw_filtration.metadata["stable_atom_metadata"] = aw_atom_metadata
    aw_interface = build_alpha_interface(
        aw_network, aw_dual, aw_filtration, 0.0, group_a, group_b,
        calculate_curvature=False,
    )
    aw_full = build_full_interface(aw_network, aw_dual, group_a, group_b)
    if aw_interface.area_A2 is None or abs(aw_interface.area_A2 - AW_AREA_REFERENCE) > 2e-5:
        raise RuntimeError(f"AW area regression failed: {aw_interface.area_A2}")
    aw_area = aw_interface.area_A2
    aw_output = Path(output_dir)
    aw_metadata = export_dual_visualization(
        aw_network, aw_dual, aw_output, filtration=aw_filtration, alpha=0.0,
        max_dimension=1, interface_full=aw_full, interface_alpha=aw_interface,
        structure_path=pdb_path,
    )
    del aw_full, aw_interface, aw_filtration, aw_dual

    power_network = _matched_power_network(aw_network)
    power_dual = build_dual(power_network)
    power_filtration = build_alpha_filtration(power_network, max_dimension=3)
    power_filtration.metadata["stable_atom_metadata"] = power_atom_metadata
    selected_power = {tuple(row.generator_tuple) for row in power_filtration.interface_at(0.0, group_a, group_b)}
    if selected_power != power_pairs:
        raise RuntimeError(f"Power alpha=0 pair identities changed: expected {len(power_pairs)}, got {len(selected_power)}")
    power_interface = build_alpha_interface(
        power_network, power_dual, power_filtration, 0.0, group_a, group_b,
        calculate_curvature=False,
    )
    power_full = build_full_interface(power_network, power_dual, group_a, group_b)
    if power_interface.area_A2 is None or abs(power_interface.area_A2 - POWER_AREA_REFERENCE) > 2e-5:
        raise RuntimeError(f"Power area regression failed: {power_interface.area_A2}")
    power_output = Path(output_dir)
    power_metadata = export_dual_visualization(
        power_network, power_dual, power_output, filtration=power_filtration,
        alpha=0.0, interface_full=power_full, interface_alpha=power_interface,
        structure_path=pdb_path,
    )
    result = {
        "system": "2KAI",
        "atom_model": "matched protein heavy atoms",
        "radius_control": "Power radii copied exactly from the validated AW-loaded network",
        "alpha_native": {"power": {"value": 0.0, "units": "A^2"},
                         "aw": {"value": 0.0, "units": "A"}},
        "expected_shared_bicolor_pairs": len(power_pairs & aw_pairs),
        "power_area_A2": power_interface.area_A2,
        "aw_area_A2": aw_area,
        "power_visualization": power_metadata,
        "aw_visualization": aw_metadata,
        "aw_note": "Only validated AW pair-level births are queried; AW triangle/tetrahedron births are not calculated or rendered.",
    }
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "2KAI_visualization_summary.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aw-network", required=True, type=Path)
    parser.add_argument("--pdb", required=True, type=Path)
    parser.add_argument("--phase5-dir", required=True, type=Path,
                        help="2KAI phase5 directory containing power/aw mappings")
    parser.add_argument("--output-dir", type=Path, default=Path("output/2KAI/visualization"))
    args = parser.parse_args(argv)
    result = run_demo(args.aw_network, args.pdb, args.phase5_dir, args.output_dir)
    print(json.dumps({key: result[key] for key in (
        "system", "expected_shared_bicolor_pairs", "power_area_A2", "aw_area_A2")}, indent=2))
    print(f"Visualization exports written to: {args.output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
