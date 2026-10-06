"""Audits distinguishing generic matched-radius Power alpha from Cazals alpha."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


def write_power_alpha0_selection_comparison(
    network, pdb_path, group_a_chains, group_b_chains, output_dir, *,
    matched_group_a=None, matched_group_b=None, system_name=None,
    stable_atom_ids=None,
):
    """Compare two unchanged Power alpha-zero selections and audit their inputs.

    The matched-radius branch uses the radii on ``network.balls``. The frozen
    Cazals branch uses its PDB atom model and probe-expanded Chothia radii.
    This function reports the differences; it does not alter either input.
    """
    from vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation.cazals_pdb import (
        prepare_cazals_pdb,
    )
    from vorpy.src.analyze.power_interface_curvature import (
        analyze_power_interface_curvature,
    )

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    ids = [int(value) for value in network.balls.index]
    stable_atom_ids = stable_atom_ids or {}
    points = np.asarray(network.balls.loc[ids, "loc"].tolist(), dtype=float)
    radii = np.asarray(network.balls.loc[ids, "rad"].tolist(), dtype=float)
    if matched_group_a is None or matched_group_b is None:
        chains = [_field(network.balls.loc[index], ("chain_name", "chain")) for index in ids]
        group_a = {index for index, chain in zip(ids, chains, strict=True)
                   if chain in set(group_a_chains)}
        group_b = {index for index, chain in zip(ids, chains, strict=True)
                   if chain in set(group_b_chains)}
    else:
        group_a, group_b = set(map(int, matched_group_a)), set(map(int, matched_group_b))
    positions = {generator_id: position for position, generator_id in enumerate(ids)}
    matched = analyze_power_interface_curvature(
        points, radii,
        {positions[value] for value in group_a},
        {positions[value] for value in group_b},
        alpha=0.0, condition_beta_m=1e300,
    )
    matched_pairs = {
        tuple(sorted((ids[facet.generator_ids[0]], ids[facet.generator_ids[1]])))
        for facet in matched.alpha_zero_facets
        if _bicolor(tuple(positions[value] for value in facet.generator_ids),
                    {positions[value] for value in group_a},
                    {positions[value] for value in group_b})
    }

    cazals = prepare_cazals_pdb(
        pdb_path, group_a_chains, group_b_chains, include_hydrogens=False,
        include_hetatm=False, strict_heavy_fallbacks=False, verbose=False,
    )
    frozen = analyze_power_interface_curvature(
        cazals.points, cazals.expanded_radii, cazals.group_a, cazals.group_b,
        alpha=0.0, condition_beta_m=5.0,
        fallback_radius_assignments=len(cazals.fallback_atoms),
    )
    frozen_pairs = {
        tuple(sorted(facet.generator_ids)) for facet in frozen.alpha_zero_facets
        if _bicolor(facet.generator_ids, set(map(int, cazals.group_a)),
                    set(map(int, cazals.group_b)))
    }

    network_atoms = {
        _normalize_key(stable_atom_ids[index]) if index in stable_atom_ids else _stable_key(network.balls.loc[index]): index
        for index in ids
    }
    pdb_atoms = {_cazals_key(atom): atom for atom in cazals.atoms}
    common = set(network_atoms) & set(pdb_atoms)
    network_only, pdb_only = set(network_atoms) - common, set(pdb_atoms) - common
    atom_rows = []
    for key in sorted(set(network_atoms) | set(pdb_atoms)):
        network_id = network_atoms.get(key)
        atom = pdb_atoms.get(key)
        xyz_network = (np.asarray(network.balls.at[network_id, "loc"], dtype=float)
                       if network_id is not None else None)
        xyz_pdb = (np.asarray((atom.x, atom.y, atom.z), dtype=float)
                   if atom is not None else None)
        xyz_equal = (bool(np.allclose(xyz_network, xyz_pdb, rtol=0.0, atol=1e-8))
                     if xyz_network is not None and xyz_pdb is not None else False)
        aw_radius = float(network.balls.at[network_id, "rad"]) if network_id is not None else None
        base_radius = atom.base_radius if atom is not None else None
        expanded = atom.expanded_radius if atom is not None else None
        atom_rows.append({
            "record_type": "atom_input_audit", "stable_atom_id": "|".join(map(str, key)),
            "network_generator_id": network_id,
            "present_in_both": key in common,
            "coordinates_equal_at_1e-8_A": xyz_equal,
            "network_x_A": xyz_network[0] if xyz_network is not None else None,
            "network_y_A": xyz_network[1] if xyz_network is not None else None,
            "network_z_A": xyz_network[2] if xyz_network is not None else None,
            "cazals_x_A": xyz_pdb[0] if xyz_pdb is not None else None,
            "cazals_y_A": xyz_pdb[1] if xyz_pdb is not None else None,
            "cazals_z_A": xyz_pdb[2] if xyz_pdb is not None else None,
            "network_matched_radius_A": aw_radius,
            "cazals_base_radius_A": base_radius,
            "cazals_probe_expansion_A": cazals.probe_radius if atom is not None else None,
            "cazals_expanded_radius_A": expanded,
            "network_power_weight_A2": aw_radius**2 if aw_radius is not None else None,
            "cazals_power_weight_A2": atom.weight if atom is not None else None,
            "network_radius_class": "AW-loaded network radius" if network_id is not None else None,
            "cazals_radius_class": atom.radius_class if atom is not None else None,
        })

    cazals_key_to_position = {_cazals_key(atom): atom.index for atom in cazals.atoms}
    pair_rows = []
    for pair in sorted(matched_pairs | frozen_pairs):
        key_i = (_normalize_key(stable_atom_ids[pair[0]]) if pair[0] in stable_atom_ids
                 else _stable_key(network.balls.loc[pair[0]]) if pair[0] in network.balls.index else None)
        key_j = (_normalize_key(stable_atom_ids[pair[1]]) if pair[1] in stable_atom_ids
                 else _stable_key(network.balls.loc[pair[1]]) if pair[1] in network.balls.index else None)
        if key_i in cazals_key_to_position and key_j in cazals_key_to_position:
            ci, cj = cazals_key_to_position[key_i], cazals_key_to_position[key_j]
            stable_pair = tuple(sorted((ci, cj)))
        else:
            stable_pair = None
        pair_rows.append({
            "record_type": "pair_selection", "generator_i": pair[0], "generator_j": pair[1],
            "stable_atom_i": "|".join(map(str, key_i)) if key_i else None,
            "stable_atom_j": "|".join(map(str, key_j)) if key_j else None,
            "in_matched_radius_power_alpha0": pair in matched_pairs,
            "in_cazals_chothia_probe_power_alpha0": stable_pair in frozen_pairs if stable_pair else False,
            "pair_mapping_to_cazals": stable_pair is not None,
        })

    same_atoms = not network_only and not pdb_only and len(common) == len(ids) == len(cazals.atoms)
    same_coordinates = same_atoms and all(row["coordinates_equal_at_1e-8_A"] for row in atom_rows)
    radii_differ = any(
        row["present_in_both"] and not np.isclose(
            row["network_matched_radius_A"], row["cazals_expanded_radius_A"],
            rtol=0.0, atol=1e-12,
        ) for row in atom_rows
    )
    summary = {
        "system": str(system_name or getattr(network, "group_name", Path(pdb_path).stem)),
        "network_atom_count": len(ids), "cazals_atom_count": len(cazals.atoms),
        "matched_atom_count": len(common), "network_only_atoms": len(network_only),
        "cazals_only_atoms": len(pdb_only), "same_atom_set": same_atoms,
        "same_coordinates_at_1e-8_A": same_coordinates,
        "network_radius_model": "exact radius supplied to the AW-loaded network; used directly as Power radius",
        "cazals_radius_model": "Chothia/Cazals base radii plus 1.4 A probe expansion",
        "different_radius_vectors": radii_differ,
        "network_power_weight_rule": "w_i = r_i^2",
        "cazals_power_weight_rule": "w_i = (r_i + 1.4 A)^2",
        "matched_radius_alpha_convention": "weighted squared power distance threshold, alpha=0 A^2",
        "cazals_alpha_convention": "same weighted squared power distance threshold, alpha=0 A^2",
        "matched_bicolor_rule": "all bicolor alpha-zero dual edges; no M filter",
        "cazals_bicolor_rule": "alpha-zero bicolor facets before M=5; M=5 final count also reported",
        "matched_alpha0_pairs": len(matched_pairs),
        "cazals_alpha0_pairs": len(frozen_pairs),
        "cazals_alpha0_m5_final_pairs": len(frozen.final_facets),
        "shared_pair_count": len(matched_pairs) - len(
            matched_pairs - _map_cazals_pairs_to_network(frozen_pairs, cazals, network_atoms)
        ),
        "matched_only_pairs": len(matched_pairs - _map_cazals_pairs_to_network(frozen_pairs, cazals, network_atoms)),
        "cazals_only_pairs": len(frozen_pairs - _map_network_pairs_to_cazals(matched_pairs, network_atoms, cazals)),
        "coordinate_or_atom_model_difference": not same_atoms or not same_coordinates,
        "selection_divergence_cause": (
            "Power weights differ because the matched branch uses AW-loaded radii directly, while frozen Cazals uses atom-specific Chothia/Cazals base radii plus the 1.4 A probe; atoms and coordinates are identical. The alpha threshold and bicolor rule are otherwise the same before M=5."
            if same_atoms and same_coordinates and radii_differ else
            "Multiple input differences are present; inspect per-atom audit rows before attributing selection changes."
        ),
    }
    _write_rows(output / "power_alpha0_selection_comparison.csv", atom_rows + pair_rows)
    (output / "power_alpha0_selection_comparison_summary.txt").write_text(
        "Power alpha=0 selection input audit\n"
        + "\n".join(f"{key}: {value}" for key, value in summary.items()) + "\n",
        encoding="utf-8",
    )
    return summary


def _stable_key(row):
    return (
        _field(row, ("chain_name", "chain")),
        _field(row, ("auth_seq_id", "res_seq", "residue_number")),
        _field(row, ("pdb_ins_code", "insertion_code", "icode")),
        _field(row, ("auth_comp_id", "res_name", "residue_name")),
        _field(row, ("auth_atom_id", "name", "atom_name")),
    )


def _cazals_key(atom):
    return (atom.chain, str(atom.residue_number), atom.insertion_code,
            atom.residue_name, atom.atom_name)


def _normalize_key(values):
    return tuple("" if value is None else str(value).strip() for value in values)


def _field(row, names):
    for name in names:
        value = row.get(name)
        if value is not None and str(value).strip() not in {"", "nan", "None"}:
            return str(value).strip()
    return ""


def _bicolor(pair, group_a, group_b):
    return (pair[0] in group_a and pair[1] in group_b) or (pair[1] in group_a and pair[0] in group_b)


def _map_cazals_pairs_to_network(cazals_pairs, prepared, network_atoms):
    reverse = {atom.index: _normalize_key(_cazals_key(atom)) for atom in prepared.atoms}
    return {tuple(sorted(network_atoms[reverse[index]] for index in pair))
            for pair in cazals_pairs if all(reverse[index] in network_atoms for index in pair)}


def _map_network_pairs_to_cazals(network_pairs, network_atoms, prepared):
    reverse = {_normalize_key(_cazals_key(atom)): atom.index for atom in prepared.atoms}
    key_by_network_id = {generator_id: stable_key
                         for stable_key, generator_id in network_atoms.items()}
    return {tuple(sorted(reverse[key_by_network_id[generator_id]] for generator_id in pair))
            for pair in network_pairs
            if all(generator_id in key_by_network_id
                   and key_by_network_id[generator_id] in reverse
                   for generator_id in pair)}


def _write_rows(path, rows):
    rows = list(rows)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        if fields:
            writer.writeheader()
            writer.writerows(rows)
