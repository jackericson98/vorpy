"""Compact validation suite for the weighted Power/Cazals interface analysis.

This module is an analysis/reporting layer. It leaves the AW implementation,
Power solver, PDB inputs, and VorPy archive format untouched.
"""

from __future__ import annotations

from collections import defaultdict, deque
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from vorpy.src.analyze.power_interface_curvature import (
    analyze_power_interface_curvature,
    analyze_power_regular_complex,
    cazals_signed_beta,
)
from vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation.cazals_pdb import (
    prepare_cazals_pdb,
)


AMINO_ACIDS = frozenset("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
WATER_NAMES = frozenset({"HOH", "WAT", "DOD", "H2O"})
DEFAULT_PROBE = 1.4
DEFAULT_M = 5.0
SCC_AREA_FRACTION = 0.075

# Heavy-atom covalent bonds in the standard amino-acid residue templates.
_SIDE_BONDS = {
    "ALA": [("CA", "CB")],
    "ARG": [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "NE"), ("NE", "CZ"), ("CZ", "NH1"), ("CZ", "NH2")],
    "ASN": [("CA", "CB"), ("CB", "CG"), ("CG", "OD1"), ("CG", "ND2")],
    "ASP": [("CA", "CB"), ("CB", "CG"), ("CG", "OD1"), ("CG", "OD2")],
    "CYS": [("CA", "CB"), ("CB", "SG")],
    "GLN": [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "NE2")],
    "GLU": [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "OE2")],
    "GLY": [],
    "HIS": [("CA", "CB"), ("CB", "CG"), ("CG", "ND1"), ("CG", "CD2"), ("ND1", "CE1"), ("CE1", "NE2"), ("NE2", "CD2")],
    "ILE": [("CA", "CB"), ("CB", "CG1"), ("CB", "CG2"), ("CG1", "CD1")],
    "LEU": [("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2")],
    "LYS": [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "CE"), ("CE", "NZ")],
    "MET": [("CA", "CB"), ("CB", "CG"), ("CG", "SD"), ("SD", "CE")],
    "PHE": [("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CD1", "CE1"), ("CD2", "CE2"), ("CE1", "CZ"), ("CE2", "CZ")],
    "PRO": [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "N")],
    "SER": [("CA", "CB"), ("CB", "OG")],
    "THR": [("CA", "CB"), ("CB", "OG1"), ("CB", "CG2")],
    "TRP": [("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CD1", "NE1"), ("NE1", "CE2"), ("CE2", "CD2"), ("CD2", "CE3"), ("CE3", "CZ3"), ("CZ3", "CH2"), ("CH2", "CZ2"), ("CZ2", "CE2")],
    "TYR": [("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CD1", "CE1"), ("CD2", "CE2"), ("CE1", "CZ"), ("CE2", "CZ"), ("CZ", "OH")],
    "VAL": [("CA", "CB"), ("CB", "CG1"), ("CB", "CG2")],
}


def pdb_inventory(path, selected_chains):
    """Count parsed PDB atom categories without changing or rewriting input."""
    path = Path(path)
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    records = [line for line in lines if line.startswith(("ATOM  ", "HETATM"))]
    atom_records = [line for line in records if line.startswith("ATOM  ")]
    hetatm = [line for line in records if line.startswith("HETATM")]
    protein = [line for line in atom_records if line[17:20].strip().upper() in AMINO_ACIDS]
    water = [line for line in hetatm if line[17:20].strip().upper() in WATER_NAMES]
    other_hetatm = [line for line in hetatm if line not in water]
    in_chains = [line for line in protein if line[21:22].strip() in set(selected_chains)]
    return {
        "path": str(path),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "total_parsed_atom_records": len(records),
        "protein_atom_records": len(protein),
        "water_atom_records": len(water),
        "other_hetatm_records": len(other_hetatm),
        "nonprotein_atom_records": len(atom_records) - len(protein),
        "protein_atoms_in_selected_chains": len(in_chains),
        "retained_atoms": len(in_chains),
        "excluded_protein_atoms_other_chains": len(protein) - len(in_chains),
        "excluded_water_atoms": len(water),
        "excluded_other_hetatm_atoms": len(other_hetatm),
        "selected_chain_atom_counts": {
            chain: sum(line[21:22].strip() == chain for line in in_chains)
            for chain in sorted(set(selected_chains))
        },
    }


def compnd_chain_mapping(path):
    """Return molecule/chain text from PDB COMPND records for identity checks."""
    sections = []
    current = None
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("COMPND"):
            for field in line[10:].strip().split(";"):
                label, sep, value = field.partition(":")
                if not sep:
                    continue
                label = label.strip().upper()
                if label == "MOL_ID":
                    if current:
                        sections.append(current)
                    current = {"molecule_id": value.strip(), "molecule": "", "chains": []}
                elif current and label == "MOLECULE":
                    current["molecule"] = value.strip()
                elif current and label == "CHAIN":
                    current["chains"].extend(piece.strip() for piece in value.split(",") if piece.strip())
    if current:
        sections.append(current)
    return sections


def covalent_heavy_bonds(prepared, pdb_path):
    """Build known protein heavy-atom bonds without a generic distance cutoff.

    Intramolecular edges come from standard amino-acid templates. Peptide C-N
    edges are added only between consecutive residues in the PDB chain order,
    and a 1.9 A check verifies that known linkage. SSBOND and CONECT records
    contribute explicit cross-residue/cross-chain edges.
    """
    atoms = prepared.atoms
    points = prepared.points
    def atom_key(atom):
        return (atom.chain, atom.residue_name, atom.residue_number,
                atom.insertion_code, atom.atom_name.strip().upper())
    key_to_id = {atom_key(atom): atom.index for atom in atoms}
    serial_to_id = {atom.pdb_serial: atom.index for atom in atoms}
    bonds = set()

    def add_keys(left, right):
        i, j = key_to_id.get(left), key_to_id.get(right)
        if i is not None and j is not None and i != j:
            bonds.add(tuple(sorted((i, j))))

    for atom in atoms:
        residue = (atom.chain, atom.residue_name, atom.residue_number, atom.insertion_code)
        backbone = [("N", "CA"), ("CA", "C"), ("C", "O"), ("C", "OXT")]
        for name1, name2 in backbone + _SIDE_BONDS.get(atom.residue_name, []):
            add_keys((*residue, name1), (*residue, name2))

    chain_residues = defaultdict(list)
    for atom in atoms:
        residue = (atom.chain, atom.residue_name, atom.residue_number, atom.insertion_code)
        if not chain_residues[atom.chain] or chain_residues[atom.chain][-1] != residue:
            chain_residues[atom.chain].append(residue)
    for residues in chain_residues.values():
        for left, right in zip(residues, residues[1:]):
            c_id = key_to_id.get((*left, "C"))
            n_id = key_to_id.get((*right, "N"))
            if c_id is not None and n_id is not None and np.linalg.norm(points[c_id] - points[n_id]) < 1.9:
                bonds.add(tuple(sorted((c_id, n_id))))

    for line in Path(pdb_path).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("CONECT"):
            serials = []
            for offset in range(6, len(line), 5):
                try:
                    serials.append(int(line[offset:offset + 5]))
                except ValueError:
                    pass
            for serial in serials[1:]:
                i, j = serial_to_id.get(serials[0]), serial_to_id.get(serial)
                if i is not None and j is not None and i != j:
                    bonds.add(tuple(sorted((i, j))))
        elif line.startswith("SSBOND"):
            try:
                add_keys((line[15].strip(), "CYS", int(line[17:21]), line[21].strip(), "SG"),
                         (line[29].strip(), "CYS", int(line[31:35]), line[35].strip(), "SG"))
            except (ValueError, IndexError):
                continue
    return bonds


def power_facet_polygon_area(result, facet, *, coordinate_tolerance=1e-7,
                             triangle_tetrahedra=None):
    """Recover a bounded dual power-facet polygon from its incident tetrahedra."""
    if triangle_tetrahedra is None:
        tri_to_tets = defaultdict(list)
        for tet in result.full_simplices[3]:
            for tri in _faces(tet, 3):
                tri_to_tets[tri].append(tet)
    else:
        tri_to_tets = triangle_tetrahedra
    incident_tets = set()
    for triangle in facet.incident_regular_triangles:
        incident = tri_to_tets.get(tuple(triangle), ())
        if len(incident) != 2:
            return {"status": "unbounded_or_open_facet", "area_A2": None,
                    "polygon_vertex_count": 0, "reason": f"facet boundary triangle has {len(incident)} incident tetrahedra"}
        incident_tets.update(incident)
    missing = sorted(tet for tet in incident_tets if tet not in result.power_vertices)
    if missing:
        return {"status": "unresolved_power_vertices", "area_A2": None,
                "polygon_vertex_count": 0, "reason": f"missing {len(missing)} power vertices"}
    positions = []
    for tet in sorted(incident_tets):
        point = np.asarray(result.power_vertices[tet].position, dtype=float)
        if not any(np.linalg.norm(point - previous) <= coordinate_tolerance for previous in positions):
            positions.append(point)
    if len(positions) < 3:
        return {"status": "degenerate_polygon", "area_A2": 0.0,
                "polygon_vertex_count": len(positions), "polygon_vertices": tuple(),
                "reason": "fewer than three unique power vertices"}

    i, j = facet.generator_ids
    normal = result.points[j] - result.points[i]
    norm = float(np.linalg.norm(normal))
    if norm <= coordinate_tolerance:
        return {"status": "degenerate_bisector", "area_A2": None,
                "polygon_vertex_count": len(positions), "polygon_vertices": tuple(),
                "reason": "coincident generator centers"}
    normal /= norm
    center = np.mean(positions, axis=0)
    axis = np.eye(3)[int(np.argmin(np.abs(normal)))]
    u = np.cross(normal, axis)
    u /= np.linalg.norm(u)
    v = np.cross(normal, u)
    angles = [math.atan2(float(np.dot(point - center, v)), float(np.dot(point - center, u)))
              for point in positions]
    polygon = [point for _, point in sorted(zip(angles, positions), key=lambda item: item[0])]
    area_vector = np.zeros(3)
    for p, q in zip(polygon, polygon[1:] + polygon[:1]):
        area_vector += np.cross(p, q)
    area = abs(float(np.dot(area_vector, normal))) / 2.0
    if not np.isfinite(area) or area <= 0:
        return {"status": "degenerate_polygon", "area_A2": area,
                "polygon_vertex_count": len(polygon),
                "polygon_vertices": tuple(tuple(map(float, p)) for p in polygon),
                "reason": "nonpositive polygon area"}
    return {"status": "finite_bounded", "area_A2": area,
            "polygon_vertex_count": len(polygon),
            "polygon_vertices": tuple(tuple(map(float, p)) for p in polygon), "reason": ""}


def _faces(simplex, size):
    from itertools import combinations
    return [tuple(sorted(face)) for face in combinations(simplex, size)]


def interface_components(result, area_by_facet):
    """Connected components of selected facets under shared finite interface edges."""
    facets = {facet.generator_ids for facet in result.final_facets}
    adjacency = {facet: set() for facet in facets}
    for edge in result.included_edges:
        f1, f2 = edge.bicolor_facets
        if f1 in facets and f2 in facets:
            adjacency[f1].add(f2)
            adjacency[f2].add(f1)
    components = []
    unseen = set(facets)
    while unseen:
        seed = min(unseen)
        queue = deque([seed])
        unseen.remove(seed)
        component = set()
        while queue:
            node = queue.popleft()
            component.add(node)
            for neighbor in adjacency[node] & unseen:
                unseen.remove(neighbor)
                queue.append(neighbor)
        area = sum(area_by_facet.get(facet, 0.0) for facet in component)
        components.append({"facets": tuple(sorted(component)), "facet_count": len(component), "area_A2": area})
    total_area = sum(component["area_A2"] for component in components)
    for index, component in enumerate(sorted(components, key=lambda item: item["area_A2"], reverse=True), 1):
        component["component_id"] = index
        component["via_fraction"] = component["area_A2"] / total_area if total_area else 0.0
        component["significant"] = component["via_fraction"] >= SCC_AREA_FRACTION
    return sorted(components, key=lambda item: item["component_id"])


def interface_atom_ids(final_facets, group_a_ids):
    """Unique A/B atom sets participating in the selected bicolor facets."""
    group_a = set(group_a_ids)
    atoms_a, atoms_b = set(), set()
    for facet in final_facets:
        for atom_id in facet.generator_ids:
            (atoms_a if atom_id in group_a else atoms_b).add(atom_id)
    return atoms_a, atoms_b


def interface_atom_asymmetry(atoms_a, atoms_b):
    na, nb = len(atoms_a), len(atoms_b)
    return max(na / nb, nb / na) if na and nb else float("nan")


def facet_area_statistics(area_values):
    """Compute the requested area diagnostics from finite bounded facets."""
    values = np.asarray(tuple(area_values), dtype=float)
    if not len(values):
        return {"count": 0, "VIA_A2": 0.0, "mean_A2": None, "median_A2": None,
                "fraction_below_1_A2": None, "fraction_above_10_A2": None, "maximum_A2": None}
    return {"count": len(values), "VIA_A2": float(values.sum()), "mean_A2": float(values.mean()),
            "median_A2": float(np.median(values)),
            "fraction_below_1_A2": float(np.mean(values < 1.0)),
            "fraction_above_10_A2": float(np.mean(values > 10.0)),
            "maximum_A2": float(values.max())}


def classify_same_group_pair(pair, covalent_bonds):
    return "covalent" if tuple(sorted(map(int, pair))) in covalent_bonds else "noncovalent"


def _write_csv(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def analyze_case(pdb_path, output_dir, group_a_chains, group_b_chains,
                 *, probe_radius=DEFAULT_PROBE, condition_beta_m=DEFAULT_M,
                 regular_complex_cache=None):
    """Run the existing weighted Power analyzer and produce case validation."""
    pdb_path, output_dir = Path(pdb_path), Path(output_dir)
    selected_chains = set(group_a_chains) | set(group_b_chains)
    inventory = pdb_inventory(pdb_path, selected_chains)
    inventory["compnd_chain_mapping"] = compnd_chain_mapping(pdb_path)
    prepared = prepare_cazals_pdb(
        pdb_path, group_a_chains, group_b_chains,
        probe_radius=probe_radius, include_hydrogens=False,
        include_hetatm=False, strict_heavy_fallbacks=False, verbose=False,
    )
    inventory["retained_atoms"] = len(prepared.atoms)
    inventory["fallback_radius_assignments"] = len(prepared.fallback_atoms)
    inventory["fallback_atoms"] = [
        f"{atom.chain}:{atom.residue_name}{atom.residue_number}{atom.insertion_code}:{atom.atom_name}"
        for atom in prepared.fallback_atoms
    ]
    inventory["retained_chain_counts"] = {
        chain: sum(atom.chain == chain for atom in prepared.atoms) for chain in sorted(selected_chains)
    }
    if regular_complex_cache is None:
        result = analyze_power_interface_curvature(
            prepared.points, prepared.expanded_radii, prepared.group_a, prepared.group_b,
            alpha=0.0, condition_beta_m=condition_beta_m,
            fallback_radius_assignments=len(prepared.fallback_atoms),
        )
    else:
        full = {d: set() for d in range(4)}
        filtration = {}
        with Path(regular_complex_cache).open("r", newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                simplex = tuple(sorted(map(int, row["simplex"].replace(",", " ").split())))
                dimension = int(row["dimension"])
                if len(simplex) != dimension + 1:
                    raise ValueError(f"Malformed simplex row in {regular_complex_cache}: {row}")
                full[dimension].add(simplex)
                filtration[simplex] = float(row["filtration"])
        result = analyze_power_regular_complex(
            prepared.points, prepared.expanded_radii, prepared.group_a, prepared.group_b,
            full, filtration, alpha=0.0, condition_beta_m=condition_beta_m,
            fallback_radius_assignments=len(prepared.fallback_atoms),
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    result.export(output_dir)
    (output_dir / "input_inventory.json").write_text(json.dumps(inventory, indent=2), encoding="utf-8")

    atoms = {atom.index: atom for atom in prepared.atoms}
    facet_area = {}
    facet_rows = []
    for facet in result.final_facets:
        i, j = facet.generator_ids
        polygon = power_facet_polygon_area(result, facet)
        facet_area[facet.generator_ids] = polygon["area_A2"] or 0.0
        facet_rows.append({
            "facet_id": f"{i}-{j}", "atom_i": i, "atom_j": j,
            "atom_i_serial": atoms[i].pdb_serial, "atom_j_serial": atoms[j].pdb_serial,
            "group_i": "A" if i in result.group_a else "B",
            "group_j": "A" if j in result.group_a else "B",
            "chain_i": atoms[i].chain, "chain_j": atoms[j].chain,
            "residue_i": f"{atoms[i].residue_name}{atoms[i].residue_number}{atoms[i].insertion_code}",
            "residue_j": f"{atoms[j].residue_name}{atoms[j].residue_number}{atoms[j].insertion_code}",
            "alpha_filtration": result.filtration.get(facet.generator_ids),
            "condition_b_status": facet.condition_beta_status,
            "condition_b_m_over_r": facet.m_over_r,
            "condition_b_accepted": facet.condition_beta_accepted,
            "area_A2": polygon["area_A2"],
            "polygon_vertex_count": polygon["polygon_vertex_count"],
            "finite_bounded_status": polygon["status"], "area_reason": polygon["reason"],
        })
    _write_csv(output_dir / "power_interface_facets.csv", facet_rows)

    valid_areas = [float(row["area_A2"]) for row in facet_rows
                   if row["finite_bounded_status"] == "finite_bounded"]
    via = float(sum(valid_areas))
    area_by_facet = {key: value for key, value in facet_area.items()}
    components = interface_components(result, area_by_facet)
    component_rows = [{
        "component_id": c["component_id"], "facet_count": c["facet_count"],
        "via_A2": c["area_A2"], "via_fraction": c["via_fraction"],
        "significant_scc": c["significant"],
        "facet_ids": ";".join(f"{i}-{j}" for i, j in c["facets"]),
    } for c in components]
    _write_csv(output_dir / "power_interface_components.csv", component_rows)

    interface_a, interface_b = interface_atom_ids(result.final_facets, result.group_a)
    na, nb = len(interface_a), len(interface_b)
    r_ab = interface_atom_asymmetry(interface_a, interface_b)
    totals = result.edge_totals
    positive_length = sum(edge.length for edge in result.included_edges if edge.beta_radians > 0)
    negative_length = sum(edge.length for edge in result.included_edges if edge.beta_radians < 0)
    c_abs = totals["unsigned"]

    lys15_ids = {atom.index for atom in prepared.atoms
                 if atom.chain == "I" and atom.residue_name == "LYS" and atom.residue_number == 15}
    lys15_facets = [facet for facet in result.final_facets if lys15_ids.intersection(facet.generator_ids)]
    lys15_facet_ids = {facet.generator_ids for facet in lys15_facets}
    lys15_edges = [edge for edge in result.included_edges
                   if set(edge.bicolor_facets) & lys15_facet_ids]
    lys15_positive = [edge for edge in lys15_edges if edge.beta_radians > 0]
    lys15_negative = [edge for edge in lys15_edges if edge.beta_radians < 0]
    lys15_c = sum(edge.length_beta for edge in lys15_edges)
    lys15_l = sum(edge.length for edge in lys15_edges)
    lys15_cu = sum(edge.length_abs_beta for edge in lys15_edges)
    _write_csv(output_dir / "lys15_interface_facets.csv",
               [row for row in facet_rows if tuple(map(int, row["facet_id"].split("-"))) in lys15_facet_ids])
    lys_edges_rows = []
    for edge in lys15_edges:
        lys_edges_rows.append({
            "regular_triangle": " ".join(map(str, edge.regular_triangle)),
            "group_pattern": edge.group_composition,
            "singleton_atom": edge.singleton_generator,
            "bicolor_facets": ";".join(f"{i}-{j}" for i, j in edge.bicolor_facets),
            "length_A": edge.length, "beta_rad": edge.beta_radians,
            "beta_deg": edge.beta_degrees, "length_beta_A_rad": edge.length_beta,
            "length_abs_beta_A_rad": edge.length_abs_beta,
        })
    _write_csv(output_dir / "lys15_curvature_edges.csv", lys_edges_rows)

    # Cazals Figure 6B population, directly re-evaluated from the triangle angle.
    bonds = covalent_heavy_bonds(prepared, pdb_path)
    distribution_rows = []
    for edge_id, edge in enumerate(result.included_edges, 1):
        tri = edge.regular_triangle
        singleton, group, beta = cazals_signed_beta(tri, prepared.points, result.group_a, result.group_b)
        pair = tuple(sorted(atom_id for atom_id in tri if atom_id != singleton))
        a0, a1, a2 = (atoms[atom_id] for atom_id in tri)
        pair_atoms = [atoms[atom_id] for atom_id in pair]
        distribution_rows.append({
            "edge_id": edge_id, "atom_ids": " ".join(map(str, tri)),
            "atom_serials": " ".join(str(atom.pdb_serial) for atom in (a0, a1, a2)),
            "residue_ids": ";".join(f"{atom.chain}:{atom.residue_name}{atom.residue_number}{atom.insertion_code}" for atom in (a0, a1, a2)),
            "group_pattern": edge.group_composition,
            "singleton_atom_id": singleton, "singleton_serial": atoms[singleton].pdb_serial,
            "singleton_residue_id": f"{atoms[singleton].chain}:{atoms[singleton].residue_name}{atoms[singleton].residue_number}{atoms[singleton].insertion_code}:{atoms[singleton].atom_name}",
            "same_group_atom_pair_ids": " ".join(map(str, pair)),
            "same_group_atom_pair_serials": " ".join(str(atom.pdb_serial) for atom in pair_atoms),
            "same_group_residue_ids": ";".join(f"{atom.chain}:{atom.residue_name}{atom.residue_number}{atom.insertion_code}:{atom.atom_name}" for atom in pair_atoms),
            "same_group_pair_bond_status": classify_same_group_pair(pair, bonds),
            "edge_length_A": edge.length, "beta_signed_deg": math.degrees(beta),
            "abs_beta_deg": abs(math.degrees(beta)), "l_beta_A_rad": edge.length * beta,
            "beta_current_deg": edge.beta_degrees,
        })
    _write_csv(output_dir / f"{pdb_path.stem.upper()}_beta_distribution.csv", distribution_rows)
    histogram = _histogram(distribution_rows)
    _write_csv(output_dir / f"{pdb_path.stem.upper()}_beta_histogram.csv", histogram)
    _plot_histogram(distribution_rows, output_dir / f"{pdb_path.stem.upper()}_beta_distribution.png",
                    title=f"{pdb_path.stem.upper()} Power-interface singleton angles")

    edge_rows = [{
        "triangle": " ".join(map(str, edge.regular_triangle)),
        "group_pattern": edge.group_composition,
        "facet_1": " ".join(map(str, edge.bicolor_facets[0])),
        "facet_2": " ".join(map(str, edge.bicolor_facets[1])),
        "length_A": edge.length, "beta_rad": edge.beta_radians,
        "beta_deg": edge.beta_degrees, "length_beta_A_rad": edge.length_beta,
        "status": edge.status,
    } for edge in result.edges]
    _write_csv(output_dir / "power_interface_edge_incidence.csv", edge_rows)

    beta_differences = np.asarray([
        float(row["beta_current_deg"]) - float(row["beta_signed_deg"])
        for row in distribution_rows
    ])
    direct_matches = np.abs(beta_differences) < 1e-8
    supplementary_matches = np.asarray([
        abs(float(row["beta_current_deg"]) - math.copysign(180.0 - abs(float(row["beta_signed_deg"])),
                                                           float(row["beta_signed_deg"]))) < 1e-8
        for row in distribution_rows
    ])
    sign_disagreements = sum(float(row["beta_current_deg"]) * float(row["beta_signed_deg"]) < 0
                             for row in distribution_rows)

    summary = {
        "pdb": str(pdb_path), "inventory": inventory,
        "groups": {"A_chains": sorted(group_a_chains), "B_chains": sorted(group_b_chains)},
        "settings": {"probe_radius_A": probe_radius, "alpha": 0.0, "M": condition_beta_m},
        "regular_counts": result.full_counts, "alpha_zero_counts": result.alpha_counts,
        "full_regular_bicolor_facets": len(result.all_bicolor_facets),
        "alpha_zero_bicolor_facets": len(result.alpha_zero_facets),
        "condition_b_rejected": sum(f.alpha_zero_selected and f.condition_beta_accepted is False for f in result.facets),
        "condition_b_unresolved": sum(f.alpha_zero_selected and f.condition_beta_accepted is None for f in result.facets),
        "final_facets": len(result.final_facets),
        "interface_atoms_A": na, "interface_atoms_B": nb, "interface_atoms_total": na + nb, "r_AB": r_ab,
        "facet_area": {"bounded_facets": len(valid_areas), "unresolved_facets": len(facet_rows) - len(valid_areas),
                       **facet_area_statistics(valid_areas)},
        "connectivity": {"components": len(components), "significant_components": sum(c["significant"] for c in components),
                         "scc_threshold_fraction": SCC_AREA_FRACTION, "component_records": component_rows},
        "curvature": {"positive_edges": totals["positive_edges"], "negative_edges": totals["negative_edges"],
                      "C_plus_A_rad": totals["positive"], "C_minus_A_rad": totals["negative"],
                      "C_signed_A_rad": totals["signed"], "C_unsigned_A_rad": totals["unsigned"],
                      "L_A": totals["length"], "s_H_signed_rad": totals["signed_mean_rad"],
                      "s_H_signed_deg": totals["signed_mean_deg"], "s_H_unsigned_rad": totals["unsigned_mean_rad"],
                      "s_H_unsigned_deg": totals["unsigned_mean_deg"], "C_H_signed_A_rad": totals["signed"] / 2.0,
                      "positive_beta_length_fraction": positive_length / totals["length"],
                      "negative_beta_length_fraction": negative_length / totals["length"],
                      "positive_abs_C_fraction": totals["positive"] / c_abs,
                      "negative_abs_C_fraction": abs(totals["negative"]) / c_abs},
        "lys15": {"residue_label": "chain I LYS 15" if pdb_path.stem.upper() == "2KAI" else None,
                  "atom_count": len(lys15_ids), "facet_count": len(lys15_facets),
                  "facet_area_A2": sum(facet_area.get(f.generator_ids, 0.0) for f in lys15_facets),
                  "edge_count": len(lys15_edges), "positive_edges": len(lys15_positive),
                  "negative_edges": len(lys15_negative), "C_signed_A_rad": lys15_c,
                  "C_unsigned_A_rad": lys15_cu, "edge_length_A": lys15_l,
                  "length_weighted_signed_beta_rad": lys15_c / lys15_l if lys15_l else None,
                  "length_weighted_signed_beta_deg": math.degrees(lys15_c / lys15_l) if lys15_l else None},
        "beta_populations": _population_stats(distribution_rows),
        "beta_audit": {"max_abs_difference_deg": float(np.max(np.abs(beta_differences))) if len(beta_differences) else None,
                       "rms_difference_deg": float(np.sqrt(np.mean(beta_differences ** 2))) if len(beta_differences) else None,
                       "direct_matches": int(np.sum(direct_matches)),
                       "supplementary_matches": int(np.sum(supplementary_matches)),
                       "sign_disagreements": int(sign_disagreements)},
    }
    (output_dir / "validation_metrics.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (output_dir / "validation_summary.txt").write_text(_case_summary(summary), encoding="utf-8")
    return summary, result


def _population_stats(rows):
    stats = {}
    for label in ("covalent", "noncovalent"):
        values = np.asarray([float(row["abs_beta_deg"]) for row in rows
                             if row["same_group_pair_bond_status"] == label])
        stats[label] = {"N": len(values)}
        if len(values):
            stats[label].update({"mean_deg": float(np.mean(values)), "median_deg": float(np.median(values)),
                                 "std_deg": float(np.std(values)),
                                 "p05_deg": float(np.percentile(values, 5)),
                                 "p25_deg": float(np.percentile(values, 25)),
                                 "p75_deg": float(np.percentile(values, 75)),
                                 "p95_deg": float(np.percentile(values, 95))})
    return stats


def _histogram(rows, bin_width=5.0):
    edges = np.arange(0.0, 180.0 + bin_width, bin_width)
    output = []
    for label in ("covalent", "noncovalent"):
        values = [float(row["abs_beta_deg"]) for row in rows if row["same_group_pair_bond_status"] == label]
        counts, _ = np.histogram(values, bins=edges)
        output.extend({"pair_class": label, "bin_left_deg": float(left),
                       "bin_right_deg": float(right), "count": int(count)}
                      for left, right, count in zip(edges[:-1], edges[1:], counts))
    return output


def _plot_histogram(rows, path, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    bins = np.arange(0.0, 180.0 + 5.0, 5.0)
    covalent = [float(row["abs_beta_deg"]) for row in rows if row["same_group_pair_bond_status"] == "covalent"]
    noncovalent = [float(row["abs_beta_deg"]) for row in rows if row["same_group_pair_bond_status"] == "noncovalent"]
    observed_max = max(covalent + noncovalent, default=90.0)
    display_max = min(180.0, 5.0 * math.ceil(observed_max / 5.0) + 5.0)
    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    ax.hist([covalent, noncovalent], bins=bins, label=[f"Covalent ({len(covalent)})", f"Noncovalent ({len(noncovalent)})"],
            histtype="step", linewidth=1.8)
    ax.set(xlabel=r"$|\beta|$ (degrees)", ylabel="Interior interface edges", title=title,
           xlim=(0, display_max))
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def _case_summary(summary):
    area, curve, lys, conn = (summary["facet_area"], summary["curvature"],
                              summary["lys15"], summary["connectivity"])
    lines = [
        f"Cazals Power-interface validation: {Path(summary['pdb']).stem.upper()}",
        f"Groups A/B: {summary['groups']['A_chains']} / {summary['groups']['B_chains']}",
        f"Input records={summary['inventory']['total_parsed_atom_records']}; protein={summary['inventory']['protein_atom_records']}; "
        f"water={summary['inventory']['water_atom_records']}; other HETATM={summary['inventory']['other_hetatm_records']}; "
        f"retained={summary['inventory']['retained_atoms']}; radius fallbacks={summary['inventory']['fallback_radius_assignments']}",
        f"Settings: probe={summary['settings']['probe_radius_A']:.2f} A, alpha=0, M={summary['settings']['M']:g}",
        f"Interface atoms A/B/total={summary['interface_atoms_A']}/{summary['interface_atoms_B']}/{summary['interface_atoms_total']}; r_AB={summary['r_AB']:.6g}",
        f"Bicolor facets full/alpha0/rejected beta/final={summary['full_regular_bicolor_facets']}/{summary['alpha_zero_bicolor_facets']}/{summary['condition_b_rejected']}/{summary['final_facets']}",
        f"Facet areas: finite={area['bounded_facets']}, unresolved={area['unresolved_facets']}, VIA={area['VIA_A2']:.9g} A^2, mean={area['mean_A2']:.9g}, median={area['median_A2']:.9g}, <1={area['fraction_below_1_A2']:.4%}, >10={area['fraction_above_10_A2']:.4%}, max={area['maximum_A2']:.9g}",
        f"Connectivity: cc={conn['components']}, scc(>=7.5% VIA)={conn['significant_components']}",
        f"Curvature edges +/-={curve['positive_edges']}/{curve['negative_edges']}; L={curve['L_A']:.9g} A",
        f"C+/C-/C_signed/C_unsigned={curve['C_plus_A_rad']:.9g}/{curve['C_minus_A_rad']:.9g}/{curve['C_signed_A_rad']:.9g}/{curve['C_unsigned_A_rad']:.9g} A rad",
        f"s_H signed={curve['s_H_signed_rad']:.9g} rad/{curve['s_H_signed_deg']:.9g} deg; unsigned={curve['s_H_unsigned_rad']:.9g} rad/{curve['s_H_unsigned_deg']:.9g} deg; C_H=C/2={curve['C_H_signed_A_rad']:.9g} A rad",
        f"Positive/negative beta edge-length fractions={curve['positive_beta_length_fraction']:.6%}/{curve['negative_beta_length_fraction']:.6%}",
        f"Positive/negative |C| fractions={curve['positive_abs_C_fraction']:.6%}/{curve['negative_abs_C_fraction']:.6%}",
    ]
    beta_audit = summary["beta_audit"]
    lines.append(f"Beta triangle audit: max/RMS difference={beta_audit['max_abs_difference_deg']:.9g}/{beta_audit['rms_difference_deg']:.9g} deg; direct/supplementary/sign errors={beta_audit['direct_matches']}/{beta_audit['supplementary_matches']}/{beta_audit['sign_disagreements']}")
    if lys["residue_label"]:
        lines.append(f"Lys15: facets={lys['facet_count']}, area={lys['facet_area_A2']:.9g} A^2, edges={lys['edge_count']} (+/- {lys['positive_edges']}/{lys['negative_edges']}), C={lys['C_signed_A_rad']:.9g}, |C|={lys['C_unsigned_A_rad']:.9g} A rad, weighted beta={lys['length_weighted_signed_beta_deg']:.9g} deg")
    for label, values in summary["beta_populations"].items():
        if values["N"]:
            lines.append(f"{label} |beta| N={values['N']}, mean/median/SD={values['mean_deg']:.6g}/{values['median_deg']:.6g}/{values['std_deg']:.6g} deg; p05/p25/p75/p95={values['p05_deg']:.6g}/{values['p25_deg']:.6g}/{values['p75_deg']:.6g}/{values['p95_deg']:.6g}")
    return "\n".join(lines) + "\n"


def main(argv=None):
    """Run the focused 1UDI / 2KAI Cazals validation suite.

    This entry point is intentionally separate from ordinary VorPy network
    solves because it builds weighted regular complexes and detailed audit
    outputs.
    """
    import argparse

    parser = argparse.ArgumentParser(description="Validate VorPy's Cazals/Power interface analysis.")
    parser.add_argument("--system", choices=("1UDI", "2KAI", "all"), default="all")
    parser.add_argument("--output-dir", type=Path, default=Path("output/cazals_validation"))
    args = parser.parse_args(argv)
    pdb_dir = Path(__file__).resolve().parents[2] / "data"
    cases = {
        "1UDI": ({"E"}, {"I"}),
        "2KAI": ({"A", "B"}, {"I"}),
    }
    requested = tuple(cases) if args.system == "all" else (args.system,)
    for name in requested:
        pdb_path = pdb_dir / f"{name}.pdb"
        if not pdb_path.is_file():
            raise SystemExit(f"Required validation PDB is missing: {pdb_path}")
        group_a, group_b = cases[name]
        summary, result = analyze_case(pdb_path, args.output_dir / name, group_a, group_b)
        print((args.output_dir / name / "validation_summary.txt").read_text(encoding="utf-8"))
        if name == "2KAI":
            from vorpy.src.analyze.cazals_edge_forensics import run_2kai_forensics

            prepared = prepare_cazals_pdb(
                pdb_path, group_a, group_b, probe_radius=DEFAULT_PROBE,
                include_hydrogens=False, include_hetatm=False,
                strict_heavy_fallbacks=False, verbose=False,
            )
            forensic_dir = args.output_dir / name / "forensics"
            run_2kai_forensics(result, prepared, forensic_dir)
            print(f"2KAI forensic audit written to {forensic_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
