"""Reproducible comparative-study preflight and atom-model construction.

Phase 1 deliberately gates expensive geometry on resolved method definitions.
No geometry, radii, or curvature mathematics is implemented here. Run with
``python -m vorpy.src.analyze.comparative_geometry_study --audit-only``.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import re
import os
import sys
import time
import shutil
import numpy as np


PROTEIN = frozenset("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL HSD HSE HSP HID HIE HIP ASH GLH LYN CYX MSE SEC PYL".split())
WATER = frozenset("HOH WAT SOL TIP3 TIP TIP4 SPC DOD H2O".split())
IONS = frozenset("NA CL K CA MG ZN FE MN CU CO NI CD CS LI RB SR BA BR I F".split())
IDENTITY_FIELDS = ("chain", "residue_number", "insertion_code", "residue_name", "atom_name")
FROZEN_GROUPS = {"A": ["A", "B"], "B": ["I"]}
# Bump this explicit pipeline revision when a scientific method definition
# changes; editorial/export fixes do not invalidate solved geometry archives.
MODULE_REVISION = "comparative_study_2026_10_05_v1"


@dataclass(frozen=True)
class Atom:
    line: str
    line_number: int
    model: int
    residue_occurrence: int
    record: str
    chain: str
    residue_number: int
    insertion_code: str
    residue_name: str
    atom_name: str
    altloc: str
    element: str
    xyz: tuple[float, float, float]
    category: str

    @property
    def identity(self):
        return (self.chain, self.residue_number, self.insertion_code,
                self.residue_name, self.atom_name)

    @property
    def hydrogen(self):
        return self.element in {"H", "D", "T"}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_csv(path, rows, fields=None):
    rows = list(rows)
    fields = list(fields or dict.fromkeys(k for row in rows for k in row))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            if any(isinstance(v, float) and not math.isfinite(v) for v in row.values()):
                raise ValueError(f"Nonfinite CSV value: {path}")
            writer.writerow({k: json.dumps(v, sort_keys=True, allow_nan=False)
                             if isinstance(v, (dict, list, tuple)) else v
                             for k, v in row.items()})


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")


def audit_record(record):
    """Preserve explicit nonfinite diagnostic sentinels; metrics must be finite."""
    from dataclasses import asdict
    return {key: ("NaN (diagnostic unavailable)" if math.isnan(value) else
                  "+Infinity (unbounded)" if value > 0 else "-Infinity (unbounded)")
            if isinstance(value, float) and not math.isfinite(value) else value
            for key, value in asdict(record).items()}


def read_pdb(path):
    """Audit all records, including solvent recorded as ATOM and wrapped resids.

    Residue occurrences are contiguous runs, not unique residue numbers: the
    supplied GROMACS waters wrap past 9999. Duplicate stable IDs remain visible.
    Multiple models and ambiguous protein alternate locations gate analysis.
    """
    atoms, headers = [], []
    model, occurrence, previous = 1, 0, None
    for line_number, raw in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), 1):
        line = raw.ljust(80)
        record = line[:6].strip()
        if record == "MODEL":
            model = int(line[10:14])
            previous = None
        if record == "TER":
            previous = None
        if record not in {"ATOM", "HETATM"}:
            headers.append(raw)
            continue
        residue_name = line[17:20].strip().upper()
        atom_name = line[12:16].strip().upper()
        category = ("protein" if residue_name in PROTEIN else
                    "water" if residue_name in WATER else
                    "ion" if residue_name in IONS else "other")
        element = line[76:78].strip().upper()
        if not element:
            name = atom_name.lstrip("0123456789")
            element = residue_name if category == "ion" else name[:1]
        residue = (model, line[21].strip(), int(line[22:26]), line[26].strip(), residue_name)
        if residue != previous:
            occurrence += 1
            previous = residue
        xyz = tuple(float(line[start:start + 8]) for start in (30, 38, 46))
        if not all(math.isfinite(v) for v in xyz):
            raise ValueError(f"Nonfinite coordinates: {path}:{line_number}")
        atoms.append(Atom(raw, line_number, model, occurrence, record, residue[1],
                          residue[2], residue[3], residue_name, atom_name,
                          line[16].strip(), element, xyz, category))
    if not atoms:
        raise ValueError(f"No atom records in {path}")
    return atoms, headers


def canonical_identity(atom, chain_map):
    # Explicit nomenclature harmonization only. Coordinates are never changed.
    name = "CD1" if atom.residue_name == "ILE" and atom.atom_name == "CD" else atom.atom_name
    return (chain_map.get(atom.chain, atom.chain), atom.residue_number,
            atom.insertion_code, atom.residue_name, name)


def heavy_atoms(atoms):
    return [a for a in atoms if a.category == "protein" and not a.hydrogen]


def infer_chain_map(atoms, reference):
    """Match chain identity sets without row numbers or geometric alignment.

    A source chain must contain every reference-chain heavy identity after the
    documented ILE alias. Extra atoms are audited separately. Matching must be
    unique and one-to-one. Blank-chain inputs are intentionally not guessed.
    """
    source, target = defaultdict(set), defaultdict(set)
    for atom in heavy_atoms(atoms):
        source[atom.chain].add(canonical_identity(atom, {})[1:])
    for atom in heavy_atoms(reference):
        target[atom.chain].add(atom.identity[1:])
    if "" in source:
        return {}, "unresolved: protein chain IDs are blank"
    mapping = {}
    for chain, identities in source.items():
        matches = [key for key, expected in target.items() if expected <= identities]
        if len(matches) != 1:
            return {}, f"unresolved: chain {chain} has {len(matches)} complete reference matches"
        mapping[chain] = matches[0]
    if len(set(mapping.values())) != len(mapping) or set(mapping.values()) != set(target):
        return {}, "unresolved: chain mapping is not bijective onto the crystal partners"
    return mapping, "unique complete heavy-identity containment, with ILE CD->CD1 alias"


def match_atoms(atoms, reference, chain_map):
    reference_ids = [a.identity for a in heavy_atoms(reference)]
    mapped = [canonical_identity(a, chain_map) for a in heavy_atoms(atoms)]
    original = [(chain_map.get(a.chain, a.chain), *a.identity[1:]) for a in heavy_atoms(atoms)]
    counts = Counter(mapped)
    ref_counts = Counter(reference_ids)
    wanted, observed = set(reference_ids), set(mapped)
    return {
        "reference_heavy_atoms": len(reference_ids), "source_heavy_atoms": len(mapped),
        "matched": len(wanted & observed), "missing": sorted(wanted - observed),
        "extra": sorted(observed - wanted),
        "duplicates": sorted(k for k, n in counts.items() if n > 1),
        "reference_duplicates": sorted(k for k, n in ref_counts.items() if n > 1),
        "matched_before_atom_alias": len(wanted & set(original)),
        "missing_before_atom_alias": sorted(wanted - set(original)),
        "extra_before_atom_alias": sorted(set(original) - wanted),
        "atom_aliases_applied": sum(a != b for a, b in zip(original, mapped)),
    }


def centroid(atoms):
    return [sum(a.xyz[i] for a in atoms) / len(atoms) for i in range(3)] if atoms else None


def biological_groups(headers):
    """Propose the first deposited two-partner identity-transform assembly.

    Uses COMPND molecule order for sign labels and author REMARK 350 assembly
    order to select a single biological complex, not all crystallographic copies.
    Nonidentity transforms are explicitly unsupported by this audit view.
    """
    compounds = " ".join(h[10:].strip() for h in headers if h.startswith("COMPND"))
    molecules = []
    for block in re.split(r"MOL_ID:\s*", compounds)[1:]:
        name = re.search(r"MOLECULE:\s*([^;]+)", block)
        chains = re.search(r"CHAIN:\s*([^;]+)", block)
        if name and chains:
            molecules.append((name[1].strip(), {c.strip() for c in chains[1].split(",")}))
    assemblies, current = {}, None
    for line in headers:
        if not line.startswith("REMARK 350"):
            continue
        if "BIOMOLECULE:" in line:
            current = line.split("BIOMOLECULE:", 1)[1].strip()
            assemblies[current] = {"chains": [], "transforms": []}
        elif current and "APPLY THE FOLLOWING TO CHAINS:" in line:
            assemblies[current]["chains"] += [c.strip() for c in line.split("CHAINS:", 1)[1].split(",")]
        elif current and "BIOMT" in line:
            values = line.split()
            assemblies[current]["transforms"].append([float(v) for v in values[4:8]])
    identity = [[1., 0., 0., 0.], [0., 1., 0., 0.], [0., 0., 1., 0.]]
    for number, assembly in assemblies.items():
        if assembly["transforms"] != identity:
            continue
        partners = [(name, sorted(chains & set(assembly["chains"]))) for name, chains in molecules]
        partners = [(name, chains) for name, chains in partners if chains]
        if len(partners) == 2:
            groups = {"A": partners[0][1], "B": partners[1][1]}
            evidence = (f"PDB COMPND and REMARK 350 biological assembly {number}, identity transform; "
                        f"group A={partners[0][0]}; group B={partners[1][0]}; "
                        f"other deposited assembly chains={ {k: v['chains'] for k, v in assemblies.items() if k != number} }")
            return groups, evidence
    return {"A": [], "B": []}, "unresolved: no two-partner identity-transform biological assembly found"


def audit_input(path, system, kind, atoms, headers, reference):
    protein = [a for a in atoms if a.category == "protein"]
    chain_map, evidence = infer_chain_map(atoms, reference) if system == "2KAI" else ({}, "unresolved: partner definitions required")
    matching = match_atoms(atoms, reference, chain_map) if system == "2KAI" else None
    residues, protein_residues = defaultdict(set), defaultdict(set)
    for a in atoms:
        residues[a.chain or "<blank>"].add(a.residue_occurrence)
        if a.category == "protein":
            protein_residues[a.chain or "<blank>"].add(a.residue_occurrence)
    center = centroid(heavy_atoms(atoms))
    distance = math.sqrt(sum(v*v for v in center)) if center else None
    model_count = len({a.model for a in atoms})
    alternate = sum(bool(a.altloc) for a in protein)
    ready = bool(chain_map and matching is not None and not matching["missing"]
                 and not matching["duplicates"] and not matching["reference_duplicates"]
                 and model_count == 1 and not alternate)
    groups = {g: sorted(k for k, v in chain_map.items() if v in cs) for g, cs in FROZEN_GROUPS.items()}
    if system != "2KAI":
        groups, evidence = biological_groups(headers)
        present_chains = {a.chain for a in protein}
        ready = bool(groups["A"] and groups["B"] and
                     set(groups["A"] + groups["B"]) <= present_chains and
                     model_count == 1 and not alternate and
                     len({a.identity for a in protein}) == len(protein))
    title = " ".join(h[10:].strip() for h in headers if h.startswith("TITLE"))
    title_time = re.search(r"\bt\s*=\s*([0-9.]+)", title)
    filename_time = re.search(r"frame_(\d+)ps", path.stem, re.I)
    time_ps = float(title_time[1]) if title_time else None
    if filename_time and (time_ps is None or time_ps != float(filename_time[1])):
        ready = False
        evidence += "; frame filename/title time mismatch"
    row = {
        "filename": path.name, "path": str(path.resolve()), "sha256": sha256(path),
        "system": system, "kind": kind, "status": "audited", "models": model_count,
        "time_ps": time_ps, "chains": sorted(residues),
        "residues_per_chain": {k: len(v) for k, v in residues.items()},
        "protein_residues_per_chain": {k: len(v) for k, v in protein_residues.items()},
        "ATOM_count": sum(a.record == "ATOM" for a in atoms),
        "HETATM_count": sum(a.record == "HETATM" for a in atoms),
        "protein_heavy_atoms": len(heavy_atoms(atoms)),
        "protein_hydrogens": sum(a.hydrogen for a in protein),
        "waters": len({a.residue_occurrence for a in atoms if a.category == "water"}),
        "water_atoms": sum(a.category == "water" for a in atoms),
        "ions": len({a.residue_occurrence for a in atoms if a.category == "ion"}),
        "ion_atoms": sum(a.category == "ion" for a in atoms),
        "other_atoms": sum(a.category == "other" for a in atoms), "total_atoms": len(atoms),
        "residue_atom_counts": dict(Counter(a.residue_name for a in atoms)),
        "protein_heavy_centroid_A": center, "all_atom_centroid_A": centroid(atoms),
        "protein_centroid_distance_from_origin_A": distance,
        "coordinates_appear_centered": "origin-centered" if distance is not None and distance <= 1.0 else "not origin-centered; box-centered status unresolved",
        "CRYST1": [h for h in headers if h.startswith("CRYST1")],
        "group_A_chains": groups["A"], "group_B_chains": groups["B"],
        "canonical_group_A_chains": FROZEN_GROUPS["A"] if chain_map else [],
        "canonical_group_B_chains": FROZEN_GROUPS["B"] if chain_map else [],
        "chain_map": chain_map, "group_evidence": evidence,
        "all_atoms_needed_for_group_matching_present": ready,
        "protein_altloc_records": alternate,
        "repeated_protein_stable_ids": len(protein) - len({a.identity for a in protein}),
        "Tier_H_matched": matching["matched"] if matching else None,
        "Tier_H_missing": len(matching["missing"]) if matching else None,
        "Tier_H_extra": len(matching["extra"]) if matching else None,
        "Tier_H_atom_aliases": matching["atom_aliases_applied"] if matching else None,
        "selected_partner_heavy_atoms": sum(a.chain in groups["A"] + groups["B"] for a in heavy_atoms(atoms)),
        "group_A_heavy_atoms": sum(a.chain in groups["A"] for a in heavy_atoms(atoms)),
        "group_B_heavy_atoms": sum(a.chain in groups["B"] for a in heavy_atoms(atoms)),
    }
    if system != "2KAI":
        row["canonical_group_A_chains"] = groups["A"]
        row["canonical_group_B_chains"] = groups["B"]
    return row, matching


def write_view(path, atoms, identities=None):
    """Write atom-model view, preserving every selected coordinate text field."""
    lines = ["REMARK 900 ANALYSIS VIEW; SOURCE COORDINATES UNCHANGED; SEE INPUT MANIFEST"]
    for index, atom in enumerate(atoms):
        line = atom.line.ljust(80)
        if identities is not None:
            chain, resid, icode, resname, name = identities[index]
            atom_field = f" {name:<3}" if len(name) < 4 else name[:4]
            line = line[:12] + atom_field + line[16:17] + f"{resname:>3}" + line[20:21] + (chain or " ") + f"{resid:4d}" + (icode or " ") + line[27:]
        lines.append(line)
    lines.append("END")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def construct_views(path, atoms, reference, audit, output):
    """H matches crystal exactly; P contains the complete protein; S all atoms.

    H+solvent is separately labeled HS_context: it retains H protein atoms plus
    the original waters/ions so solvent deltas do not include protein H atoms
    or terminal OXT additions. It is an input view, not a defined geometry.
    """
    if not audit["all_atoms_needed_for_group_matching_present"]:
        return [], []
    mapping = audit["chain_map"]
    by_id = {canonical_identity(a, mapping): a for a in heavy_atoms(atoms)}
    expected = [a.identity for a in heavy_atoms(reference)]
    h_atoms = [by_id[key] for key in expected]
    p_atoms = [a for a in atoms if a.category == "protein"]
    solvent = [a for a in atoms if a.category in {"water", "ion"}]
    variants = {
        "H": (h_atoms, expected, "crystal-matched protein heavy atoms; excludes extra OXT, H, water, ions"),
        "P": (p_atoms, [canonical_identity(a, mapping) for a in p_atoms], "complete simulated protein including protein H and terminal OXT; no solvent"),
        "S": (atoms, [canonical_identity(a, mapping) if a.category == "protein" else a.identity for a in atoms], "complete simulated protein + original solvent and ions; not paired directly with H"),
        "HS_context": (h_atoms + solvent, expected + [a.identity for a in solvent], "exact Tier-H protein coordinates/identities plus original solvent; paired solvent-only contrast"),
    }
    manifest = []
    for tier, (selected, ids, description) in variants.items():
        destination = output / "atom_models" / path.stem / f"{tier}.pdb"
        write_view(destination, selected, ids)
        written, _ = read_pdb(destination)
        if [a.identity for a in written] != ids or [a.xyz for a in written] != [a.xyz for a in selected]:
            raise ValueError(f"Written atom-model identity/coordinate verification failed: {destination}")
        manifest.append({"filename": path.name, "tier": tier, "path": str(destination.resolve()),
                         "atoms": len(selected), "sha256": sha256(destination),
                         "description": description, "coordinates_modified": False,
                         "written_identities_and_coordinates_verified": True,
                         "geometry_status": "not computed; scientific definitions pending"})
    matching_rows = []
    for atom, key in zip(h_atoms, expected):
        matching_rows.append({"filename": path.name, "status": "matched", "canonical_id": key,
                              "source_id": atom.identity, "source_line": atom.line_number,
                              "atom_alias_applied": atom.atom_name != key[-1]})
    return manifest, matching_rows


def discover_inputs(data_dir, frames_dir, systems):
    paths, missing = [], []
    pdbs = sorted(p for p in data_dir.rglob("*") if p.is_file() and p.suffix.lower() == ".pdb")
    for system in systems:
        found = [p for p in pdbs if p.stem.casefold() == system.casefold()]
        if len(found) > 1:
            raise ValueError(f"Multiple candidate static files for {system}: {found}")
        if found:
            paths.append((found[0], system.upper(), "static"))
        else:
            missing.append(system.upper())
    for path in pdbs:
        if path.stem.casefold().startswith("2kai_md_"):
            paths.append((path, "2KAI", "centered_MD_structure"))
    if frames_dir.is_dir():
        for path in sorted(frames_dir.iterdir()):
            if path.is_file() and path.suffix.lower() == ".pdb":
                paths.append((path, "2KAI", "MD_frame"))
    return paths, missing


def run_audit(data_dir, frames_dir, output_dir, systems=("2KAI", "1AK4", "1BRC")):
    data_dir, frames_dir, output = map(Path, (data_dir, frames_dir, output_dir))
    output.mkdir(parents=True, exist_ok=True)
    inputs, missing = discover_inputs(data_dir, frames_dir, systems)
    crystal = next((p for p, s, k in inputs if s == "2KAI" and k == "static"), None)
    if crystal is None:
        candidates = [p for p in data_dir.glob("*.pdb") if p.stem.upper() == "2KAI"]
        crystal = candidates[0] if len(candidates) == 1 else None
    if crystal is None:
        raise ValueError("The frozen 2KAI crystal reference is required for identity auditing")
    reference, _ = read_pdb(crystal)
    audits, manifests, matches, details = [], [], [], {}
    for path, system, kind in inputs:
        atoms, headers = read_pdb(path)
        row, matching = audit_input(path, system, kind, atoms, headers, reference)
        audits.append(row)
        details[str(path.resolve())] = matching
        if matching:
            for category in ("missing", "extra", "duplicates"):
                matches.extend({"filename": path.name, "status": category, "canonical_id": identity}
                               for identity in matching[category])
        if kind == "MD_frame":
            views, matched = construct_views(path, atoms, reference, row, output)
            manifests.extend(views)
            matches.extend(matched)
        elif kind == "static" and row["all_atoms_needed_for_group_matching_present"]:
            selected = [a for a in heavy_atoms(atoms) if a.record == "ATOM" and
                        a.chain in row["group_A_chains"] + row["group_B_chains"]]
            destination = output / "atom_models" / path.stem / "H.pdb"
            write_view(destination, selected)
            written, _ = read_pdb(destination)
            if [a.identity for a in written] != [a.identity for a in selected] or [a.xyz for a in written] != [a.xyz for a in selected]:
                raise ValueError(f"Static analysis view verification failed: {path}")
            manifests.append({"filename": path.name, "tier": "H", "path": str(destination.resolve()),
                              "atoms": len(selected), "sha256": sha256(destination),
                              "description": "selected biological partners, ATOM protein heavy atoms; excludes other crystallographic copies and solvent",
                              "coordinates_modified": False, "written_identities_and_coordinates_verified": True,
                              "geometry_status": "not computed; scientific definitions pending"})
    # Missing requested inputs are explicit rows, never silently omitted.
    audits.extend({"filename": f"{system}.pdb", "system": system, "kind": "static",
                   "status": "missing", "all_atoms_needed_for_group_matching_present": False,
                   "group_evidence": "not determinable: input absent"} for system in missing)
    write_csv(output / "input_audit.csv", audits)
    write_csv(output / "atom_model_manifest.csv", manifests)
    write_csv(output / "atom_identity_matching.csv", matches)
    write_json(output / "atom_matching_details.json", details)
    blockers = []
    blockers.extend(f"Missing structure: {system}.pdb; biological groups cannot be determined from the unavailable file." for system in missing)
    blockers.extend(f"{r['filename']}: {r['group_evidence']}" for r in audits
                    if r.get("status") == "audited" and not r["all_atoms_needed_for_group_matching_present"])
    frames = [r for r in audits if r.get("kind") == "MD_frame"]
    if not frames:
        blockers.append("No MD frame PDB files discovered")
    frame_times = [r["time_ps"] for r in frames]
    if len(set(frame_times)) != len(frame_times):
        blockers.append("Duplicate or missing MD frame times")
    groups = {r["filename"]: {k: r.get(k) for k in ("group_A_chains", "group_B_chains", "chain_map", "group_evidence")}
              for r in audits}
    write_json(output / "group_definitions.json", groups)
    write_json(output / "method_definition_audit.json", {
        "alpha0_selected_power": {"implementation": "power_interface_cli.run_power_interface_pdb",
                                  "alpha": 0.0, "condition_beta_m": 5.0, "probe_radius_A": 1.4},
        "frozen_comparison_power": {"implementation": "vpy_cmnd.Command._run_power_vs_aw -> run_power_interface_pdb",
                                     "distinct_from_alpha0_branch": False,
                                     "evidence": "power_interface_cli defaults alpha=0.0/M=5.0; comparison passes the same power_settings"},
        "full_power": {"status": "user defined", "selection": "all finite/bounded bicolor facets; no alpha or M filtering",
                       "implementation": "power_interface_selection.power_selection_views", "intermediate": "alpha0_power_before_M"},
        "full_power_aw_radii": {"status": "user authorized radius-matched control", "selection": "all bounded bicolor facets",
                                "radii": "exact frozen AW loader radii; weights=r_AW^2; no added probe",
                                "contrast": "AW minus this control isolates weighting geometry under each method's stated eligibility rules"},
        "solvent_context": {
            "power_adapter_limit": "run_power_interface_pdb prepares selected protein chain atoms only; its PDB adapter cannot be used to claim a solvent-present tessellation",
            "low_level_groups": "Power and AW analyzers accept groups that are subsets of all generators; neutral generators must be retained in the full network explicitly",
            "paired_atom_models": "H versus HS_context holds every protein coordinate and atom identity fixed; P versus S is a separate hydrogen-inclusive contrast",
            "boundary_and_radii": "User approved finite explicit-solvent context, no periodic images or artificial boundary; actual AW-loaded radii shared with matched Power",
        },
    })
    source_dir = Path(__file__).resolve().parent
    sources = [source_dir / name for name in ("comparative_geometry_study.py", "power_vs_aw.py", "power_interface_cli.py", "power_interface_curvature.py", "aw_interface_curvature.py")]
    write_json(output / "provenance.json", {
        "python": platform.python_version(), "platform": platform.platform(),
        "data_dir": str(data_dir.resolve()), "frames_dir": str(frames_dir.resolve()),
        "requested_systems": list(systems), "source_sha256": {str(p): sha256(p) for p in sources},
        "identity_fields": IDENTITY_FIELDS,
        "harmonization": ["unique chain heavy-identity mapping", "ILE CD -> CD1", "H excludes extra terminal OXT; P/S retain them"],
        "centering_definition": "protein heavy-atom arithmetic centroid within 1 Angstrom of origin; box centering not inferred from filename or CRYST1",
        "residue_count_definition": "contiguous residue occurrences, including wrapped solvent IDs; protein identity duplicates reported separately",
        "coordinate_policy": "no alignment, recentering, unwrapping, or coordinate edits",
    })
    write_json(output / "phase_status.json", {"phase": 1, "status": "audited; downstream blocked", "blockers": blockers,
                                               "geometry_runs_launched": 0, "frame_count": len(frames), "frame_times_ps": frame_times})
    summary = ["Comparative geometry study — phase 1 input audit", "",
               "OBSERVATIONS", f"Audited {sum(r.get('status') == 'audited' for r in audits)} available structures, including {len(frames)} MD frames.",
               f"Frame title times (ps): {frame_times}",
               "2KAI crystal partners: group A = kallikrein chains A+B; group B = BPTI chain I (PDB COMPND).",
               *[f"{r['filename']}: A={r['group_A_chains']}, B={r['group_B_chains']}. {r['group_evidence']}"
                 for r in audits if r.get("kind") == "static" and r.get("status") == "audited" and r["system"] != "2KAI"],
               "MD frame partner assignments use a unique full heavy-identity chain match; see group_definitions.json.",
               "All atom matching uses (chain, residue number, insertion code, residue name, atom name), not row number.",
               "Raw and harmonized differences are retained in atom_matching_details.json and atom_identity_matching.csv.",
               "H/P/S and paired HS_context views are inputs only; no geometry has been computed.", "",
               "INTERPRETATION LIMITS", "Centered filenames do not prove origin centering, molecule wholeness, or periodic boundary treatment.",
               "Trajectory frames are correlated observations, not independent biological replicates.", "",
               "SCIENTIFIC/INPUT BLOCKERS (stop before dependent geometry)", *[f"- {b}" for b in blockers], "",
               "Phases 2–6 and numerical/figure exports have not been run. No scientific metrics or placeholder figures were fabricated."]
    (output / "comparative_study_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    return audits, blockers


def loaded_atom_id(row):
    seq = row["auth_seq_id"] if str(row["auth_seq_id"]).strip() else row["res_seq"]
    return (str(row["chain_name"]).strip(), int(seq),
            str(row["pdb_ins_code"]).strip(), str(row["res_name"]).strip(), str(row["name"]).strip())


def radii_audit(system, prepared, destination):
    """Join actual loader/analyzer radii by stable atom identity, never row."""
    loaded = list(system.balls.iterrows())
    used = set()
    rows = []
    for atom in prepared.atoms:
        identity = (atom.chain, atom.residue_number, atom.insertion_code, atom.residue_name, atom.atom_name)
        # Some legacy PDB reader paths lose insertion codes. Resolve those
        # cases by the remaining stable identity fields plus exact coordinates;
        # never use dataframe row order as an atom match.
        coarse = [int(index) for index, row in loaded
                  if str(row["chain_name"]).strip() == atom.chain
                  and int(row["auth_seq_id"] if str(row["auth_seq_id"]).strip() else row["res_seq"]) == atom.residue_number
                  and str(row["res_name"]).strip() == atom.residue_name
                  and str(row["name"]).strip() == atom.atom_name]
        matches = [index for index in coarse if tuple(map(float, system.balls.at[index, "loc"])) == (atom.x, atom.y, atom.z)]
        if len(matches) != 1 or matches[0] in used:
            raise ValueError(f"AW/Power stable-identity/coordinate join is not unique for {identity}: {len(matches)} candidates")
        index = matches[0]
        used.add(index)
        aw = system.balls.loc[index]
        rows.append({"atom_id": identity, "power_generator_id": atom.index, "aw_generator_id": index,
                     "power_base_radius_A": atom.base_radius, "power_expanded_radius_A": atom.expanded_radius,
                     "aw_radius_A": float(aw["rad"]), "aw_minus_power_base_A": float(aw["rad"]) - atom.base_radius,
                     "join_note": "exact full identity" if str(aw["pdb_ins_code"]).strip() == atom.insertion_code else "loader insertion code absent; uniquely disambiguated by exact coordinates"})
    if len(rows) != len(loaded) or len(used) != len(loaded):
        raise ValueError("AW/Power atom sets differ")
    write_csv(destination, rows)
    return rows


def metric_row(method, area, atoms, groups, interface_atoms, elements, candidates, eligible,
               components, edge, smooth, combined, exclusions, unresolved, radius_model):
    unsigned, signed = edge["unsigned"], edge["signed"]
    return {
        "method": method, "radius_model": radius_model, "atom_count": atoms,
        "group_A_atoms": groups[0], "group_B_atoms": groups[1],
        "interface_atom_count": interface_atoms, "interface_element_count": elements,
        "candidate_elements": candidates, "eligible_elements": eligible, "selected_elements": elements,
        "excluded_elements": candidates - elements, "excluded_reasons": exclusions,
        "connected_components": components, "interface_area_A2": area,
        "area_per_interface_atom_A2": area / interface_atoms if interface_atoms else None,
        "signed_raw_edge_A_rad": signed, "unsigned_raw_edge_A_rad": unsigned,
        "positive_edge_A_rad": edge["positive"], "negative_edge_A_rad": edge["negative"],
        "smooth_H_A": smooth, "conventional_edge_H_A": 0.5 * signed, "combined_H_A": combined,
        "H_per_area_inv_A": combined / area if area else None,
        "unsigned_edge_per_area_rad_per_A": unsigned / area if area else None,
        "cancellation_fraction": 1.0 - abs(signed) / unsigned if unsigned else None,
        "unresolved_quadratures": unresolved,
        "normalization_status": "defined" if area and interface_atoms and unsigned else "zero denominator; blank is undefined, not dropped",
        "edge_domain": "interior seams shared by two selected interface patches; open-boundary curvature not assigned",
    }


def extract_power_metrics(method, result, polygons):
    from vorpy.src.analyze.cazals_validation import interface_components, interface_atom_ids
    selected = result.final_facets
    if any(polygons[f.generator_ids]["area_A2"] is None for f in selected):
        raise ValueError(f"Unresolved selected Power facet area in {method}; cannot silently drop")
    areas = {f.generator_ids: polygons[f.generator_ids]["area_A2"] for f in selected}
    a, b = interface_atom_ids(selected, result.group_a)
    edge = result.edge_totals
    return metric_row(method, sum(areas.values()), len(result.points),
                      (len(result.group_a), len(result.group_b)), len(a | b), len(selected), len(result.facets),
                      sum(p["status"] == "finite_bounded" for p in polygons.values()),
                      len(interface_components(result, areas)), edge, result.smooth_surface_curvature,
                      0.5 * edge["signed"], dict(Counter(f.reason for f in result.facets if not f.final_selected)),
                      0, "frozen_Cazals_Chothia_R=r+1.4")


def extract_aw_metrics(result):
    from vorpy.src.analyze.power_vs_aw import _aw_surface_components
    return metric_row("aw", result.interface_area, len(result.network.balls),
                      (len(result.group_a), len(result.group_b)),
                      len({i for s in result.selected_surfaces for i in s.generator_ids}),
                      len(result.selected_surfaces), len(result.surfaces),
                      sum(s.bounded and s.feature_complete and s.area is not None for s in result.surfaces),
                      len(_aw_surface_components(result)), result.edge_curvature, result.surface_curvature,
                      result.combined_curvature, dict(Counter(s.reason for s in result.surfaces if not s.included)),
                      sum(e.integration_status == "maximum_order_reached" for e in result.edges),
                      "frozen_VorPy_AW_loader")


def build_or_load_aw(view, output):
    from vorpy.src.system import System
    from vorpy.src.group import Group
    from vorpy.src.io import load_network, save_network
    source = Path(__file__).resolve().parents[1]
    geometry_files = sorted(p for folder in ("network", "calculations", "inputs", "chemistry", "objects", "group", "system")
                            for p in (source / folder).rglob("*.py"))
    settings = {"net_type": "aw", "max_vert": 5.0, "surf_res": 0.2,
                "box_size": 1.25, "num_splits": 24, "build_type": "all"}
    key = {"input_sha256": sha256(view), "settings": dict(settings),
           "geometry_source_sha256": {str(p.relative_to(source)): sha256(p) for p in geometry_files}}
    archive, metadata = output / "full_aw.vpy", output / "full_aw_cache.json"
    if archive.exists() and metadata.exists():
        saved = json.loads(metadata.read_text())
        saved["settings"] = {k: saved["settings"].get(k) for k in settings}
        if saved == key:
            print(f"Reusing solved AW archive: {archive}", flush=True)
            return load_network(archive)
    system = System(file=str(view), make_dir=False)
    system.files["dir"] = str(output)
    group = Group(sys=system, name="comparative_full_system", atoms=list(system.balls.index),
                  settings=dict(settings), make_net=True, build_net=False)
    cwd = Path.cwd()
    try:
        group.build()
        save_network(group.net, archive)
    finally:
        os.chdir(cwd)
    write_json(metadata, key)
    return group.net


def validate_reversal(power, aw, views, polygons):
    from vorpy.src.analyze.power_interface_curvature import analyze_power_regular_complex
    from vorpy.src.analyze.power_interface_selection import power_selection_views
    from vorpy.src.analyze.aw_interface_curvature import analyze_aw_interface_curvature
    reverse = analyze_power_regular_complex(power.points, power.expanded_radii, power.group_b, power.group_a,
                                            power.full_simplices, power.filtration)
    reverse_views, reverse_polygons = power_selection_views(reverse)
    checks = []
    for method, forward in views.items():
        f, r = extract_power_metrics(method, forward, polygons), extract_power_metrics(method, reverse_views[method], reverse_polygons)
        for key in ("interface_area_A2", "unsigned_raw_edge_A_rad", "signed_raw_edge_A_rad", "conventional_edge_H_A", "combined_H_A"):
            sign = -1 if key in {"signed_raw_edge_A_rad", "conventional_edge_H_A", "combined_H_A"} else 1
            passed = math.isclose(f[key], sign * r[key], rel_tol=1e-10, abs_tol=1e-8)
            checks.append({"method": method, "check": f"AB_reversal_{key}", "passed": passed})
    reversed_aw = analyze_aw_interface_curvature(aw.network, aw.group_b, aw.group_a, selection_mode=aw.selection_mode)
    f, r = extract_aw_metrics(aw), extract_aw_metrics(reversed_aw)
    for key in ("interface_area_A2", "unsigned_raw_edge_A_rad", "smooth_H_A", "signed_raw_edge_A_rad", "conventional_edge_H_A", "combined_H_A"):
        sign = 1 if key in {"interface_area_A2", "unsigned_raw_edge_A_rad"} else -1
        checks.append({"method": "aw", "check": f"AB_reversal_{key}",
                       "passed": math.isclose(f[key], sign * r[key], rel_tol=1e-10, abs_tol=1e-8)})
    return checks


def coordinate_index(atom_ids, balls):
    """Build Power-row IDs from full-network IDs after exact coordinate checks."""
    ids = sorted(map(int, balls.index))
    coordinates = np.asarray(balls.loc[ids, "loc"].tolist(), dtype=float)
    keys = {}
    for index, point in enumerate(coordinates):
        keys.setdefault(tuple(map(float, point)), []).append(index)
    return ids, coordinates, keys


def solvated_geometry_case(frame_audit, output, crystal_atoms):
    """Finite explicit-solvent H vs H+S pair, with protein A/B and water separate."""
    from collections import defaultdict
    from vorpy.src.analyze.aw_interface_curvature import analyze_aw_interface_curvature
    from vorpy.src.analyze.power_interface_curvature import analyze_power_interface_curvature
    from vorpy.src.analyze.power_interface_selection import power_selection_views
    from vorpy.src.analyze.power_vs_aw import _aw_surface_components

    path = Path(frame_audit["path"])
    atoms, _ = read_pdb(path)
    h_atoms = heavy_atoms(atoms)
    mapping, reason = infer_chain_map(atoms, crystal_atoms)
    if not mapping:
        raise ValueError(f"Solvent frame chain mapping unresolved: {path.name}: {reason}")
    canonical = [canonical_identity(a, mapping) for a in h_atoms]
    index_by_id = {identity: index for index, identity in enumerate(canonical)}
    protein_a = {i for i, key in enumerate(canonical) if key[0] in FROZEN_GROUPS["A"]}
    protein_b = {i for i, key in enumerate(canonical) if key[0] in FROZEN_GROUPS["B"]}
    protein = protein_a | protein_b
    solvent = [a for a in atoms if a.category in {"water", "ion"}]
    waters = [a for a in solvent if a.category == "water" and a.element == "O"]
    if not waters:
        raise ValueError(f"No water oxygens available for solvent boundary audit: {path.name}")
    hs_path = Path(output) / "atom_models" / path.stem / "HS_context.pdb"
    views = [a for a in h_atoms] + solvent
    # Explicitly verify the paired protein rows, stable IDs, coordinates and radii.
    reread, _ = read_pdb(hs_path)
    if len(reread) != len(views) or [a.xyz for a in reread[:len(h_atoms)]] != [a.xyz for a in h_atoms]:
        raise ValueError(f"H+S view does not preserve its paired H protein: {path.name}")
    frame_out = Path(output) / "solvent" / path.stem
    case_rows, excluded_rows, boundary_rows = [], [], []
    for tier, view_path in (("H", Path(output) / "atom_models" / path.stem / "H.pdb"),
                            ("H+S", hs_path)):
        case = Path(output) / "cases" / path.stem / tier
        case.mkdir(parents=True, exist_ok=True)
        if tier == "H":
            # Reuse the already validated Tier-H analysis from Phase 4. The
            # paired solvent study must change context only, not repeat the
            # protein-only geometry calculation.
            baseline = read_csv(case / "metrics.csv")
            case_rows.extend({**row, "tier": "H", "context": "desolvated",
                              "interface_kind": "protein_protein"}
                             for row in baseline
                             if row.get("method") in {"aw", "full_power_aw_radii"})
            continue
        net = build_or_load_aw(view_path, case)
        all_ids, points, coordinate_lookup = coordinate_index(None, net.balls)
        protein_ids = {}
        # Match every protein generator by stable molecular identity and exact coordinates.
        source_by_key = {canonical_identity(a, mapping): a for a in h_atoms}
        for network_id in all_ids:
            row = net.sys.balls.loc[network_id]
            if str(row["res_name"]).upper() not in PROTEIN:
                continue
            chain = mapping.get(str(row["chain_name"]).strip(), str(row["chain_name"]).strip())
            seq = int(row["auth_seq_id"] if str(row["auth_seq_id"]).strip() else row["res_seq"])
            resname, atomname = str(row["res_name"]).strip().upper(), str(row["name"]).strip().upper()
            atomname = "CD1" if resname == "ILE" and atomname == "CD" else atomname
            xyz = tuple(map(float, net.balls.at[network_id, "loc"]))
            candidates = [key for key, atom in source_by_key.items()
                          if key[0] == chain and key[1] == seq and key[3] == resname and key[4] == atomname
                          and atom.xyz == xyz]
            if len(candidates) != 1:
                raise ValueError(f"Full-system AW stable-identity/coordinate join is not unique for {(chain, seq, resname, atomname)}")
            key = candidates[0]
            protein_ids[network_id] = index_by_id[key]
        if len(protein_ids) != len(h_atoms):
            raise ValueError("H and H+S protein generators do not match one-to-one")
        aw_a = {n for n, h in protein_ids.items() if h in protein_a}
        aw_b = {n for n, h in protein_ids.items() if h in protein_b}
        water_ids = {i for i in all_ids if str(net.sys.balls.at[i, "res_name"]).upper() in WATER}
        ions = {i for i in all_ids if str(net.sys.balls.at[i, "res_name"]).upper() in IONS}
        if tier == "H+S" and (not water_ids or not ions):
            raise ValueError("Expected both water and ion surroundings in H+S AW network")
        aw_pp = analyze_aw_interface_curvature(net, aw_a, aw_b, selection_mode="supported")
        aw_context = {"system": "2KAI", "filename": path.name, "time_ps": frame_audit["time_ps"],
                      "tier": tier, "context": "finite_explicit_solvent_context" if tier == "H+S" else "desolvated", "method": "aw"}
        case_rows.append({**aw_context, **extract_aw_metrics(aw_pp)})
        aw.export_csv(case / "aw_protein_protein")
        if tier == "H+S":
            aw_pw = analyze_aw_interface_curvature(net, aw_a | aw_b, water_ids, selection_mode="supported")
            pwrow = {**aw_context, "interface_kind": "protein_water", "method": "aw_protein_water",
                     **extract_aw_metrics(aw_pw)}
            case_rows.append(pwrow)
            aw_pw.export_csv(case / "aw_protein_water")
            for record in aw_pw.surfaces:
                if not record.included:
                    excluded_rows.append({**aw_context, "method": "aw_protein_water", "feature": "surface",
                                          "feature_id": record.surface_id, "reason": record.reason})
        from vorpy.src.analyze.power_interface_cli import prepare_cazals_pdb
        # The H/H+S power branch uses exactly the loaded AW radii, no probe addition.
        ordered_points = np.asarray(net.balls.loc[all_ids, "loc"].tolist(), dtype=float)
        radii = np.asarray([float(net.balls.at[i, "rad"]) for i in all_ids], dtype=float)
        power_index = {network_id: row_index for row_index, network_id in enumerate(all_ids)}
        p_a, p_b = {power_index[i] for i in aw_a}, {power_index[i] for i in aw_b}
        power_pp = analyze_power_interface_curvature(ordered_points, radii, p_a, p_b)
        ppviews, polygons = power_selection_views(power_pp)
        pfull = ppviews["full_power"]
        row = extract_power_metrics("full_power_aw_radii", pfull, polygons)
        row["radius_model"] = "exact_frozen_AW_radii; w=r_AW^2; no_probe"
        case_rows.append({**aw_context, "method": "full_power_aw_radii", **row})
        pfull.export(case / "power_protein_protein", export_filtration=False)
        if tier == "H+S":
            pwater = {power_index[i] for i in water_ids}
            pwpower = analyze_power_interface_curvature(ordered_points, radii, p_a | p_b, pwater)
            pwviews, pwp = power_selection_views(pwpower)
            pwfull = pwviews["full_power"]
            pwrow = extract_power_metrics("full_power_aw_radii_protein_water", pwfull, pwp)
            pwrow["radius_model"] = "exact_frozen_AW_radii; water group B; ions neutral"
            case_rows.append({**aw_context, "interface_kind": "protein_water",
                              "method": "full_power_aw_radii_protein_water", **pwrow})
            pwfull.export(case / "power_protein_water", export_filtration=False)
            for feature, polygon in pwp.items():
                if polygon["status"] != "finite_bounded":
                    excluded_rows.append({**aw_context, "method": "full_power_aw_radii_protein_water",
                                          "feature": "facet", "feature_id": feature, "reason": polygon["status"]})
        for record in aw_pp.surfaces:
            if not record.included:
                excluded_rows.append({**aw_context, "method": "aw", "feature": "surface",
                                      "feature_id": record.surface_id, "reason": record.reason})
        solvent_centers = np.asarray([a.xyz for a in waters], dtype=float)
        from vorpy.src.analyze.finite_solvent_boundary import solvent_hull, hull_clearances
        _, equations = solvent_hull(solvent_centers)
        if tier == "H+S":
            samples, reaches = {}, {}
            max_solvent_radius = max(float(net.balls.at[i, "rad"]) for i in water_ids | ions)
            for record in aw_pp.selected_surfaces:
                surface = net.surfs.loc[record.surface_id]
                raw_points = surface.get("points")
                sample = np.asarray(raw_points, dtype=float)
                if sample.ndim == 3 and sample.shape[-1] == 3:
                    sample = sample.reshape((-1, 3))
                if sample.ndim != 2 or sample.shape[1] != 3 or not len(sample):
                    sample = np.asarray([net.balls.at[i, "loc"] for i in record.generator_ids], dtype=float)
                vertex_ids = surface.get("verts")
                if vertex_ids is not None:
                    extra = [net.verts.at[int(i), "loc"] for i in vertex_ids if int(i) in net.verts.index]
                    if extra:
                        sample = np.vstack((sample, np.asarray(extra, dtype=float)))
                a_id, b_id = record.generator_ids
                arow, brow = net.balls.loc[a_id], net.balls.loc[b_id]
                distance = np.linalg.norm(sample - np.asarray(arow["loc"], dtype=float), axis=1)
                q = np.max(distance - float(arow["rad"]))
                samples[f"aw_surface_{record.surface_id}"] = sample
                reaches[f"aw_surface_{record.surface_id}"] = max(0., float(q)) + max_solvent_radius
            from vorpy.src.analyze.finite_solvent_boundary import audit_feature_samples
            boundary_rows.extend({**aw_context, "method": "aw", **r} for r in audit_feature_samples(samples, equations, reaches))
            # Power facet polygons provide exact sampled dual-polygon vertices.
            p_samples, p_reaches = {}, {}
            max_water_power_radius = max(radii[power_index[i]] for i in water_ids | ions)
            triangle_to_tets = defaultdict(set)
            for tet in power_pp.full_simplices[3]:
                for tri in combinations(tet, 3):
                    triangle_to_tets[tri].add(tet)
            for pair, polygon in polygons.items():
                if polygon["status"] != "finite_bounded" or not polygon.get("polygon_vertices"):
                    continue
                sample = np.asarray(polygon["polygon_vertices"], dtype=float)
                local_ids = next((f.incident_regular_triangles for f in power_pp.facets if f.generator_ids == pair), ())
                point_ids = pair
                incident_tets = {tet for tri in local_ids for tet in triangle_to_tets[tri]}
                power_levels = [power_pp.power_vertices[t].power_value for t in incident_tets if t in power_pp.power_vertices]
                reach = max((np.sqrt(max(0., value + max_water_power_radius**2)) for value in power_levels), default=float("inf"))
                aindex, bindex = pair
                stablepair = (all_ids[aindex], all_ids[bindex])
                p_samples[f"power_facet_{stablepair[0]}_{stablepair[1]}"] = sample
                p_reaches[f"power_facet_{stablepair[0]}_{stablepair[1]}"] = float(reach)
            boundary_rows.extend({**aw_context, "method": "full_power_aw_radii", **r}
                                 for r in audit_feature_samples(p_samples, equations, p_reaches))
            for method in ("aw", "full_power_aw_radii"):
                own = [r for r in boundary_rows if r["time_ps"] == frame_audit["time_ps"] and r["method"] == method]
                boundary_rows.append({**aw_context, "method": method, "feature_id": "FRAME_BOUNDARY_SUMMARY",
                                      "sample_count": sum(r["sample_count"] for r in own),
                                      "minimum_hull_distance_A": min((r["minimum_hull_distance_A"] for r in own if r["minimum_hull_distance_A"] is not None), default=None),
                                      "minimum_boundary_margin_A": min((r["minimum_boundary_margin_A"] for r in own if r["minimum_boundary_margin_A"] is not None), default=None),
                                      "boundary_status": "boundary-sensitive" if any(r["boundary_status"] != "interior_at_sampled_geometry" for r in own) else "interior_at_sampled_geometry"})
    return case_rows, excluded_rows, boundary_rows


def run_geometry_case(audit, output_dir, *, validate_frozen=False):
    """Run unchanged frozen baselines plus bounded/alpha-only Power selections."""
    from dataclasses import asdict
    from vorpy.src.analyze.power_interface_cli import run_power_interface_pdb
    from vorpy.src.analyze.power_interface_selection import power_selection_views
    from vorpy.src.analyze.aw_interface_curvature import analyze_aw_interface_curvature
    output_dir = Path(output_dir).resolve()
    name = Path(audit["filename"]).stem
    view = output_dir / "atom_models" / name / "H.pdb"
    case = output_dir / "cases" / name / "H"
    case.mkdir(parents=True, exist_ok=True)
    completion = case / "complete.json"
    if completion.exists():
        saved = json.loads(completion.read_text())
        if saved.get("module_revision") == MODULE_REVISION and saved.get("input_sha256") == audit["sha256"]:
            print(f"CASE {name}: reusing validated metrics", flush=True)
            return read_csv(case / "metrics.csv")
    source = Path(audit["path"]) if audit["kind"] == "static" else view
    print(f"CASE {name}: frozen Cazals geometry", flush=True)
    run = run_power_interface_pdb(source, set(audit["canonical_group_A_chains"]),
                                  set(audit["canonical_group_B_chains"]), case / "cazals_alpha0_power")
    power, prepared = run["result"], run["prepared"]
    views, polygons = power_selection_views(power)
    power.export(case / "regular_geometry")
    rows = [extract_power_metrics(method, result, polygons) for method, result in views.items()]
    excluded, memberships = [], []
    stable = {a.index: (a.chain, a.residue_number, a.insertion_code, a.residue_name, a.atom_name) for a in prepared.atoms}
    for method, result in views.items():
        method_dir = case / method
        method_dir.mkdir(exist_ok=True)
        write_csv(method_dir / "facet_selection.csv", [{"stable_atom_ids": [stable[i] for i in f.generator_ids],
                   "selected": f.final_selected, "area_A2": polygons[f.generator_ids]["area_A2"],
                   "geometry_status": polygons[f.generator_ids]["status"], **audit_record(f)} for f in result.facets])
        write_csv(method_dir / "measured_edges.csv", [audit_record(e) for e in result.edges])
        for f in result.facets:
            if not f.final_selected:
                excluded.append({"method": method, "feature": "facet", "id": f.generator_ids,
                                 "stable_atom_ids": [stable[i] for i in f.generator_ids], "reason": f.reason,
                                 "geometry_status": polygons[f.generator_ids]["status"]})
        for e in result.edges:
            if e.status != "included":
                excluded.append({"method": method, "feature": "edge", "id": e.regular_triangle,
                                 "reason": e.reason, "geometry_status": e.geometry_status})
    print(f"CASE {name}: full-system AW geometry", flush=True)
    network = build_or_load_aw(view, case)
    radii = radii_audit(network.sys, prepared, case / "radii_audit.csv")
    mapping = {r["power_generator_id"]: r["aw_generator_id"] for r in radii}
    aw = analyze_aw_interface_curvature(network, {mapping[i] for i in power.group_a},
                                        {mapping[i] for i in power.group_b}, selection_mode="supported")
    aw.export_csv(case / "aw")
    rows.append(extract_aw_metrics(aw))
    from vorpy.src.analyze.power_interface_curvature import analyze_power_interface_curvature
    control_radii = [float(network.balls.at[mapping[i], "rad"]) for i in range(len(prepared.atoms))]
    control_base = analyze_power_interface_curvature(prepared.points, control_radii, power.group_a, power.group_b)
    control_views, control_polygons = power_selection_views(control_base)
    control = control_views["full_power"]
    control_row = extract_power_metrics("full_power_aw_radii", control, control_polygons)
    control_row["radius_model"] = "exact_frozen_AW_radii; w=r_AW^2; no_probe"
    rows.append(control_row)
    control_dir = case / "full_power_aw_radii"
    control_dir.mkdir(exist_ok=True)
    write_csv(control_dir / "facet_selection.csv", [{"stable_atom_ids": [stable[i] for i in f.generator_ids],
              "area_A2": control_polygons[f.generator_ids]["area_A2"],
              "geometry_status": control_polygons[f.generator_ids]["status"], **audit_record(f)} for f in control.facets])
    write_csv(control_dir / "measured_edges.csv", [audit_record(e) for e in control.edges])
    for f in control.facets:
        if not f.final_selected:
            excluded.append({"method": "full_power_aw_radii", "feature": "facet", "id": f.generator_ids,
                             "stable_atom_ids": [stable[i] for i in f.generator_ids], "reason": f.reason,
                             "geometry_status": control_polygons[f.generator_ids]["status"]})
    for e in control.edges:
        if e.status != "included":
            excluded.append({"method": "full_power_aw_radii", "feature": "edge", "id": e.regular_triangle,
                             "reason": e.reason, "geometry_status": e.geometry_status})
    selection_sets = {method: {f.generator_ids for f in result.final_facets} for method, result in views.items()}
    write_csv(case / "selection_membership.csv", [{"stable_atom_ids": [stable[i] for i in f.generator_ids],
              "generator_ids": f.generator_ids, "area_A2": polygons[f.generator_ids]["area_A2"],
              **{method: f.generator_ids in selected for method, selected in selection_sets.items()}}
              for f in power.facets])
    for s in aw.surfaces:
        if not s.included:
            excluded.append({"method": "aw", "feature": "surface", "id": s.surface_id,
                             "reason": s.reason, **audit_record(s)})
    for e in aw.edges:
        if e.status != "included":
            excluded.append({"method": "aw", "feature": "edge", "id": e.edge_id,
                             "reason": e.reason, **audit_record(e)})
    context = {"system": audit["system"], "filename": audit["filename"], "kind": audit["kind"],
               "tier": "H", "time_ps": audit["time_ps"], "context": "protein_only",
               "input_sha256": audit["sha256"]}
    rows = [{**context, **r} for r in rows]
    write_csv(case / "metrics.csv", rows)
    write_csv(case / "excluded_geometry.csv", [{**context, **r} for r in excluded])
    checks = validate_reversal(power, aw, views, polygons)
    # The control uses the same validated sign/measurement code, with its own
    # solved regular geometry; verify it independently as well.
    from vorpy.src.analyze.power_interface_curvature import analyze_power_regular_complex
    reverse_control_base = analyze_power_regular_complex(control_base.points, control_base.expanded_radii,
        control_base.group_b, control_base.group_a, control_base.full_simplices, control_base.filtration)
    reverse_control = power_selection_views(reverse_control_base)[0]["full_power"]
    checks.extend([
        {"method": "full_power_aw_radii", "check": "AB_reversal_signed", "passed": math.isclose(control.edge_totals["signed"], -reverse_control.edge_totals["signed"], abs_tol=1e-8)},
        {"method": "full_power_aw_radii", "check": "AB_reversal_unsigned", "passed": math.isclose(control.edge_totals["unsigned"], reverse_control.edge_totals["unsigned"], abs_tol=1e-8)},
        {"method": "full_power_aw_radii", "check": "AB_reversal_membership", "passed": {f.generator_ids for f in control.final_facets} == {f.generator_ids for f in reverse_control.final_facets}},
    ])
    if validate_frozen:
        expected = {"cazals_alpha0_power": {"interface_area_A2": 873.762685, "signed_raw_edge_A_rad": -127.050015,
                    "unsigned_raw_edge_A_rad": 607.364973, "combined_H_A": -63.5250077},
                    "aw": {"interface_area_A2": 949.385953, "smooth_H_A": 1.88673524,
                    "signed_raw_edge_A_rad": -105.41778802, "unsigned_raw_edge_A_rad": 686.17775580,
                    "conventional_edge_H_A": -52.70889401, "combined_H_A": -50.82215877,
                    "selected_elements": 404, "candidate_elements": 480}}
        for row in rows:
            for key, value in expected.get(row["method"], {}).items():
                checks.append({"method": row["method"], "check": f"frozen_{key}", "expected": value, "actual": row[key],
                               "passed": math.isclose(row[key], value, rel_tol=0, abs_tol=1e-6)})
        previous = Path(__file__).resolve().parents[3] / "output/2KAI-orientation-report/2KAI/power_interface"
        for path in (case / "cazals_alpha0_power").glob("power_interface_*.csv"):
            checks.append({"method": "cazals_alpha0_power", "check": f"byte_identity_{path.name}",
                           "passed": (previous / path.name).exists() and sha256(path) == sha256(previous / path.name)})
    write_csv(case / "validation.csv", checks)
    if not all(c["passed"] for c in checks):
        raise ValueError(f"Phase validation failed for {name}; see {case / 'validation.csv'}")
    if audit["kind"] == "static" or audit["filename"].casefold() == "frame_0100ps.pdb":
        create_visualization_exports(audit, output_dir, rows, power_run=run,
                                     network=network, aw_result=aw,
                                     control_power=control_base,
                                     control_polygons=control_polygons)
    print(f"CASE {name}: validation passed", flush=True)
    write_json(completion, {"module_revision": MODULE_REVISION, "input_sha256": audit["sha256"], "validated": True})
    return rows


def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        for key, value in tuple(row.items()):
            if value and value[:1] in "[{":
                try:
                    row[key] = json.loads(value)
                except json.JSONDecodeError:
                    pass
    return rows


STUDY_METRICS = ("interface_area_A2", "interface_atom_count", "selected_elements", "connected_components",
                 "signed_raw_edge_A_rad", "unsigned_raw_edge_A_rad", "smooth_H_A", "conventional_edge_H_A",
                 "combined_H_A", "H_per_area_inv_A", "unsigned_edge_per_area_rad_per_A",
                 "positive_edge_A_rad", "negative_edge_A_rad", "cancellation_fraction", "unresolved_quadratures")


def aggregate_study(output):
    """Aggregate saved CSVs only; missing/undefined values are explicit counts."""
    import statistics
    output = Path(output)
    rows, exclusions = [], []
    for path in sorted((output / "cases").glob("*/H/metrics.csv")):
        rows.extend(read_csv(path))
        exclusions.extend(read_csv(path.with_name("excluded_geometry.csv")))
    static = [r for r in rows if r["kind"] == "static"]
    trajectory = [r for r in rows if r["kind"] == "MD_frame"]
    write_csv(output / "static_method_comparison.csv", static)
    write_csv(output / "trajectory_metrics.csv", trajectory)
    write_csv(output / "excluded_geometry_audit.csv", exclusions)
    summaries = []
    for method in sorted({r["method"] for r in trajectory}):
        selected = [r for r in trajectory if r["method"] == method]
        crystal = next(r for r in static if r["system"] == "2KAI" and r["method"] == method)
        for metric in STUDY_METRICS:
            defined = [float(r[metric]) for r in selected if r[metric] != ""]
            if not all(math.isfinite(v) for v in defined):
                raise ValueError(f"Nonfinite metric {method}/{metric}")
            if len(defined) != len(selected):
                summaries.append({"method": method, "metric": metric, "n_frames": len(selected),
                                  "n_undefined": len(selected) - len(defined), "status": "incomplete; no partial aggregate"})
                continue
            mean = statistics.mean(defined)
            sd = statistics.stdev(defined) if len(defined) > 1 else None
            baseline = float(crystal[metric]) if crystal[metric] != "" else None
            nonnegative = all(v >= 0 for v in defined)
            summaries.append({"method": method, "metric": metric, "n_frames": len(selected), "n_undefined": 0,
                              "mean": mean, "sd": sd, "median": statistics.median(defined),
                              "min": min(defined), "max": max(defined),
                              "coefficient_of_variation": sd / mean if sd is not None and mean > 0 and nonnegative else None,
                              "CV_status": "defined for nonnegative metric with positive mean" if sd is not None and mean > 0 and nonnegative else "not meaningful for signed/zero-mean or single-frame metric",
                              "crystal": baseline, "MD_mean_minus_crystal": mean - baseline if baseline is not None else None,
                              "status": "descriptive correlated trajectory; no biological significance test"})
    write_csv(output / "trajectory_summary.csv", summaries)
    effects, selections = [], []
    for key in sorted({(r["system"], r["filename"]) for r in rows}):
        case_rows = {r["method"]: r for r in rows if (r["system"], r["filename"]) == key}
        if "aw" not in case_rows or "full_power_aw_radii" not in case_rows:
            raise ValueError(f"Missing radius-matched control for {key}")
        aw, control = case_rows["aw"], case_rows["full_power_aw_radii"]
        for metric in STUDY_METRICS:
            effect = float(aw[metric]) - float(control[metric]) if aw[metric] != "" and control[metric] != "" else None
            effects.append({"system": key[0], "filename": key[1], "time_ps": aw["time_ps"],
                            "tier": "H", "metric": metric, "contrast": "AW_minus_full_power_aw_radii",
                            "delta": effect, "interpretation": "same atoms/coordinates/radii; weighting and stated eligibility rules differ"})
        membership = read_csv(output / "cases" / Path(key[1]).stem / "H/selection_membership.csv")
        full = {json.dumps(r["stable_atom_ids"], separators=(",", ":")) for r in membership if r["full_power"] == "True"}
        for method in ("alpha0_power_before_M", "cazals_alpha0_power"):
            chosen = {json.dumps(r["stable_atom_ids"], separators=(",", ":")) for r in membership if r[method] == "True"}
            full_atoms = {tuple(a) for pair in full for a in json.loads(pair)}
            chosen_atoms = {tuple(a) for pair in chosen for a in json.loads(pair)}
            comparison = {"system": key[0], "filename": key[1], "selected_method": method,
                          "full_facets": len(full), "selected_facets": len(chosen),
                          "retained_fraction": len(chosen) / len(full) if full else None,
                          "shared_facets": len(full & chosen), "full_only_facets": len(full - chosen),
                          "selected_only_facets": len(chosen - full),
                          "removed_stable_facets": [json.loads(v) for v in sorted(full - chosen)],
                          "removed_interface_atoms": sorted(full_atoms - chosen_atoms),
                          "added_interface_atoms": sorted(chosen_atoms - full_atoms)}
            for metric in STUDY_METRICS:
                comparison[f"full_{metric}"] = case_rows["full_power"][metric]
                comparison[f"selected_{metric}"] = case_rows[method][metric]
            selections.append(comparison)
    write_csv(output / "power_vs_aw_effects.csv", effects)
    write_csv(output / "alpha0_vs_full_power.csv", selections)
    paired = []
    solv_file = output / "solvent_paired_case_metrics.csv"
    if solv_file.exists():
        solvent_rows = read_csv(solv_file)
        by_key = {(r["filename"], r["tier"], r["interface_kind"] if "interface_kind" in r else "protein_protein", r["method"]): r
                  for r in solvent_rows}
        for filename, tier, interface, method in sorted(by_key):
            if tier != "H+S":
                continue
            base_method = method.removesuffix("_protein_water") if interface == "protein_water" else method
            source_interface = "protein_water" if interface == "protein_water" else "protein_protein"
            solv = by_key[(filename, "H+S", source_interface, method)]
            hrow = by_key.get((filename, "H", "protein_protein", base_method))
            if hrow is None:
                continue
            for metric in STUDY_METRICS:
                a, b = solv.get(metric, ""), hrow.get(metric, "")
                if a == "" or b == "":
                    raise ValueError(f"Undefined paired solvent quantity {filename}/{method}/{metric}")
                paired.append({"system": "2KAI", "filename": filename, "time_ps": solv["time_ps"],
                               "tier": "H vs H+S", "interface_kind": interface, "method": method,
                               "metric": metric, "desolvated_H": b, "finite_explicit_solvent_HS": a,
                               "delta_solvent_HS_minus_H": float(a) - float(b),
                               "coordinate_identity": "same protein stable IDs, radii, and coordinates; waters/ions added",
                               "ions": "neutral surroundings; not in either interface partner group"})
    write_csv(output / "solvent_paired_effects.csv", paired)

    decomposition = []
    for system in sorted({r["system"] for r in static}):
        stat = {r["method"]: r for r in static if r["system"] == system}
        if not {"aw", "full_power_aw_radii", "full_power", "cazals_alpha0_power"} <= set(stat):
            raise ValueError(f"Static methods incomplete for {system}")
        for method in stat:
            row = stat[method]
            decomposition.append({"system": system, "period": "static",
                                  "method": method, **{k: row[k] for k in STUDY_METRICS},
                                  "geometry_effect_comparator": "aw_minus_full_power_aw_radii on matched atom radii"})
    for frame in trajectory:
        case_rows = [r for r in trajectory if r["filename"] == frame["filename"]]
        indexed = {r["method"]: r for r in case_rows}
        for method, row in indexed.items():
            decomposition.append({"system": "2KAI", "period": frame["filename"], "time_ps": row["time_ps"],
                                  "method": method, **{k: row[k] for k in STUDY_METRICS},
                                  "geometry_effect_AW_minus_matched_Power": float(indexed["aw"]["interface_area_A2"]) - float(indexed["full_power_aw_radii"]["interface_area_A2"]) if method in {"aw", "full_power_aw_radii"} else "",
                                  "MD_effect_AW_minus_crystal": float(row["combined_H_A"]) - float(next(r for r in static if r["system"] == "2KAI" and r["method"] == "aw")["combined_H_A"]) if method == "aw" else ""})
    write_csv(output / "method_effect_decomposition.csv", decomposition)


def create_visualization_exports(audit, output, rows, *, power_run=None, network=None,
                                 aw_result=None, control_power=None,
                                 control_polygons=None):
    """Export selected atom PDBs and sampled interface surfaces as PyMOL CGO."""
    from vorpy.src.analyze.power_interface_curvature import analyze_power_interface_curvature
    from vorpy.src.analyze.power_interface_selection import power_selection_views
    from vorpy.src.analyze.aw_interface_curvature import analyze_aw_interface_curvature
    from vorpy.src.analyze.power_interface_cli import run_power_interface_pdb
    from vorpy.src.analyze.power_vs_aw import _power_area
    from dataclasses import asdict

    output = Path(output)
    name = Path(audit["filename"]).stem
    case = output / "cases" / name / "H"
    source = Path(audit["path"]) if audit["kind"] == "static" else output / "atom_models" / name / "H.pdb"
    run = power_run or run_power_interface_pdb(
        source, set(audit["canonical_group_A_chains"]),
        set(audit["canonical_group_B_chains"]), case / "viz_power_replay")
    power, prepared = run["result"], run["prepared"]
    views, polygons = power_selection_views(power)
    network = network or build_or_load_aw(output / "atom_models" / name / "H.pdb", case)
    radii = radii_audit(network.sys, prepared, case / "viz_radii_audit.csv")
    mapping = {r["power_generator_id"]: r["aw_generator_id"] for r in radii}
    aw = aw_result or analyze_aw_interface_curvature(
        network, {mapping[i] for i in power.group_a},
        {mapping[i] for i in power.group_b}, selection_mode="supported")
    if control_power is None:
        controls = analyze_power_interface_curvature(prepared.points,
            [float(network.sys.balls.at[mapping[i], "rad"]) for i in range(len(prepared.atoms))],
            power.group_a, power.group_b)
        control_power, control_polygons = power_selection_views(controls)
        control_power = control_power["full_power"]
    else:
        # Geometry and facets were solved in run_geometry_case; reuse them.
        _, solved_control_polygons = power_selection_views(control_power)
        control_polygons = control_polygons or solved_control_polygons
    if control_power is not None and hasattr(control_power, "final_facets"):
        control_views = {"full_power": control_power}
    selection_by_method = {**views, "full_power_aw_radii": control_views["full_power"]}
    polys_by_method = {**{k: polygons for k in views}, "full_power_aw_radii": control_polygons}
    mapping_identity = {a.index: (a.chain, a.residue_number, a.insertion_code, a.residue_name, a.atom_name)
                        for a in prepared.atoms}
    export_root = output / "visualization" / name
    export_root.mkdir(parents=True, exist_ok=True)
    view_path = export_root / "protein_heavy.pdb"
    if not view_path.exists():
        shutil.copy2(source, view_path)
    pml = [f"load {view_path.resolve().as_posix()}, protein", "hide everything, protein",
           "show cartoon, protein and polymer", "orient protein", "view = cmd.get_view()", "bg_color white"]
    scripts = []
    methods = {**selection_by_method, "aw": aw}
    palette = {"full_power": (0.12, 0.38, 0.75), "alpha0_power_before_M": (0.5, 0.55, 0.6),
               "cazals_alpha0_power": (0.2, 0.65, 0.2), "full_power_aw_radii": (0.55, 0.35, 0.75),
               "aw": (0.9, 0.28, 0.12)}
    for method, result in methods.items():
        cgopath = export_root / f"{method}.json"
        surfaces = []
        if method == "aw":
            for surface in result.selected_surfaces:
                raw = network.surfs.loc[surface.surface_id].get("points")
                vertices = np.asarray(raw, dtype=float)
                if vertices.ndim == 3 and vertices.shape[-1] == 3:
                    vertices = vertices.reshape((-1, 3))
                if vertices.ndim != 2 or vertices.shape[1] != 3 or len(vertices) < 3:
                    continue
                # Render the existing sampled AW polygon; no new tessellation is performed.
                surfaces.append({"feature_id": surface.surface_id, "points": vertices.tolist(),
                                 "color": palette[method], "generator_ids": surface.generator_ids})
            atom_ids = sorted({i for s in result.selected_surfaces for i in s.generator_ids})
            records = [network.sys.balls.loc[i] for i in atom_ids]
            atoms = [{"chain": str(r["chain_name"]), "residue": int(r["auth_seq_id"] if str(r["auth_seq_id"]).strip() else r["res_seq"]),
                      "residue_name": str(r["res_name"]), "atom": str(r["name"]), "xyz": list(map(float, r["loc"]))}
                     for r in records]
        else:
            selected = {f.generator_ids for f in result.final_facets}
            for ids in sorted(selected):
                polygon = polys_by_method[method][ids]
                if polygon["status"] == "finite_bounded" and polygon.get("polygon_vertices"):
                    surfaces.append({"feature_id": ids, "points": polygon["polygon_vertices"],
                                     "color": palette[method], "generator_ids": [mapping_identity[i] for i in ids]})
            atoms_ids = sorted({i for pair in selected for i in pair})
            atoms = [{"chain": mapping_identity[i][0], "residue": mapping_identity[i][1],
                      "residue_name": mapping_identity[i][3], "atom": mapping_identity[i][4],
                      "xyz": prepared.points[i].tolist()} for i in atoms_ids]
        write_json(cgopath, {"method": method, "system": name, "surfaces": surfaces, "interface_atoms": atoms,
                             "atom_view_sha256": sha256(view_path),
                             "source_structure_sha256": audit["sha256"]})
        cmd_path = export_root / f"{method}_surface.py"
        cmd_path.write_text(
            "import json\nfrom pymol import cmd\nfrom pymol.cgo import BEGIN, TRIANGLES, COLOR, VERTEX, END\n"
            f"data=json.load(open(r'{cgopath.resolve()}'))\nobj=[]\n"
            "for s in data['surfaces']:\n"
            " pts=s['points']; c=s['color']\n"
            " if len(pts)<3: continue\n"
            " obj += [BEGIN,TRIANGLES,COLOR,*c]\n"
            " for i in range(1,len(pts)-1):\n"
            "  obj += [VERTEX,*pts[0],VERTEX,*pts[i],VERTEX,*pts[i+1]]\n"
            " obj += [END]\n"
            # The paired PML establishes and retains one camera before these
            # CGO objects are loaded; adding an object does not need a view var.
            f"cmd.load_cgo(obj, '{method}_interface')\n", encoding="utf-8")
        scripts.append(cmd_path)
    # A matched camera script loads the same coordinates once and toggles methods.
    (export_root / "paired_interfaces.pml").write_text(
        f"load {view_path.resolve().as_posix()}, protein\nhide everything, protein\nshow cartoon, protein and polymer\norient protein\n"
        + "\n".join(f"run {p.resolve().as_posix()}" for p in scripts) + "\n",
        encoding="utf-8")
    # Stable Lys15 patch export, using the validated 2KAI residue identity.
    if name == "2KAI":
        lys_ids = [i for i, identity in mapping_identity.items()
                   if identity[0] == "I" and identity[3] == "LYS" and identity[1] == 15]
        write_json(export_root / "lys15_localized_region.json", {
            "atom_identity": "chain I, LYS 15", "generator_ids": lys_ids,
            "power_facets": [{"generator_ids": list(f.generator_ids), "area_A2": polygons[f.generator_ids]["area_A2"]}
                             for f in selection_by_method["cazals_alpha0_power"].final_facets
                             if set(f.generator_ids) & set(lys_ids)],
            "aw_surfaces": [s.surface_id for s in aw.selected_surfaces if set(s.generator_ids) & set(mapping[i] for i in lys_ids)]})


def make_figures_from_csv(output):
    from vorpy.src.analyze.comparative_study_figures import make_figures
    return make_figures(output)


def summarize_study(output):
    output = Path(output)
    stat = read_csv(output / "static_method_comparison.csv")
    summary = ["Comparative geometry study", "", "OBSERVATIONS"]
    for system in ("2KAI", "1AK4", "1BRC"):
        rows = {r["method"]: r for r in stat if r["system"] == system}
        if rows:
            summary.append(f"{system}: static H atoms; methods={sorted(rows)}.")
            for method, values in sorted(rows.items()):
                summary.append(f"  {method}: N={values['atom_count']}, A={float(values['interface_area_A2']):.9g} A^2, "
                               f"signed/raw-edge={float(values['signed_raw_edge_A_rad']):.9g}, "
                               f"unsigned/raw-edge={float(values['unsigned_raw_edge_A_rad']):.9g} A rad, "
                               f"smooth={float(values['smooth_H_A']):.9g} A, "
                               f"C_H={float(values['combined_H_A']):.9g} A, components={values['connected_components']}.")
    traj = read_csv(output / "trajectory_summary.csv")
    if traj:
        summary.extend(["", "MD DESCRIPTIVES (time-correlated trajectory; no independent-sample significance test)"])
        for row in traj:
            if row.get("metric") in {"interface_area_A2", "combined_H_A", "unsigned_edge_per_area_rad_per_A"} and row.get("mean"):
                summary.append(f"  {row['method']} {row['metric']}: mean={float(row['mean']):.9g}; "
                               f"SD={float(row['sd']):.9g}; range=[{float(row['min']):.9g}, {float(row['max']):.9g}]; "
                               f"crystal={row['crystal']}; MD mean minus crystal={row['MD_mean_minus_crystal']}.")
    paired_path = output / "solvent_paired_effects.csv"
    if paired_path.exists():
        paired = read_csv(paired_path)
        summary.extend(["", "FINITE EXPLICIT-SOLVENT CONTEXT", f"Paired solvent effects: {len(paired)} rows."])
        for method in sorted({r["method"] for r in paired}):
            picked = [float(r["delta_solvent_HS_minus_H"]) for r in paired
                      if r["method"] == method and r["metric"] == "interface_area_A2"]
            if picked:
                summary.append(f"  {method} protein-protein area solvent deltas: "
                               f"mean={sum(picked)/len(picked):.9g} A^2, {len(picked)} trajectory frames.")
    summary.extend(["", "INTERPRETATIONS", "Geometry contrasts compare AW with full Power using the exact same AW-loaded radii.",
                    "Frozen Chothia Cazals Power versus frozen AW also differs in radius model; it is not a pure weighting-only contrast.",
                    "Full Power to alpha0_power_before_M isolates alpha selection; alpha0_power_before_M to cazals_alpha0_power isolates the M=5 filter.",
                    "A positive cancellation fraction means opposing signs cancel; it is 1-|signed|/unsigned.",
                    "AW selected smooth, edge, and combined signs follow the A->B interface normal.",
                    "Solvent results describe a finite explicit-solvent crop; convex-hull margin flags finite-cloud boundary sensitivity.",
                    "Power edge curvature covers interior seams shared by selected facets; open-boundary edges are not assigned a turning term.",
                    "Empty metric cells mean mathematically undefined with a recorded denominator status, not removed observations."])
    (output / "comparative_study_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")


def phase_run(output, data_dir, frames_dir, target_phase):
    output = Path(output)
    audits = read_csv(output / "input_audit.csv")
    for audit in audits:
        for field in ("group_A_chains", "group_B_chains", "canonical_group_A_chains",
                      "canonical_group_B_chains", "chain_map"):
            if isinstance(audit.get(field), str) and audit[field][:1] in "[{":
                audit[field] = json.loads(audit[field])
        if audit.get("time_ps") not in {None, ""}:
            audit["time_ps"] = float(audit["time_ps"])
    static = [r for r in audits if r.get("status") == "audited" and r.get("kind") == "static"]
    frames = sorted([r for r in audits if r.get("status") == "audited" and r.get("kind") == "MD_frame"],
                    key=lambda r: float(r["time_ps"]))
    crystal = next(r for r in static if r["system"] == "2KAI")
    # The supplied unlabelled centered copy is audited but is never relabeled by
    # row order; MD frames have complete identity-based partner mappings.
    phase_metrics = []
    tasks = []
    if target_phase >= 2:
        tasks.append((2, [crystal], True))
    if target_phase >= 3:
        tasks.append((3, [r for r in static if r["system"] in {"1AK4", "1BRC"}], False))
    if target_phase >= 4:
        tasks.append((4, frames, False))
    for phase, cases, frozen in tasks:
        for audit in cases:
            if audit.get("all_atoms_needed_for_group_matching_present") in {"False", False}:
                raise ValueError(f"Input matching failed; phase {phase} gated: {audit['filename']}")
            print(f"PHASE {phase} INPUT {audit['filename']}", flush=True)
            rows = run_geometry_case(audit, output, validate_frozen=frozen)
            phase_metrics.extend(rows)
            write_json(output / "phase_status.json", {"phase": phase, "last_completed_input": audit["filename"],
                                                       "status": "running", "validated_rows": len(phase_metrics),
                                                       "geometry_source_revision": MODULE_REVISION})
            aggregate_study(output)
    if target_phase >= 5:
        if not frames:
            raise ValueError("No 2KAI frame files are available for paired solvent study")
        reference, _ = read_pdb(crystal["path"])
        paired_rows, excluded, boundary = [], [], []
        for audit in frames:
            print(f"PHASE 5 INPUT {audit['filename']}", flush=True)
            data = solvated_geometry_case(audit, output, reference)
            paired_rows.extend(data[0]); excluded.extend(data[1]); boundary.extend(data[2])
            h_file = output / "cases" / Path(audit["filename"]).stem / "H/metrics.csv"
            baseline = read_csv(h_file)
            for filename, additions in (("solvent_paired_case_metrics.csv", data[0]),
                                        ("solvent_excluded_geometry.csv", data[1]),
                                        ("solvent_boundary_audit.csv", data[2])):
                destination = output / filename
                write_csv(destination, _merge_rows(read_csv(destination) if destination.exists() else [], additions))
            aggregate_study(output)
        boundary_rows = read_csv(output / "solvent_boundary_audit.csv")
        unsafe = [r for r in boundary_rows if r.get("feature_id") == "FRAME_BOUNDARY_SUMMARY" and r.get("boundary_status") != "interior_at_sampled_geometry"]
        write_json(output / "phase_status.json", {"phase": 5, "status": "complete with boundary flags" if unsafe else "complete",
                                                   "boundary_sensitive_frames": [r["filename"] for r in unsafe],
                                                   "protocol": "finite explicit-solvent context"})
    if target_phase >= 6:
        aggregate_study(output)
        make_figures_from_csv(output)
        static2 = next(r for r in static if r["system"] == "2KAI")
        representative = frames[0]
        summarize_study(output)
        write_json(output / "phase_status.json", {"phase": 6, "status": "complete",
            "frames": len(frames), "figure_source": "saved CSV data only",
            "representative_MD_frame": representative["filename"]})


def _merge_rows(first, second):
    unique = {}
    for row in [*first, *second]:
        unique[(row.get("filename"), row.get("tier"), row.get("method"), row.get("interface_kind"), row.get("feature_id"))] = row
    return list(unique.values())


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--audit-only", action="store_true", help="Audit and build documented atom views; default")
    mode.add_argument("--compute", action="store_true", help="Compute the selected phases")
    mode.add_argument("--plots-only", action="store_true", help="Regenerate plots and summary from saved CSV data")
    parser.add_argument("--phase", type=int, choices=range(1, 7), default=None,
                        help="Compute through phase 1..6; default is audit-only unless --compute is set (then 6)")
    parser.add_argument("--data-dir", type=Path, default=Path("vorpy/data"))
    parser.add_argument("--frames-dir", type=Path, default=Path("vorpy/data/2KAI_Frames"))
    parser.add_argument("--output-dir", type=Path, default=Path("output/comparative_study"))
    parser.add_argument("--systems", nargs="+", default=["2KAI", "1AK4", "1BRC"])
    args = parser.parse_args(argv)
    if args.plots_only:
        aggregate_study(args.output_dir)
        make_figures_from_csv(args.output_dir)
        summarize_study(args.output_dir)
        print(f"Plots regenerated from saved CSV files in {args.output_dir / 'figures'}")
        return
    audits, blockers = run_audit(args.data_dir, args.frames_dir, args.output_dir, args.systems)
    print(f"Input audit: {args.output_dir / 'input_audit.csv'} ({len(audits)} rows including missing inputs)")
    print(f"Phase 1 complete; {len(blockers)} documented blockers.")
    if args.compute:
        phase_run(args.output_dir, args.data_dir, args.frames_dir, args.phase or 6)


if __name__ == "__main__":
    main()
