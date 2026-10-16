"""Read-only preflight for the four-run 1BNA-DB1976 benchmark.

This script parses the two PDB inputs, validates the dry/wet molecular pair,
checks the selected source universes without solving geometry, and verifies the
four CLI plans.  It never creates output directories or launches VorPy.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
DATA = ROOT / "vorpy" / "data"
MANIFEST = ROOT / "docs" / "benchmarks" / "1bna_db1976" / "benchmark_manifest.json"
DEFAULT_OUTPUT_ROOT = ROOT / "output" / "benchmarks" / "1bna_DB1976"
WATER = {"HOH", "WAT", "SOL", "H2O", "TIP3", "TIP3P", "TIP4", "TIP4P", "SPC", "SPCE"}
IONS = {
    "LI", "NA", "NA+", "SOD", "K", "K+", "POT", "RB", "CS", "MG", "MG2", "MG2+",
    "CA", "CA2", "CA2+", "SR", "BA", "ZN", "ZN2", "ZN2+", "CD", "FE", "FE2", "FE3",
    "MN", "CU", "CU1", "CU2", "CO", "NI", "AL", "F", "CL", "CL-", "BR", "I", "IOD",
    "SO4", "PO4", "NH4", "CLA", "ION", "OUT", "MN2",
}
EXPECTED = {
    "wet": {"atoms": 10494, "molecular_atoms": 805, "water_atoms": 9669, "ion_atoms": 20},
    "dry": {"atoms": 805, "molecular_atoms": 805, "water_atoms": 0, "ion_atoms": 0},
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def pdb_records(path: Path):
    records = []
    models = 0
    for raw in path.read_text(encoding="ascii", errors="replace").splitlines():
        tag = raw[:6].strip().upper()
        if tag == "MODEL":
            models += 1
        if tag not in {"ATOM", "HETATM"}:
            continue
        if len(raw) < 54:
            raise ValueError(f"short atom record in {path}: {raw!r}")
        records.append({
            "serial": int(raw[6:11]),
            "name": raw[12:16].strip(),
            "res_name": raw[17:20].strip().upper(),
            "chain": raw[21:22].strip(),
            "res_seq": raw[22:26].strip(),
            "xyz": tuple(float(raw[start:start + 8]) for start in (30, 38, 46)),
            "element": raw[76:78].strip().upper() if len(raw) >= 78 else "",
        })
    return records, models or 1


def _identity(record):
    return (record["serial"], record["name"], record["res_name"], record["chain"], record["res_seq"])


def source_counts(records):
    water = [row for row in records if row["res_name"] in WATER]
    ions = [row for row in records if row["res_name"] in IONS and row["res_name"] not in WATER]
    return {
        "atoms": len(records),
        "molecular_atoms": len(records) - len(water) - len(ions),
        "water_atoms": len(water),
        "ion_atoms": len(ions),
        "water_residue_count": len({(r["res_name"], r["res_seq"]) for r in water}),
    }


def check_groups(records):
    group_a = [r for r in records if r["chain"] in {"A", "B"}]
    group_b = [r for r in records if r["chain"] == "C" and r["res_name"] == "MOL" and r["res_seq"] == "25"]
    if len(group_a) != 758 or len(group_b) != 47:
        raise AssertionError(f"group counts are A+B={len(group_a)}, C/MOL 25={len(group_b)}")
    return {
        "group_a_selector": "c A and c B",
        "group_b_selector": "c C",
        "group_a_atoms": len(group_a),
        "group_b_atoms": len(group_b),
        "group_a_chains": ["A", "B"],
        "group_b_chain": "C",
        "group_b_residue": "MOL 25",
    }


def compare_pair(wet, dry):
    if len(wet) < len(dry):
        raise AssertionError("wet source has fewer molecular records than dry source")
    mismatches = []
    max_delta = 0.0
    for index, (left, right) in enumerate(zip(wet[: len(dry)], dry)):
        if _identity(left) != _identity(right):
            mismatches.append({"index": index, "wet": _identity(left), "dry": _identity(right)})
            continue
        delta = math.sqrt(sum((left["xyz"][axis] - right["xyz"][axis]) ** 2 for axis in range(3)))
        max_delta = max(max_delta, delta)
        if delta > 1e-6:
            mismatches.append({"index": index, "coordinate_delta_A": delta})
    return {"matched_atoms": len(dry) - len(mismatches), "max_coordinate_delta_A": max_delta, "mismatches": mismatches[:5]}


def _system_audit(path: Path, mode_name: str):
    from vorpy.src.boundary import SolventMode, prepare_network_geometry
    from vorpy.src.system.system import System

    mode = SolventMode.FULL_SOLVENT if mode_name == "wet" else SolventMode.MOLECULE_ONLY
    system = System(file=str(path), make_dir=False, print_actions=False, solvent_mode=mode)
    geometry = prepare_network_geometry(system, mode)
    active = [index for index in geometry.system_indices if index is not None]
    stable_ids = [str(value) for value in geometry.stable_ids]
    if len(active) != EXPECTED[mode_name]["atoms"]:
        raise AssertionError(f"{mode_name} active generator count is {len(active)}")
    if len(set(stable_ids)) != len(stable_ids) or not all(stable_ids):
        raise AssertionError(f"{mode_name} stable generator identities are not unique")

    balls = system.balls
    ligand = balls[
        (balls["chain_name"].astype(str).str.strip() == "C")
        & (balls["res_name"].astype(str).str.strip().str.upper() == "MOL")
        & (balls["res_seq"].astype(str).str.strip() == "25")
    ]
    if mode_name == "wet":
        selenium = ligand[ligand["element"].astype(str).str.strip().str.upper() == "SE"]
        if len(selenium) != 1 or float(selenium["rad"].iloc[0]) != 1.90:
            raise AssertionError("SE21 is not resolved to supported SE radius 1.90 A")
    inferred = []
    from vorpy.src.objects.atom import get_element
    for row in ligand.itertuples(index=False):
        expected_element = get_element(atom_name=str(row.name))
        if str(row.element).upper() != expected_element.upper():
            inferred.append((row.name, row.element, expected_element))
        if not math.isfinite(float(row.rad)) or float(row.rad) <= 0:
            raise AssertionError(f"invalid radius for ligand atom {row.name}")
    if inferred:
        raise AssertionError(f"ligand element/radius inference mismatch: {inferred[:3]}")

    import numpy as np
    group_a = balls[balls["chain_name"].astype(str).str.strip().isin(["A", "B"])]
    coords_a = np.asarray(group_a["loc"].tolist(), dtype=float)
    coords_b = np.asarray(ligand["loc"].tolist(), dtype=float)
    distances = np.linalg.norm(coords_a[:, None, :] - coords_b[None, :, :], axis=2)
    min_distance = float(distances.min())
    radius_gap = min_distance - float(group_a.iloc[int(np.unravel_index(distances.argmin(), distances.shape)[0])]["rad"]) - float(ligand.iloc[int(np.unravel_index(distances.argmin(), distances.shape)[1])]["rad"])
    return {
        "active_generators": len(active),
        "stable_generator_ids": len(stable_ids),
        "selenium_radius_A": 1.90 if mode_name == "wet" else None,
        "ligand_element_fallback_compatible": True,
        "min_group_center_distance_A": min_distance,
        "min_group_vdw_gap_A": radius_gap,
        "source_mode": mode.value,
    }


def _powershell_processes():
    if os.name != "nt":
        return []
    command = (
        "Get-CimInstance Win32_Process | "
        "Where-Object { $_.CommandLine -and "
        "($_.CommandLine -match '(?i)(-m\\s+vorpy\\b|(?:^|\\s)vorpy(?:\\.exe)?\\s+)') } | "
        "Select-Object ProcessId,Name,CommandLine | ConvertTo-Json -Compress"
    )
    completed = subprocess.run(
        ["powershell.exe", "-NoProfile", "-Command", command],
        capture_output=True, text=True, check=False,
    )
    if not completed.stdout.strip():
        return []
    data = json.loads(completed.stdout)
    return data if isinstance(data, list) else [data]


def resource_snapshot():
    if os.name == "nt":
        command = "[math]::Round((Get-CimInstance Win32_LogicalDisk -Filter \"DeviceID='C:'\").FreeSpace / 1GB, 2)"
        disk = subprocess.run(["powershell.exe", "-NoProfile", "-Command", command], capture_output=True, text=True, check=False).stdout.strip()
        command = "[math]::Round((Get-Counter '\\Memory\\Available Bytes').CounterSamples[0].CookedValue / 1GB, 2)"
        ram = subprocess.run(["powershell.exe", "-NoProfile", "-Command", command], capture_output=True, text=True, check=False).stdout.strip()
    else:
        usage = shutil.disk_usage(ROOT)
        disk = round(usage.free / 1024**3, 2)
        ram = round(os.sysconf("SC_AVPHYS_PAGES") * os.sysconf("SC_PAGE_SIZE") / 1024**3, 2)
    return {"available_ram_GiB": float(ram), "free_disk_GiB": float(disk)}


def run_specs(output_root=DEFAULT_OUTPUT_ROOT):
    wet = (DATA / "1bna_DB1976.pdb").resolve().as_posix()
    dry = (DATA / "1bna_DB1976_no_sol.pdb").resolve().as_posix()
    root = Path(output_root).expanduser().resolve().as_posix()
    return {
        "aw_wet": ["--fail-if-exists", wet, "-g", "c A and c B", "-g", "c C", "-i", "-c", "system", "-s", "nt", "aw", "and", "sr", "0.2", "and", "mv", "5", "--boundary", "none", "--solvent-mode", "FULL_SOLVENT", "--profile", "analysis", "--timing-report", f"{root}/aw_wet/timing.json", "-e", "dir", f"{root}/aw_wet", "--save-network", f"{root}/aw_wet/1bna_DB1976/1bna_DB1976_aw_wet.vpy"],
        "pow_wet": ["--fail-if-exists", wet, "-g", "c A and c B", "-g", "c C", "-i", "-c", "system", "-s", "nt", "pow", "and", "sr", "0.2", "and", "mv", "5", "--boundary", "none", "--solvent-mode", "FULL_SOLVENT", "--profile", "analysis", "--timing-report", f"{root}/pow_wet/timing.json", "-e", "dir", f"{root}/pow_wet", "--save-network", f"{root}/pow_wet/1bna_DB1976/1bna_DB1976_pow_wet.vpy"],
        "aw_dry": ["--fail-if-exists", dry, "-g", "c A and c B", "-g", "c C", "-i", "-c", "system", "-s", "nt", "aw", "and", "sr", "0.2", "and", "mv", "5", "--boundary", "none", "--solvent-mode", "MOLECULE_ONLY", "--profile", "analysis", "--timing-report", f"{root}/aw_dry/timing.json", "-e", "dir", f"{root}/aw_dry", "--save-network", f"{root}/aw_dry/1bna_DB1976_no_sol/1bna_DB1976_aw_dry.vpy"],
        "pow_dry": ["--fail-if-exists", dry, "-g", "c A and c B", "-g", "c C", "-i", "-c", "system", "-s", "nt", "pow", "and", "sr", "0.2", "and", "mv", "5", "--boundary", "none", "--solvent-mode", "MOLECULE_ONLY", "--profile", "analysis", "--timing-report", f"{root}/pow_dry/timing.json", "-e", "dir", f"{root}/pow_dry", "--save-network", f"{root}/pow_dry/1bna_DB1976_no_sol/1bna_DB1976_pow_dry.vpy"],
    }


def verify_cli_specs(specs):
    from vorpy.src.command.calculation_plan import build_calculation_plan
    from vorpy.src.command.vpy_cmnd import boundary_cli_options, execution_profile_cli_options, fail_if_exists_cli_options, solvent_cli_options

    verified = {}
    for name, argv in specs.items():
        args = list(argv)
        args, profile = execution_profile_cli_options(args)
        args, fail_if_exists = fail_if_exists_cli_options(args)
        args, boundary, boundary_explicit = boundary_cli_options(args)
        args, solvent_mode = solvent_cli_options(args)
        input_source = args[0]
        plan = build_calculation_plan(args, input_source=input_source, analysis_profile=profile)
        if not fail_if_exists or not boundary_explicit or boundary.mode != "none":
            raise AssertionError(f"{name}: protection/boundary options did not parse as required")
        if not plan.interface_requests or not plan.full_system_request:
            raise AssertionError(f"{name}: missing one-source interface plan")
        if not plan.solver_settings or not any(setting[:2] == ("nt", name.split("_")[0]) for setting in plan.solver_settings):
            raise AssertionError(f"{name}: method setting did not parse")
        expected_mode = "FULL_SOLVENT" if name.endswith("wet") else "MOLECULE_ONLY"
        if getattr(solvent_mode, "value", solvent_mode) != expected_mode:
            raise AssertionError(f"{name}: wrong solvent mode")
        verified[name] = {
            "input": input_source,
            "profile": profile,
            "method": name.split("_")[0],
            "solvent_mode": expected_mode,
            "boundary": boundary.mode,
            "full_system_request": plan.full_system_request,
            "interface_request_count": len(plan.interface_requests),
        }
    return verified


def audit(manifest_path=MANIFEST, output_root=DEFAULT_OUTPUT_ROOT, repeat_id=None):
    wet_path = DATA / "1bna_DB1976.pdb"
    dry_path = DATA / "1bna_DB1976_no_sol.pdb"
    wet_records, wet_frames = pdb_records(wet_path)
    dry_records, dry_frames = pdb_records(dry_path)
    if wet_frames != 1 or dry_frames != 1:
        raise AssertionError(f"expected one frame, got wet={wet_frames}, dry={dry_frames}")
    for label, records in (("wet", wet_records), ("dry", dry_records)):
        serials = [record["serial"] for record in records]
        if serials != list(range(1, len(records) + 1)):
            raise AssertionError(f"{label} atom serial IDs are not stable 1..N")
    wet_counts = source_counts(wet_records)
    dry_counts = source_counts(dry_records)
    if wet_counts != {**EXPECTED["wet"], "water_residue_count": 3223}:
        raise AssertionError(f"wet composition mismatch: {wet_counts}")
    if dry_counts != {**EXPECTED["dry"], "water_residue_count": 0}:
        raise AssertionError(f"dry composition mismatch: {dry_counts}")
    groups = check_groups(dry_records)
    pair = compare_pair(wet_records, dry_records)
    if pair["mismatches"]:
        raise AssertionError(f"dry/wet molecular correspondence failed: {pair['mismatches']}")
    systems = {"wet": _system_audit(wet_path, "wet"), "dry": _system_audit(dry_path, "dry")}
    output_root = Path(output_root).expanduser().resolve()
    specs = run_specs(output_root)
    cli = verify_cli_specs(specs)
    roots = {name: output_root / name for name in specs}
    collisions = [str(path) for path in roots.values() if path.exists()]
    archives = [str(Path(spec[-1])) for spec in specs.values() if Path(spec[-1]).exists()]
    return {
        "status": "READY_PENDING_RESOURCE_GATE" if not collisions else "BLOCKED_OUTPUT_COLLISION",
        "repeat_id": repeat_id,
        "output_root": str(output_root),
        "inputs": {
            "wet": {"path": wet_path.as_posix(), "sha256": sha256(wet_path), "frames": wet_frames, **wet_counts},
            "dry": {"path": dry_path.as_posix(), "sha256": sha256(dry_path), "frames": dry_frames, **dry_counts},
        },
        "groups": groups,
        "matched_molecular_pair": pair,
        "generator_audit": systems,
        "cli": cli,
        "output_collisions": collisions + archives,
        "resources": resource_snapshot(),
        "competing_vorpy_processes": _powershell_processes(),
        "manifest": str(Path(manifest_path).resolve()),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument(
        "--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT,
        help="fresh root containing the four run directories; historical roots are not consulted",
    )
    parser.add_argument("--repeat-id", default=None, help="optional repetition identity to echo in JSON")
    parser.add_argument("--json", action="store_true", help="emit machine-readable JSON")
    args = parser.parse_args(argv)
    try:
        result = audit(args.manifest, args.output_root, args.repeat_id)
    except Exception as error:
        print(f"PREFLIGHT BLOCKED: {type(error).__name__}: {error}", file=sys.stderr)
        return 2
    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print("1BNA-DB1976 four-run preflight")
        if result.get("repeat_id"):
            print(f"Repetition: {result['repeat_id']}")
        print(f"Output root: {result['output_root']}")
        print(f"Status: {result['status']}")
        for label in ("wet", "dry"):
            item = result["inputs"][label]
            print(f"{label}: {item['atoms']:,} atoms, SHA-256 {item['sha256']}")
        print(f"Matched dry/wet molecular atoms: {result['matched_molecular_pair']['matched_atoms']:,}; max coordinate delta {result['matched_molecular_pair']['max_coordinate_delta_A']:.6g} A")
        print(f"Groups: A+B {result['groups']['group_a_atoms']:,}; C/MOL 25 {result['groups']['group_b_atoms']:,}")
        for name, item in result["generator_audit"].items():
            print(f"{name} source: {item['active_generators']:,} active generators; min group distance {item['min_group_center_distance_A']:.3f} A; vdw gap {item['min_group_vdw_gap_A']:.3f} A")
        print(f"Resources: {result['resources']['available_ram_GiB']:.2f} GiB RAM, {result['resources']['free_disk_GiB']:.2f} GiB disk")
        print(f"Competing VorPy processes: {len(result['competing_vorpy_processes'])}")
        print("Commands: four parser-valid, fail-if-exists protected, one-source interface plans")
    return 0 if result["status"] == "READY_PENDING_RESOURCE_GATE" else 2


if __name__ == "__main__":
    raise SystemExit(main())
