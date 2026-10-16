"""Read-only validator for completed 1BNA-DB1976 interface runs.

Pass the four run roots after manual execution.  The validator never loads a
VorPy archive, rebuilds geometry, changes outputs, or treats missing metrics as
zero.  It prints a comparison table and returns nonzero for incomplete or
unresolved runs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "docs" / "benchmarks" / "1bna_db1976" / "benchmark_manifest.json"
MAX_READ_BYTES = 256 * 1024 * 1024
REQUIRED = (
    "interface_surfs.off", "interface_edges.off", "interface_verts.off",
    "interface_info.txt", "interface_surfaces.csv", "interface_edges.csv",
    "interface_vertices.csv", "group_a_atoms.pdb", "group_b_atoms.pdb",
    "surr_atoms.pdb", "surr_residues.pdb",
)


def _path(value):
    path = Path(value).expanduser()
    return path if path.is_absolute() else (ROOT / path)


def _interface_dir(run_root: Path):
    if (run_root / "interface_info.txt").is_file():
        return run_root
    candidates = []
    for info in run_root.rglob("interface_info.txt"):
        try:
            depth = len(info.parent.relative_to(run_root).parts)
        except ValueError:
            continue
        if depth <= 4:
            candidates.append(info.parent)
    return sorted(candidates)[0] if len(candidates) == 1 else None


def _off(path: Path):
    result = {"path": str(path), "valid": False, "vertices": None, "faces": None, "edges": None, "reason": None}
    if not path.is_file():
        result["reason"] = "missing"
        return result
    if path.stat().st_size > MAX_READ_BYTES:
        result["reason"] = f"refused to read mesh larger than {MAX_READ_BYTES // 1024**2} MiB"
        return result
    try:
        lines = path.read_text(encoding="utf-8", errors="strict").splitlines()
        if not lines or lines[0].strip() != "OFF":
            raise ValueError("missing OFF header")
        index = 1
        while index < len(lines) and (not lines[index].strip() or lines[index].lstrip().startswith("#")):
            index += 1
        counts = [int(value) for value in lines[index].split()[:3]]
        if len(counts) != 3 or min(counts) < 0:
            raise ValueError("invalid OFF counts")
        n_vertices, n_faces, n_edges = counts
        index += 1
        vertices = []
        for _ in range(n_vertices):
            values = [float(value) for value in lines[index].split()[:3]]
            if len(values) != 3 or not all(math.isfinite(value) for value in values):
                raise ValueError("nonfinite OFF vertex")
            vertices.append(values)
            index += 1
        for _ in range(n_faces):
            values = [int(value) for value in lines[index].split()]
            if not values or values[0] < 3 or len(values) < values[0] + 1:
                raise ValueError("invalid OFF face")
            if any(value < 0 or value >= n_vertices for value in values[1:values[0] + 1]):
                raise ValueError("OFF face index outside vertex range")
            index += 1
        result.update(valid=True, vertices=n_vertices, faces=n_faces, edges=n_edges)
    except (OSError, ValueError, IndexError, UnicodeError) as error:
        result["reason"] = str(error)
    return result


def _csv_table(path: Path, required_columns=()):
    result = {"path": str(path), "valid": False, "rows": None, "reason": None, "area_A2": None}
    if not path.is_file():
        result["reason"] = "missing"
        return result
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            fields = reader.fieldnames or []
            missing = [field for field in required_columns if field not in fields]
            if missing:
                raise ValueError(f"missing columns: {missing}")
            rows = list(reader)
        area = None
        area_field = next((field for field in fields if "Surface Area" in field), None)
        if area_field is not None:
            values = []
            for row in rows:
                raw = (row.get(area_field) or "").strip()
                if raw:
                    values.append(float(raw))
            if values:
                area = sum(values)
        result.update(valid=True, rows=len(rows), area_A2=area)
    except (OSError, ValueError, UnicodeError, csv.Error) as error:
        result["reason"] = str(error)
    return result


def _pdb_count(path: Path):
    if not path.is_file():
        return {"present": False, "atoms": None, "reason": "missing"}
    if path.stat().st_size > MAX_READ_BYTES:
        return {"present": False, "atoms": None, "reason": "PDB exceeds read guard"}
    try:
        atoms = sum(1 for line in path.read_text(encoding="ascii", errors="strict").splitlines() if line[:6].strip() in {"ATOM", "HETATM"})
    except (OSError, UnicodeError) as error:
        return {"present": False, "atoms": None, "reason": str(error)}
    return {"present": True, "atoms": atoms, "reason": None}


def _info(path: Path):
    result = {"present": False, "reason": "missing", "text": ""}
    if not path.is_file():
        return result
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as error:
        result["reason"] = str(error)
        return result
    result.update(present=True, reason=None, text=text)
    patterns = {
        "solvent_mode": r"^\s*Solvent mode:\s*(\S+)",
        "source_generators": r"^\s*Source generators:\s*(\d+)",
        "group_a_atoms": r"^\s*Group 1 atoms:\s*(\d+)",
        "group_b_atoms": r"^\s*Group 2 atoms:\s*(\d+)",
        "vertices": r"^\s*Vertices:\s*(\d+)",
        "edges": r"^\s*Edges:\s*(\d+)",
        "surfaces": r"^\s*Surfaces:\s*(\d+)",
        "direct_surfaces": r"^\s*Direct Group 1.*surfaces:\s*(\d+)",
        "area_A2": r"^\s*Interface area \([^)]*\):\s*([-+0-9.eE]+|Not available)",
        "completeness": r"^\s*Completeness:\s*(.+)$",
    }
    for key, pattern in patterns.items():
        match = re.search(pattern, text, re.IGNORECASE | re.MULTILINE)
        if match:
            value = match.group(1).strip()
            result[key] = float(value) if key == "area_A2" and value.lower() != "not available" else int(value) if key not in {"solvent_mode", "completeness", "area_A2"} else value
        else:
            result[key] = None
    result["curvature_status"] = "NOT_CALCULATED" if "NOT_CALCULATED" in text.upper() else None
    return result


def _runtime(run_root: Path):
    phase_timings = {}
    archive_seconds = None
    timing_files = []
    wall_seconds = None
    cpu_seconds = None
    peak_rss_bytes = None
    memory_source = None
    candidates = sorted(set(run_root.rglob("*.timing.json")) | set(run_root.rglob("summary.json")))
    for path in candidates:
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError):
            continue
        timing_files.append(str(path))
        if isinstance(payload.get("wall_seconds"), (int, float)):
            wall_seconds = float(payload["wall_seconds"])
        if isinstance(payload.get("cpu_seconds"), (int, float)):
            cpu_seconds = float(payload["cpu_seconds"])
        if isinstance(payload.get("peak_rss_bytes"), (int, float)):
            peak_rss_bytes = int(payload["peak_rss_bytes"])
        if payload.get("memory_source"):
            memory_source = payload["memory_source"]
        timing = payload.get("timing", {}) or {}
        total = timing.get("total")
        if isinstance(total, dict):
            if isinstance(total.get("wall_seconds"), (int, float)):
                wall_seconds = float(total["wall_seconds"])
            if isinstance(total.get("cpu_seconds"), (int, float)):
                cpu_seconds = float(total["cpu_seconds"])
        else:
            if isinstance(total, (int, float)):
                archive_seconds = float(total)
        if isinstance(timing, dict):
            for key in ("build", "solve", "source", "analysis", "export", "archive"):
                value = timing.get(key)
                if isinstance(value, (int, float)):
                    phase_timings[key] = float(value)
                elif isinstance(value, dict) and isinstance(value.get("wall_seconds"), (int, float)):
                    phase_timings[key] = float(value["wall_seconds"])
        phases = payload.get("phases", {})
        if isinstance(phases, dict):
            for key, value in phases.items():
                if isinstance(value, (int, float)):
                    phase_timings[key] = float(value)
                elif isinstance(value, dict) and isinstance(value.get("wall_seconds"), (int, float)):
                    phase_timings[key] = float(value["wall_seconds"])
    return {
        "source": timing_files[0] if timing_files else None,
        "wall_seconds": wall_seconds,
        "cpu_seconds": cpu_seconds,
        "archive_seconds": archive_seconds,
        "phase_timings": phase_timings,
        "peak_rss_bytes": peak_rss_bytes,
        "memory_source": memory_source,
        "timing_files": timing_files,
    }


def _archive_audit(run_root: Path, spec, manifest):
    """Read archive metadata only; never deserialize arrays or solve geometry."""
    archives = sorted(
        path for path in run_root.rglob("*.vpy")
        if path.is_file() and not path.name.endswith(".geometry.vpy")
    )
    result = {
        "path": str(archives[0]) if archives else None,
        "sha256": None,
        "bytes": None,
        "members": None,
        "metadata": {},
        "source_sha256": None,
        "method": None,
        "solvent_mode": None,
        "errors": [],
    }
    if not archives:
        result["errors"].append("no reusable non-checkpoint .vpy archive")
        return result
    archive = archives[0]
    digest = hashlib.sha256()
    with archive.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    result["sha256"] = digest.hexdigest().upper()
    result["bytes"] = archive.stat().st_size
    try:
        with zipfile.ZipFile(archive) as bundle:
            result["members"] = len(bundle.infolist())
            result["metadata"] = json.loads(bundle.read("metadata.json"))
    except (OSError, ValueError, KeyError, zipfile.BadZipFile, json.JSONDecodeError) as error:
        result["errors"].append(f"archive metadata unreadable: {error}")
        return result
    metadata = result["metadata"]
    result["source_sha256"] = str(metadata.get("source_sha256") or "").upper() or None
    result["method"] = str(metadata.get("scheme") or "").lower() or None
    settings = metadata.get("solve_settings", {})
    if isinstance(settings, dict):
        for key, value in settings.get("map", []):
            if key == "solvent_mode":
                result["solvent_mode"] = str(value)
    expected_input = manifest.get("inputs", {}).get(spec.get("input"), {})
    expected_hash = str(expected_input.get("sha256") or "").upper() or None
    if result["source_sha256"] != expected_hash:
        result["errors"].append(
            f"source hash {result['source_sha256']!r} != expected {expected_hash!r}"
        )
    if result["method"] != str(spec.get("method") or "").lower():
        result["errors"].append(
            f"archive method {result['method']!r} != expected {spec.get('method')!r}"
        )
    if result["solvent_mode"] != spec.get("source_mode"):
        result["errors"].append(
            f"archive solvent mode {result['solvent_mode']!r} != expected {spec.get('source_mode')!r}"
        )
    return result


def validate_run(run_id, run_root: Path, spec, manifest):
    interface = _interface_dir(run_root)
    result = {
        "run": run_id,
        "run_root": str(run_root),
        "interface_dir": str(interface) if interface else None,
        "method": spec.get("method"),
        "source_mode": spec.get("source_mode"),
        "settings_fingerprint": spec.get("settings_fingerprint"),
        "revision": manifest.get("worktree", {}).get("head"),
        "missing": [],
        "warnings": [],
    }
    if interface is None:
        result.update(status="UNRESOLVED", reason="could not identify exactly one interface directory")
        result["archive"] = _archive_audit(run_root, spec, manifest)
        result["runtime"] = _runtime(run_root)
        return result
    result["files"] = {name: (interface / name).is_file() for name in REQUIRED}
    result["missing"] = [name for name, present in result["files"].items() if not present]
    result["off"] = {name: _off(interface / name) for name in ("interface_surfs.off", "interface_edges.off", "interface_verts.off")}
    result["csv"] = {
        "surfaces": _csv_table(interface / "interface_surfaces.csv", ("Surface ID", "Surface Area (Å²)")),
        "edges": _csv_table(interface / "interface_edges.csv", ("Edge ID",)),
        "vertices": _csv_table(interface / "interface_vertices.csv", ("Vertex ID",)),
    }
    result["pdb"] = {name[:-4]: _pdb_count(interface / name) for name in ("group_a_atoms.pdb", "group_b_atoms.pdb", "surr_atoms.pdb", "surr_residues.pdb")}
    result["info"] = _info(interface / "interface_info.txt")
    result["archive"] = _archive_audit(run_root, spec, manifest)
    result["archive_paths"] = [result["archive"]["path"]] if result["archive"]["path"] else []
    result["runtime"] = _runtime(run_root)
    result["forbidden_geometry"] = [str(path) for path in interface.rglob("*.geometry.vpy")]
    result["forbidden_visualization"] = [str(path) for path in interface.rglob("visualization_metadata.json")]
    result["forbidden_archives"] = [str(path) for path in interface.rglob("*.vpy")]

    blockers = list(result["missing"])
    blockers.extend(f"{name}: {item['reason']}" for name, item in result["off"].items() if not item["valid"])
    blockers.extend(f"{name}: {item['reason']}" for name, item in result["csv"].items() if not item["valid"])
    if not result["info"].get("present"):
        blockers.append("interface_info.txt unreadable")
    blockers.extend(result["archive"].get("errors", []))
    if result["forbidden_geometry"] or result["forbidden_visualization"] or result["forbidden_archives"]:
        blockers.append("unrelated archive/visualization artifact inside interface directory")
    info = result["info"]
    expected_mode = spec["source_mode"]
    if info.get("solvent_mode") != expected_mode:
        blockers.append(f"source mode is {info.get('solvent_mode')!r}, expected {expected_mode}")
    expected_generators = 10494 if expected_mode == "FULL_SOLVENT" else 805
    if info.get("source_generators") != expected_generators:
        blockers.append(f"source generator count is {info.get('source_generators')!r}, expected {expected_generators}")
    if info.get("group_a_atoms") != 758 or info.get("group_b_atoms") != 47:
        blockers.append("group membership counts are not 758/47")
    csv_area = result["csv"]["surfaces"].get("area_A2")
    info_area = info.get("area_A2")
    if csv_area is None or info_area is None:
        result["warnings"].append("interface area is unavailable; no zero substituted")
    elif not math.isclose(csv_area, info_area, rel_tol=1e-6, abs_tol=1e-3):
        blockers.append(f"CSV/info area mismatch: {csv_area} vs {info_area} A^2")
    for name, item in result["pdb"].items():
        if not item["present"]:
            blockers.append(f"{name} PDB unavailable")
    result["topology"] = {key: info.get(key) for key in ("vertices", "edges", "surfaces", "direct_surfaces")}
    result["area_A2"] = info_area
    result["completeness"] = info.get("completeness")
    result["curvature_status"] = info.get("curvature_status")
    if result["completeness"] not in {"COMPLETE", "complete"}:
        result["warnings"].append(f"completeness remains {result['completeness']!r}")
    if result["curvature_status"] != "CALCULATED":
        result["warnings"].append("curvature is not certified by this interface-profile run")
    result["status"] = "UNRESOLVED" if blockers else "PARTIAL" if result["warnings"] else "VALIDATED"
    result["reason"] = "; ".join(blockers or result["warnings"] or ["all requested checks passed"])
    return result


def compare_table(results):
    headers = [item["run"] for item in results]
    rows = [
        ("status", [item.get("status") for item in results]),
        ("method", [item.get("method") for item in results]),
        ("revision", [item.get("revision") for item in results]),
        ("settings fingerprint", [item.get("settings_fingerprint") for item in results]),
        ("interface directory", [item.get("interface_dir") for item in results]),
        ("source mode", [item.get("info", {}).get("solvent_mode") for item in results]),
        ("source generators", [item.get("info", {}).get("source_generators") for item in results]),
        ("archive source SHA-256", [item.get("archive", {}).get("source_sha256") for item in results]),
        ("archive SHA-256", [item.get("archive", {}).get("sha256") for item in results]),
        ("topology V/E/F", ["{}/{}/{}".format(item.get("topology", {}).get("vertices"), item.get("topology", {}).get("edges"), item.get("topology", {}).get("surfaces")) if item.get("topology") else None for item in results]),
        ("direct surfaces", [item.get("topology", {}).get("direct_surfaces") for item in results]),
        ("unique interface area (A^2)", [item.get("area_A2") for item in results]),
        ("completeness", [item.get("completeness") for item in results]),
        ("curvature status", [item.get("curvature_status") for item in results]),
        ("runtime seconds", [item.get("runtime", {}).get("wall_seconds") for item in results]),
        ("archive seconds", [item.get("runtime", {}).get("archive_seconds") for item in results]),
        ("phase timings", [item.get("runtime", {}).get("phase_timings") or "UNAVAILABLE" for item in results]),
        ("peak RSS bytes", [item.get("runtime", {}).get("peak_rss_bytes") for item in results]),
        ("reason", [item.get("reason") for item in results]),
    ]
    lines = ["| Metric | " + " | ".join(headers) + " |", "|---|" + "---|" * len(headers)]
    lines.extend("| " + key + " | " + " | ".join("" if value is None else str(value).replace("|", "\\|") for value in values) + " |" for key, values in rows)
    return "\n".join(lines)


def historical_table(results, manifest):
    baselines = manifest.get("historical_comparison", {})
    rows = []
    for run_id in ("pow_wet", "aw_wet"):
        key = f"{run_id}_wall_seconds"
        baseline = baselines.get(key)
        repeat = next((item.get("runtime", {}).get("wall_seconds") for item in results if item.get("run") == run_id), None)
        delta = (repeat - baseline) if repeat is not None and baseline is not None else None
        percent = (100.0 * delta / baseline) if delta is not None and baseline else None
        rows.append((run_id, baseline, repeat, delta, percent,
                     "AVAILABLE" if repeat is not None else "UNAVAILABLE"))
    lines = ["Historical wet baseline comparison (phase timings are independent and remain unavailable unless recorded):",
             "| Run | Historical wall s | Repeat wall s | Delta s | Delta % | Status |",
             "|---|---:|---:|---:|---:|---|"]
    for run_id, baseline, repeat, delta, percent, status in rows:
        values = [run_id, baseline, repeat, delta, None if percent is None else round(percent, 3), status]
        lines.append("| " + " | ".join("UNAVAILABLE" if value is None else str(value) for value in values) + " |")
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    for name in ("aw-wet", "pow-wet", "aw-dry", "pow-dry"):
        parser.add_argument(f"--{name}", required=False, type=Path)
    parser.add_argument("--output-root", type=Path, help="infer four run roots beneath this repetition root")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(argv)
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    names = ("aw-wet", "pow-wet", "aw-dry", "pow-dry")
    if args.output_root:
        root = _path(args.output_root)
        roots = {name.replace("-", "_"): root / name.replace("-", "_") for name in names}
    else:
        missing = [name for name in names if getattr(args, name.replace("-", "_")) is None]
        if missing:
            parser.error("provide --output-root or all four run paths; missing " + ", ".join(missing))
        roots = {name.replace("-", "_"): _path(getattr(args, name.replace("-", "_"))) for name in names}
    results = [validate_run(run_id, roots[run_id], manifest["runs"][run_id], manifest) for run_id in ("aw_wet", "pow_wet", "aw_dry", "pow_dry")]
    if args.json:
        print(json.dumps({"results": results, "historical": manifest.get("historical_comparison", {})}, indent=2, sort_keys=True))
    else:
        print(compare_table(results))
        print()
        print(historical_table(results, manifest))
    return 0 if all(item.get("status") == "VALIDATED" for item in results) else 2


if __name__ == "__main__":
    raise SystemExit(main())
