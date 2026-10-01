"""Run every structure file in a frame directory with both AW and POW."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime
from pathlib import Path
import re
import subprocess
import sys
import time


REPO = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = Path(r"E:\Kalliklein\2KAI_gromacs\final\analysis_frames")
EXTENSIONS = {".pdb", ".gro", ".cif", ".mmcif"}


def natural_key(path):
    return [int(part) if part.isdigit() else part.casefold()
            for part in re.split(r"(\d+)", str(path))]


def build_command(source, mode, destination, boundary, max_vert, exports):
    command = [sys.executable, str(REPO / "vorpy"), str(source)]
    if source.suffix.lower() == ".pdb":
        command.append("--all-frames")
    command.extend(["--boundary", boundary, "-s", "nt", mode])
    if max_vert is not None:
        command.extend(["and", "mv", str(max_vert)])
    command.extend(["-e", "dir", str(destination), "and"])
    command.extend(["only", "logs"] if exports == "logs" else [exports])
    return command


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path, default=DEFAULT_INPUT,
                        help=f"Frame directory (default: {DEFAULT_INPUT}).")
    parser.add_argument("--output", type=Path, default=REPO / "output" / "2kai_analysis_frames",
                        help="Output parent; a fresh timestamped run folder is created.")
    parser.add_argument("--recursive", action="store_true", help="Also search subdirectories.")
    parser.add_argument("--boundary", choices=["none", "shell", "explicit"], default="none")
    parser.add_argument("--max-vert", type=float, help="Optional mv setting; otherwise use VorPy's default.")
    parser.add_argument("--exports", choices=["large", "medium", "small", "all", "logs"],
                        default="large", help="VorPy export preset (default: large).")
    parser.add_argument("--dry-run", action="store_true", help="List commands without running or writing files.")
    args = parser.parse_args(argv)
    source_dir = args.input.resolve()
    if not source_dir.is_dir():
        parser.error(f"Frame directory not found: {source_dir}. Check that the drive is mounted.")
    if args.max_vert is not None and not (0 < args.max_vert <= 5000):
        parser.error("--max-vert must be positive and no greater than 5000.")
    candidates = source_dir.rglob("*") if args.recursive else source_dir.iterdir()
    sources = sorted((p for p in candidates if p.is_file() and p.suffix.lower() in EXTENSIONS),
                     key=natural_key)
    if not sources:
        parser.error("No PDB, GRO, CIF, or mmCIF frame files found.")
    run_dir = args.output.resolve() / datetime.now().strftime("run_%Y%m%d_%H%M%S_%f")
    print(f"{len(sources)} files; {2*len(sources)} AW/POW jobs, run sequentially.", flush=True)
    print(f"Output: {run_dir}", flush=True)
    if not args.dry_run:
        run_dir.mkdir(parents=True, exist_ok=False)
    failures = 0
    manifest = None
    try:
        if not args.dry_run:
            manifest = (run_dir / "summary.csv").open("w", newline="", encoding="utf-8")
            writer = csv.writer(manifest)
            writer.writerow(["input", "mode", "status", "exit_code", "seconds", "output", "log", "command"])
            manifest.flush()
        for index, source in enumerate(sources):
            for mode_index, mode in enumerate(("aw", "pow")):
                # Preserve relative directories and extensions to avoid stem collisions.
                destination = run_dir / mode / source.relative_to(source_dir)
                command = build_command(source, mode, destination, args.boundary,
                                        args.max_vert, args.exports)
                display = subprocess.list2cmdline(command)
                print(f"[{2*index+mode_index+1}/{2*len(sources)}] {mode.upper()}: {source.name}", flush=True)
                if args.dry_run:
                    print(display)
                    continue
                destination.mkdir(parents=True, exist_ok=True)
                log_path = destination / "console.log"
                started = time.monotonic()
                status, exit_code = "failed", -1
                with log_path.open("w", encoding="utf-8") as log:
                    log.write(display + "\n\n")
                    log.flush()
                    try:
                        result = subprocess.run(command, cwd=REPO, stdin=subprocess.DEVNULL,
                                                stdout=log, stderr=subprocess.STDOUT, check=False)
                        exit_code = result.returncode
                        status = "ok" if exit_code == 0 else "failed"
                    except OSError as error:
                        log.write(f"Could not start VorPy: {error}\n")
                    except KeyboardInterrupt:
                        status = "interrupted"
                        raise
                    finally:
                        elapsed = time.monotonic()-started
                        writer.writerow([str(source), mode, status, exit_code, f"{elapsed:.2f}",
                                         str(destination), str(log_path), display])
                        manifest.flush()
                failures += status != "ok"
                print(f"  {status.upper()} ({elapsed:.1f}s); log: {log_path}", flush=True)
    except KeyboardInterrupt:
        print("\nStopped. Completed jobs and summary.csv are preserved.", file=sys.stderr)
        return 130
    finally:
        if manifest is not None:
            manifest.close()
    if not args.dry_run:
        print(f"Finished: {2*len(sources)-failures} succeeded; {failures} failed.")
        print(f"Summary: {run_dir / 'summary.csv'}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
