#!/usr/bin/env python3
"""Benchmark vertex search with stable 5-A defaults and detailed counters."""
import argparse
import contextlib
import io
import json
import signal
import time
from pathlib import Path

from vorpy.src.network import Network
from vorpy.src.network import fast
from vorpy.workbench.services.structure_loader import load_pdb


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", type=Path)
    parser.add_argument("--max-vert", type=float, default=5.0)
    parser.add_argument("--residues", type=int, default=10, help="first N residue records; 0 means all atoms")
    parser.add_argument("--full-context", action="store_true", help="keep all atoms as spatial context while selecting the first residues")
    parser.add_argument("--network", choices=("aw", "pow", "prm", "all"), default="all")
    parser.add_argument("--timeout", type=int, default=0, help="seconds per network; 0 means unlimited")
    args = parser.parse_args()

    atoms = load_pdb(args.pdb).atoms
    residue_ids = []
    for atom in atoms:
        if atom.residue_sequence not in residue_ids:
            residue_ids.append(atom.residue_sequence)
    selected = set(residue_ids[:args.residues]) if args.residues else set(residue_ids)
    selected_indices = [index for index, atom in enumerate(atoms) if atom.residue_sequence in selected]
    chosen = list(atoms) if args.full_context else [atoms[index] for index in selected_indices]

    locs = [atom.position for atom in chosen]
    rads = [atom.radius for atom in chosen]
    group = selected_indices if args.full_context else list(range(len(chosen)))
    print(f"atoms={len(atoms)} context={len(chosen)} selected={len(group)} max_vert={args.max_vert:g}", flush=True)

    network_types = ("aw", "pow", "prm") if args.network == "all" else (args.network,)
    for net_type in network_types:
        for key, value in fast.POW_PRM_METRICS.items():
            fast.POW_PRM_METRICS[key] = 0.0 if isinstance(value, float) else 0
        settings = {
            "surf_res": 0.2, "box_size": 1.25, "max_vert": args.max_vert,
            "build_type": "all", "net_type": net_type, "num_splits": 1,
            "print_metrics": False, "foam_box": None, "sys_dir": ".",
            "ball_type": "normal",
        }
        with contextlib.redirect_stdout(io.StringIO()):
            net = Network(locs=locs, rads=rads, group=group, settings=settings, sort_balls=False)
            started = time.perf_counter()
            net.sort_balls()
            if args.timeout:
                signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError("benchmark timeout")))
                signal.alarm(args.timeout)
            try:
                net.find_verts()
            finally:
                if args.timeout:
                    signal.alarm(0)
            elapsed = time.perf_counter() - started
        vertices = len(net.verts) if net.verts is not None else 0
        timing = dict(getattr(net, "vert_timing", {}))
        timing.update({
            "measured_total": elapsed,
            "avg_seconds_per_vertex": elapsed / vertices if vertices else None,
            "detail_metrics": dict(fast.POW_PRM_METRICS),
        })
        print(json.dumps({"network": net_type, "vertices": vertices, "timing": timing}), flush=True)


if __name__ == "__main__":
    main()
