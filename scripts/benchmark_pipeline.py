#!/usr/bin/env python3
"""Bounded CLI performance audit; no persistent production instrumentation.

Each trial uses a fresh process and output directory. The uninstrumented control
warms on-disk Numba/OS caches; measured trials still include per-process JIT cache
loading. Pipeline timing excludes imports, instrumentation setup and validation.
Use --warm-kernels to run an additional uninstrumented CLI warm-up in each
measured worker before its timed build. GNU time peak RSS covers the complete
worker, including imports, validation, and the optional warm-up.
"""
import argparse
from contextlib import nullcontext
import csv
import gc
import hashlib
import importlib
import importlib.metadata
import io
import json
import math
import os
from pathlib import Path
import platform
import re
import signal
import statistics
import subprocess
import sys

from performance_timing import Recorder, instrument

REPO = Path(__file__).resolve().parents[1]
CATEGORIES = ('input_setup', 'spatial_indexing', 'network_construction',
              'network_analysis', 'geometry_preparation', 'file_export',
              'orchestration')


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def canonical(value):
    """Keep exact numeric values and sequence order; canonicalize unordered sets."""
    import numpy as np
    if isinstance(value, np.ndarray):
        return canonical(value.tolist())
    if isinstance(value, np.generic):
        return canonical(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return {'nonfinite': str(value)}
    if isinstance(value, dict):
        return {str(k): canonical(v) for k, v in value.items()}
    if isinstance(value, (set, frozenset)):
        return {'set': sorted((canonical(v) for v in value), key=lambda v: json.dumps(v, sort_keys=True))}
    if isinstance(value, (tuple, list)):
        return [canonical(v) for v in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    raise TypeError(f'Unsupported snapshot value: {type(value).__name__}')


def snapshot(system):
    groups = []
    for group in system.groups:
        net = group.net
        frames = {}
        for name in ('balls', 'verts', 'edges', 'surfs'):
            frame = getattr(net, name)
            if frame is None:
                raise ValueError(f'Incomplete build: missing {name}')
            frames[name] = {
                'index': canonical(frame.index.tolist()),
                'columns': {str(col): digest(canonical(frame[col].tolist())) for col in frame.columns},
            }
        groups.append({
            'name': group.name, 'membership': canonical(group.ball_ndxs),
            'frames': frames,
            'selected_balls': len(group.ball_ndxs),
            'counts': {**{name: len(getattr(net, name)) for name in frames},
                       'complete_cells': int(net.balls['complete'].sum()),
                       'incomplete_cells': int((~net.balls['complete'].astype(bool)).sum())},
        })
    if not groups:
        raise ValueError('CLI produced no groups')
    return groups


def output_manifest(root):
    """Normalize only trial directory, known log timestamps/times, info timing rows."""
    manifest = {}
    for path in sorted(root.rglob('*')):
        if not path.is_file():
            continue
        content = path.read_text().replace(str(root.parent), '<TRIAL>')
        if path.name.endswith('_logs.csv'):
            rows = list(csv.reader(io.StringIO(content)))
            for index, row in enumerate(rows[:-1]):
                if 'Completion Date' in row and 'Vertex Time' in row:
                    for column, label in enumerate(row):
                        if label in {'Completion Date', 'Total Time', 'Vertex Time', 'Connect Time',
                                     'Surface Building Time', 'Analysis time'}:
                            rows[index + 1][column] = '<VOLATILE>'
            content = json.dumps(rows)
        elif path.name.endswith('_info.txt') or path.name == 'info.txt':
            content = re.sub(r'(?m)^(?:Vertex Time|Connection Time|Surface Building Time|Analysis Time|Total Time):[^\n]*$',
                             '<TIMING>', content)
        manifest[str(path.relative_to(root))] = hashlib.sha256(content.encode()).hexdigest()
    return manifest


def worker(args):
    warmup = None
    if args.warm_kernels and not args.control:
        warmup_directory = args.output / 'warmup'
        warmup_directory.mkdir()
        warmup_args = argparse.Namespace(**{**vars(args), 'output': warmup_directory,
                                           'control': True, 'warm_kernels': False})
        worker(warmup_args)
        warmup = json.loads((warmup_directory / 'result.json').read_text())
        gc.collect()
    # Import before pipeline timing, exactly as the normal CLI does.
    from vorpy.src.command.vpy_cmnd import Command
    from vorpy.src.system.system import System
    from vorpy.src.group.group import Group
    from vorpy.src.network.network import Network
    from vorpy.src.network import fast
    vertex_module = importlib.import_module('vorpy.src.network.find_net_verts')
    trial = args.output.resolve()
    artifacts = trial / 'artifacts'
    original_set_files = System.set_files

    def isolated_files(self, *positional, **kwargs):
        original_set_files(self, *positional, **kwargs)
        # Redirect the CLI's ordinary output allocation before input loading.
        self.files.update(vpy_dir=str(artifacts), root_dir=str(artifacts), dir=str(artifacts))

    hooks = [
        (System, '__init__', 'system_input', 'input_setup'),
        *[(Command, name, name, 'input_setup') for name in
          ('parse_commands', 'load_files', 'apply_settings', 'create_groups')],
        (Network, 'sort_balls', 'sort_balls', 'spatial_indexing'),
        (System, 'apply_spatial_index', 'reuse_spatial_index', 'spatial_indexing'),
        (Network, 'build', 'network_build', 'network_construction'),
        (Network, 'find_verts', 'vertex_search', 'network_construction'),
        (Network, 'connect', 'topology', 'network_construction'),
        (Network, 'build_edges', 'edge_geometry', 'network_construction'),
        (Network, 'build_surfaces', 'surface_geometry', 'network_construction'),
        (Network, 'analyze', 'analysis', 'network_analysis'),
        (Command, 'run_exports', 'export_orchestration', 'file_export'),
        (System, 'exports', 'system_exports', 'file_export'),
        (Group, 'exports', 'group_exports', 'file_export'),
        (Group, 'get_info', 'group_info', 'geometry_preparation'),
        (Group, 'get_layers', 'group_layers', 'geometry_preparation'),
        (vertex_module, 'write_verts', 'vertex_file', 'file_export'),
    ]
    recorder = Recorder()
    cli = ['vorpy', str(args.pdb.resolve()), '-s', 'nt', args.network,
           '-s', 'mv', '5', '-s', 'sr', '0.2', '-s', 'bs', '1.25', '-e', 'small']
    old_argv, old_cwd = sys.argv, Path.cwd()
    System.set_files = isolated_files
    sys.argv = cli
    for key, value in fast.POW_PRM_METRICS.items():
        fast.POW_PRM_METRICS[key] = 0.0 if isinstance(value, float) else 0
    try:
        os.chdir(trial)
        with instrument(recorder, hooks) if not args.control else nullcontext():
            with recorder.span('pipeline_total', 'orchestration'):
                command = Command()
                command.run()
        result = {
            'schema_version': 1, 'control': args.control, 'cli': cli,
            'same_process_warmup': warmup is not None,
            'timing': recorder.report(), 'snapshot': snapshot(command.sys),
            'output_manifest': output_manifest(artifacts),
            'settings': canonical(command.sys.groups[0].net.settings),
            'legacy_details': {
                name: canonical(getattr(command.sys.groups[0].net, name, {}))
                for name in ('vert_timing', 'surface_timing', 'build_surf_timing', 'analysis_timing')
            },
            'vertex_counters': canonical(fast.POW_PRM_METRICS),
            'environment': {
                'python': sys.version, 'platform': platform.platform(),
                'logical_cpus': os.cpu_count(),
                'cpu_model': next((line.split(':', 1)[1].strip() for line in
                                   Path('/proc/cpuinfo').read_text().splitlines()
                                   if line.startswith('model name')), platform.processor()),
                'hash_seed': os.environ.get('PYTHONHASHSEED', 'random (Python default)'),
                'packages': {name: importlib.metadata.version(name) for name in
                             ('numpy', 'scipy', 'pandas', 'numba', 'shapely')},
            },
        }
        if warmup is not None:
            result['equivalent_to_same_process_warmup'] = {
                'network': result['snapshot'] == warmup['snapshot'],
                'normalized_files': result['output_manifest'] == warmup['output_manifest'],
            }
            if not all(result['equivalent_to_same_process_warmup'].values()):
                raise RuntimeError('Same-process warm-up correctness mismatch')
        (trial / 'result.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    finally:
        System.set_files = original_set_files
        sys.argv = old_argv
        os.chdir(old_cwd)


def run_trial(args, directory, control=False):
    directory.mkdir()
    command = ['/usr/bin/time', '-v', '-o', str(directory / 'resource.txt'),
               sys.executable, str(Path(__file__).resolve()), str(args.pdb.resolve()),
               '--worker', '--output', str(directory), '--network', args.network]
    if control:
        command.append('--control')
    if args.warm_kernels:
        command.append('--warm-kernels')
    env = dict(os.environ, PYTHONPATH=str(REPO), LC_ALL='C')
    with (directory / 'console.log').open('w') as log:
        process = subprocess.Popen(command, cwd=directory, env=env, stdout=log,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            status = process.wait(timeout=args.timeout)
        except BaseException as error:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
            if isinstance(error, subprocess.TimeoutExpired):
                raise RuntimeError(f'Trial exceeded {args.timeout}s; see {directory / "console.log"}') from error
            raise
    if status:
        raise RuntimeError(f'Trial failed ({status}); see {directory / "console.log"}')
    result = json.loads((directory / 'result.json').read_text())
    resource = (directory / 'resource.txt').read_text()
    result['peak_rss_kib'] = int(re.search(r'Maximum resident set size \(kbytes\): (\d+)', resource)[1])
    result['worker_user_seconds'] = float(re.search(r'User time \(seconds\): ([0-9.]+)', resource)[1])
    result['worker_system_seconds'] = float(re.search(r'System time \(seconds\): ([0-9.]+)', resource)[1])
    result['worker_elapsed'] = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.+)', resource)[1]
    (directory / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


def historical_baseline_comparison(fixture):
    path = REPO / 'docs/performance/network_build_baselines.csv'
    with path.open(newline='') as stream:
        rows = list(csv.reader(stream))
    header = rows[0]
    invalid = [{'line': index, 'fields': len(row)}
               for index, row in enumerate(rows[1:], 2) if len(row) != len(header)]
    matching = [dict(zip(header, row)) for row in rows[1:]
                if len(row) == len(header) and row[0].casefold() == fixture.stem.casefold()]
    return {
        'source': str(path.relative_to(REPO)),
        'source_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'invalid_rows': invalid, 'matching_fixture_rows': matching,
        'speedup_comparison': None,
        'reason': ('Existing measurements do not establish an equivalent CLI pipeline baseline; '
                   'fixture, selected group, context, loader, stage scope and export preset must all match.'),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pdb', type=Path)
    parser.add_argument('--output', type=Path, required=True, help='new results directory')
    parser.add_argument('--network', choices=('pow', 'prm'), default='pow')
    parser.add_argument('--repetitions', type=int, default=3)
    parser.add_argument('--timeout', type=int, default=120, help='whole worker seconds, including imports/validation')
    parser.add_argument('--warm-kernels', action='store_true',
                        help='also warm Numba in each measured worker with an uninstrumented CLI build')
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--control', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.timeout <= 0 or args.repetitions < 1:
        parser.error('timeout and repetitions must be positive')
    if args.worker:
        worker(args)
        return
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    print('Running uninstrumented control / disk-cache warm-up', flush=True)
    control = run_trial(args, args.output / 'control', control=True)
    results = []
    for index in range(args.repetitions):
        print(f'Running measured repetition {index + 1}/{args.repetitions}', flush=True)
        result = run_trial(args, args.output / f'run-{index + 1}')
        result['equivalent_to_control'] = {
            'network': result['snapshot'] == control['snapshot'],
            'normalized_files': result['output_manifest'] == control['output_manifest'],
        }
        (args.output / f'run-{index + 1}' / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
        results.append(result)
        if not all(result['equivalent_to_control'].values()):
            (args.output / 'mismatch.json').write_text(json.dumps(result['equivalent_to_control'], indent=2))
            raise RuntimeError(f'Correctness mismatch; inspect control/run-{index + 1} results')
    rows = []
    for index, result in enumerate(results, 1):
        total = result['timing']['events'][0]
        row = {'run': index, 'wall_seconds': total['wall_seconds'],
               'cpu_seconds': total['cpu_seconds'], 'peak_rss_kib': result['peak_rss_kib'],
               **result['snapshot'][0]['counts']}
        for category in CATEGORIES:
            times = result['timing']['exclusive_categories'].get(category, {})
            for clock in ('wall', 'cpu'):
                row[f'{category}_{clock}_seconds'] = times.get(f'{clock}_seconds', 0.0)
        rows.append(row)
    with (args.output / 'summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        'schema_version': 1, 'fixture': str(args.pdb.resolve()), 'network': args.network,
        'input_sha256': hashlib.sha256(args.pdb.read_bytes()).hexdigest(),
        'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip(),
        'harness_sha256': {name: hashlib.sha256((REPO / 'scripts' / name).read_bytes()).hexdigest()
                           for name in ('benchmark_pipeline.py', 'performance_timing.py')},
        'selected_balls': results[0]['snapshot'][0]['selected_balls'],
        'counts': results[0]['snapshot'][0]['counts'],
        'environment': results[0]['environment'],
        'repetitions': args.repetitions, 'export_preset': 'tiny (CLI: -e small)',
        'warmup': ('one uninstrumented control plus one uninstrumented same-process warm-up per measured worker'
                   if args.warm_kernels else
                   'one uninstrumented fresh-process control; disk caches warm, process state fresh'),
        'all_equivalent': True,
        'historical_baseline_comparison': historical_baseline_comparison(args.pdb),
        'wall_median_seconds': statistics.median(row['wall_seconds'] for row in rows),
        'wall_min_seconds': min(row['wall_seconds'] for row in rows),
        'wall_max_seconds': max(row['wall_seconds'] for row in rows),
        'runs': rows,
    }
    (args.output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
