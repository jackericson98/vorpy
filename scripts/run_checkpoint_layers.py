"""Run a portable checkpoint manifest, largest first, with bounded concurrency."""
from __future__ import annotations

import argparse
from collections import deque
import csv
from datetime import datetime
import os
from pathlib import Path
import subprocess
import sys
import time


REPO = Path(__file__).resolve().parents[1]


def load_jobs(manifest, root):
    with manifest.open(newline='', encoding='utf-8-sig') as stream:
        rows = list(csv.DictReader(stream))
    jobs = []
    outputs = set()
    for row in rows:
        source, output = root / row['checkpoint'], root / row['output']
        if not source.is_file() or source.suffix.lower() != '.vpy':
            raise ValueError(f'Checkpoint not found: {source}')
        if output.resolve() in outputs:
            raise ValueError(f'Duplicate output directory: {output}')
        outputs.add(output.resolve())
        jobs.append((source, output, source.stat().st_size))
    if not jobs:
        raise ValueError('The manifest contains no checkpoints.')
    # shortcut: archive bytes estimate runtime; use measured durations if scheduling tails stay large.
    return sorted(jobs, key=lambda job: (-job[2], str(job[0])))


def build_command(source, output):
    return [sys.executable, '-u', str(REPO / 'vorpy'), str(source),
            '-c', 'l', '4', '-e', 'small', '-e', 'dir', str(output)]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path, help='CSV with checkpoint and output paths relative to its folder.')
    parser.add_argument('--root', type=Path, help='Data root override; defaults to the manifest folder.')
    parser.add_argument('--workers', type=int, default=20, help='Concurrent single-thread jobs (default: 20).')
    parser.add_argument('--dry-run', action='store_true', help='Validate and print commands without running them.')
    args = parser.parse_args(argv)
    if args.workers < 1:
        parser.error('--workers must be positive')
    manifest = args.manifest.resolve()
    root = (args.root or manifest.parent).resolve()
    try:
        jobs = load_jobs(manifest, root)
    except (OSError, ValueError, KeyError) as error:
        parser.error(str(error))
    print(f'{len(jobs)} checkpoints; up to {args.workers} concurrent jobs, largest first.', flush=True)
    if args.dry_run:
        for source, output, _ in jobs:
            print(subprocess.list2cmdline(build_command(source, output)))
        return 0

    env = os.environ.copy()
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS'):
        env[key] = '1'
    log_dir = manifest.parent / 'vorpy_layer4_runlogs' / datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    log_dir.mkdir(parents=True)
    print(f'Console logs and summary: {log_dir}', flush=True)
    pending = deque(enumerate(jobs, 1))
    active = []
    failures = finished = 0
    with (log_dir / 'summary.csv').open('w', newline='', encoding='utf-8') as summary:
        writer = csv.writer(summary)
        writer.writerow(['checkpoint', 'output', 'exit_code', 'seconds', 'console_log'])
        summary.flush()
        try:
            while pending or active:
                while pending and len(active) < args.workers:
                    number, (source, output, _) = pending.popleft()
                    log_path = log_dir / f'{number:04d}.log'
                    log = log_path.open('w', encoding='utf-8')
                    command = build_command(source, output)
                    log.write(subprocess.list2cmdline(command) + '\n\n')
                    log.flush()
                    try:
                        process = subprocess.Popen(command, cwd=REPO, env=env,
                            stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT)
                    except OSError:
                        log.close()
                        raise
                    active.append((process, log, source, output, log_path, time.monotonic()))
                    print(f'Start {number}/{len(jobs)}: {source.relative_to(root)}', flush=True)
                for job in active[:]:
                    process, log, source, output, log_path, started = job
                    code = process.poll()
                    if code is None:
                        continue
                    log.close()
                    active.remove(job)
                    finished += 1
                    failures += code != 0
                    writer.writerow([str(source), str(output), code,
                                     round(time.monotonic() - started, 2), str(log_path)])
                    summary.flush()
                    print(f'Done {finished}/{len(jobs)} (exit {code}): {source.relative_to(root)}', flush=True)
                if active:
                    time.sleep(0.2)
        finally:
            for process, log, *_ in active:
                if process.poll() is None:
                    process.terminate()
                process.wait()
                log.close()
    print(f'Completed {finished} jobs; {failures} failed. Summary: {log_dir / "summary.csv"}', flush=True)
    return int(failures > 0)


if __name__ == '__main__':
    raise SystemExit(main())
