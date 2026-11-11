"""Run a portable checkpoint manifest, largest first, with bounded concurrency."""
from __future__ import annotations

import argparse
from collections import deque
import csv
from datetime import datetime
import os
from pathlib import Path
import re
import shutil
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
            '-c', 'l', '4', '-e', 'small', '-e', 'dir', str(output), '-v']


def enable_live_console():
    if not sys.stdout.isatty():
        return False
    if os.name == 'nt':
        import ctypes
        import msvcrt
        handle = ctypes.c_void_p(msvcrt.get_osfhandle(sys.stdout.fileno()))
        mode = ctypes.c_ulong()
        kernel = ctypes.windll.kernel32
        return bool(kernel.GetConsoleMode(handle, ctypes.byref(mode))
                    and kernel.SetConsoleMode(handle, mode.value | 0x0004))
    return True


def latest_output(path):
    # Read only the tail: full logs can be large and stay available on disk.
    with path.open('rb') as stream:
        stream.seek(0, os.SEEK_END)
        stream.seek(max(0, stream.tell() - 8192))
        text = stream.read().decode('utf-8', errors='replace')
    lines = re.split(r'[\r\n]+', text)
    for line in reversed(lines):
        line = re.sub(r'\x1b\[[0-?]*[ -/]*[@-~]', '', line).strip()
        if line:
            return ''.join(char for char in line if char.isprintable())
    return 'Starting...'


def progress_message(message):
    match = re.fullmatch(r'\d+:\d+:\d+(?:\.\d+)? - (.+) - ([\d.]+) %', message)
    if not match:
        return message
    context, percentage = match.groups()
    network, separator, process = context.rpartition(' - ')
    stage, _, counts = (process if separator else context).partition(':')
    details = f' ({counts.strip()})' if counts else ''
    return f'{stage} {percentage}%{details}' + (f' | {network}' if separator else '')


class ProgressDisplay:
    def __init__(self, root, workers, total):
        self.root, self.workers, self.total = root, workers, total
        self.live = enable_live_console()
        self.rows = 0
        self.last_refresh = 0.
        self.previous = {}

    def render(self, active, finished, failures, event, *, force=False):
        now = time.monotonic()
        if not force and now - self.last_refresh < 0.5:
            return
        self.last_refresh = now
        lines = [f'VorPy | {finished}/{self.total} completed | {len(active)} running | {failures} failed']
        for _, _, source, _, log_path, started in active:
            label = source.parent.relative_to(self.root).as_posix().replace('/Checkpoints/', '/')
            message = progress_message(latest_output(log_path))
            number = log_path.stem.lstrip('0') or '0'
            row = f'Job {number}/{self.total} | {message} | {int(now - started)}s | {label}'
            lines.append(row)
            if not self.live and self.previous.get(log_path) != message:
                print(row, flush=True)
                self.previous[log_path] = message
        if not self.live:
            return
        lines.extend('Idle' for _ in range(self.workers - len(active)))
        lines.append(event)
        width = max(1, shutil.get_terminal_size((120, 30)).columns - 1)
        redraw = f'\x1b[{self.rows}F' if self.rows else ''
        sys.stdout.write(redraw + '\x1b[J' + '\n'.join(line[:width] for line in lines) + '\n')
        sys.stdout.flush()
        self.rows = len(lines)


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
    env['VORPY_PROGRESS_STREAM'] = '1'
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS'):
        env[key] = '1'
    log_dir = manifest.parent / 'vorpy_layer4_runlogs' / datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    log_dir.mkdir(parents=True)
    print(f'Console logs and summary: {log_dir}', flush=True)
    pending = deque(enumerate(jobs, 1))
    active = []
    failures = finished = 0
    display = ProgressDisplay(root, min(args.workers, len(jobs)), len(jobs))
    event = 'Starting workers...'
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
                    log.write(subprocess.list2cmdline(command) + '\n\nLoading checkpoint...\n')
                    log.flush()
                    try:
                        process = subprocess.Popen(command, cwd=REPO, env=env,
                            stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT)
                    except OSError:
                        log.close()
                        raise
                    active.append((process, log, source, output, log_path, time.monotonic()))
                    event = f'Start {number}/{len(jobs)}: {source.relative_to(root)}'
                    if not display.live:
                        print(event, flush=True)
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
                    event = f'Done {finished}/{len(jobs)} (exit {code}): {source.relative_to(root)}'
                    if not display.live:
                        print(event, flush=True)
                display.render(active, finished, failures, event, force=not active)
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
