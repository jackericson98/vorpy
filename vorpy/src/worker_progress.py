"""Live progress shared by checkpoint and frame worker queues."""
import os
import re
import shutil
import sys
import time


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
    try:
        stream = path.open('rb')
    except FileNotFoundError:
        return 'Starting worker...'
    with stream:
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
    def __init__(self, workers, total, unit='Job'):
        self.workers, self.total, self.unit = workers, total, unit
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
        for number, label, log_path, started in active:
            message = progress_message(latest_output(log_path))
            row = f'{self.unit} {number}/{self.total} | {message} | {int(now - started)}s | {label}'
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

    def announce(self, message):
        if self.live and self.rows:
            sys.stdout.write(f'\x1b[{self.rows}F\x1b[J')
            self.rows = 0
        print(message, flush=True)

