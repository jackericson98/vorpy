"""Streaming input with bounded progress notifications, independent of any GUI."""
from pathlib import Path


def iter_input_lines(path, progress=None, *, errors="strict"):
    path = Path(path)
    size = max(path.stat().st_size, 1)
    consumed = 0
    label = f"Reading {path.name}"
    if progress:
        progress(label, 0)
    with path.open("rb") as handle:
        for number, raw in enumerate(handle, 1):
            consumed += len(raw)
            yield raw.decode("utf-8", errors=errors)
            if progress and number % 4096 == 0:
                progress(f"{label} · {number:,} lines", min(99, consumed * 100 // size))
    if progress:
        progress(label, 100)
