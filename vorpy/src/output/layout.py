"""Small, deterministic helpers for the on-disk output layout."""

import re
from pathlib import Path


def slug(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("._") or "network"


def interface_folder_name(group1, group2):
    first = slug(getattr(group1, "name", group1) or "group_1")
    second = slug(getattr(group2, "name", group2) or "surrounding")
    return f"interface_{first}_{second}"


def unique_output_directory(root, name, occupied=()):
    """Return a non-colliding child path without creating it."""
    root = Path(root)
    occupied = {Path(path).resolve() for path in occupied}
    candidate = root / slug(name)
    suffix = 2
    while candidate.resolve() in occupied or candidate.exists():
        candidate = root / f"{slug(name)}_{suffix}"
        suffix += 1
    return candidate
