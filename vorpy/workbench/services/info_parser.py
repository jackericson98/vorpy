"""Parser for the human-readable group ``info.txt`` output."""

from __future__ import annotations

from pathlib import Path


def parse_group_info(path: Path) -> dict[str, list[tuple[str, str]]]:
    """Return output-info headings and their labelled values."""
    sections: dict[str, list[tuple[str, str]]] = {}
    current = "Overview"
    sections[current] = []
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return sections
    for raw in lines:
        line = raw.strip()
        if not line or set(line) <= {"=", "-"}:
            continue
        if ":" not in line and line.isupper() and len(line) < 80:
            current = line.title()
            sections.setdefault(current, [])
            continue
        if ":" in line:
            key, value = line.split(":", 1)
            key, value = key.strip(), value.strip()
            if key:
                sections.setdefault(current, []).append((key, value))
    return {name: values for name, values in sections.items() if values}
