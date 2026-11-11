"""Small, filename-only labels for the Workbench and CLI input context."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re


@dataclass(frozen=True)
class SourceContext:
    dataset_id: str
    role: str
    evidence: str
    frame_index: int | None = None
    frame_count: int | None = None
    paired_path: Path | None = None

    @property
    def role_label(self) -> str:
        return {
            "original": "Original pre-MD",
            "solvated": "Solvated post-MD",
            "dry": "Dry post-MD",
            "frame": "Frame series",
        }.get(self.role, "Structure")


def describe_source(
    path: str | Path | None,
    *,
    frame_index: int | None = None,
    frame_count: int | None = None,
) -> SourceContext | None:
    """Describe the repository's X_og/X/X_no_sol/X_frames naming convention."""
    if path is None:
        return None
    source = Path(path)
    stem = source.stem
    parent = source.parent.name
    frame_match = re.fullmatch(r"(.+)_f(\d+)", stem, re.IGNORECASE)
    frame_directory = re.fullmatch(r"(.+)_frames", parent, re.IGNORECASE)
    if frame_match and frame_directory and frame_match.group(1).lower() == frame_directory.group(1).lower():
        dataset_id = frame_directory.group(1)
        role = "frame"
        derived_frame = int(frame_match.group(2))
        paired_path = source.parent.parent / f"{dataset_id}.pdb"
        evidence = "filename convention"
    elif stem.lower().endswith("_og"):
        dataset_id = stem[:-3]
        role = "original"
        paired_path = source.with_name(f"{dataset_id}.pdb")
        evidence = "paired filename" if paired_path.is_file() else "filename convention"
        derived_frame = None
    elif stem.lower().endswith("_no_sol"):
        dataset_id = stem[:-7]
        role = "dry"
        paired_path = source.with_name(f"{dataset_id}.pdb")
        evidence = "paired filename" if paired_path.is_file() else "filename convention"
        derived_frame = None
    else:
        dataset_id = stem
        role = "solvated"
        paired_path = source.with_name(f"{stem}_no_sol.pdb")
        evidence = "paired filename" if paired_path.is_file() else "filename convention"
        derived_frame = None
    return SourceContext(
        dataset_id=dataset_id,
        role=role,
        evidence=evidence,
        frame_index=frame_index if frame_index is not None else derived_frame,
        frame_count=frame_count,
        paired_path=paired_path,
    )


def source_display_label(
    path: str | Path | None,
    *,
    frame_index: int | None = None,
    frame_count: int | None = None,
) -> str:
    """Return a compact, explicit label suitable for a selector or status bar."""
    context = describe_source(path, frame_index=frame_index, frame_count=frame_count)
    if context is None:
        return "Unknown input"
    label = f"{context.dataset_id} · {context.role_label} ({context.evidence})"
    if context.frame_index is not None:
        count = context.frame_count if context.frame_count is not None else "?"
        label += f" · frame {context.frame_index}/{count}"
    return label
