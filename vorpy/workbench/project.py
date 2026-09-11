"""Versioned, Qt-free Workbench project state and JSON persistence."""

from __future__ import annotations

import hashlib
import json
import numpy as np
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

from vorpy.workbench.domain import AnalysisResult, Atom, Bond, GeometryLayer

PROJECT_SCHEMA_VERSION = 1
PROJECT_SUFFIX = ".vpyworkbench.json"


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass(frozen=True)
class AtomKey:
    """Stable atom identity using the fields available from the PDB loader."""

    chain: str
    residue_sequence: str
    residue_name: str
    atom_name: str

    @classmethod
    def from_atom(cls, atom: Atom) -> AtomKey:
        return cls(
            chain=atom.chain,
            residue_sequence=atom.residue_sequence,
            residue_name=atom.residue_name,
            atom_name=atom.name,
        )


@dataclass
class StructureSource:
    id: str
    name: str
    source_path: Path
    fingerprint: str = ""

    @classmethod
    def from_path(cls, path: Path, name: str | None = None) -> StructureSource:
        resolved = path.expanduser().resolve()
        return cls(
            id=str(uuid4()),
            name=name or resolved.stem,
            source_path=resolved,
            fingerprint=fingerprint_path(resolved),
        )


@dataclass
class GroupDefinition:
    id: str
    name: str
    atom_keys: tuple[AtomKey, ...]
    color: str = "#55d6c2"
    description: str = ""

    @classmethod
    def create(cls, name: str, atoms: list[Atom]) -> GroupDefinition:
        return cls(
            id=str(uuid4()),
            name=name,
            atom_keys=tuple(AtomKey.from_atom(atom) for atom in atoms),
        )


@dataclass
class InterfaceDefinition:
    id: str
    name: str
    group_a: str
    group_b: str


@dataclass
class Project:
    id: str = field(default_factory=lambda: str(uuid4()))
    name: str = "Untitled Project"
    created_at: str = field(default_factory=_now)
    modified_at: str = field(default_factory=_now)
    structure: StructureSource | None = None
    groups: list[GroupDefinition] = field(default_factory=list)
    interfaces: list[InterfaceDefinition] = field(default_factory=list)
    result_state: dict | None = None
    view_state: dict = field(default_factory=dict)
    backend_settings: dict = field(default_factory=dict)


def fingerprint_path(path: Path) -> str:
    """Return a content fingerprint for a file or a stable path marker for a directory."""
    if not path.is_file():
        return f"directory:{path.name}"
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return f"sha256:{digest.hexdigest()}"


def save_project(project: Project, destination: Path) -> None:
    """Atomically write a project while keeping large scientific files external."""
    destination = destination.expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    project.modified_at = _now()
    payload = {
        "schema_version": PROJECT_SCHEMA_VERSION,
        "project": {
            "id": project.id,
            "name": project.name,
            "created_at": project.created_at,
            "modified_at": project.modified_at,
            "structure": _structure_to_json(project.structure, destination.parent),
            "groups": [
                {
                    "id": group.id,
                    "name": group.name,
                    "color": group.color,
                    "description": group.description,
                    "atom_keys": [asdict(key) for key in group.atom_keys],
                }
                for group in project.groups
            ],
            "interfaces": [
                {
                    "id": interface.id,
                    "name": interface.name,
                    "group_a": interface.group_a,
                    "group_b": interface.group_b,
                }
                for interface in project.interfaces
            ],
            "result_state": project.result_state,
            "view_state": project.view_state,
            "backend_settings": project.backend_settings,
        },
    }
    temporary = destination.with_name(f".{destination.name}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(destination)


def load_project(source: Path) -> Project:
    source = source.expanduser().resolve()
    payload = json.loads(source.read_text(encoding="utf-8"))
    version = payload.get("schema_version")
    if version != PROJECT_SCHEMA_VERSION:
        raise ValueError(
            f"Unsupported project schema {version!r}; expected {PROJECT_SCHEMA_VERSION}."
        )
    data = payload.get("project")
    if not isinstance(data, dict):
        raise TypeError("Project file does not contain a project object.")
    structure_data = data.get("structure")
    structure = None
    if structure_data is not None:
        path = Path(structure_data["source_path"])
        if not path.is_absolute():
            path = (source.parent / path).resolve()
        structure = StructureSource(
            id=structure_data["id"],
            name=structure_data["name"],
            source_path=path,
            fingerprint=structure_data.get("fingerprint", ""),
        )
    groups = [
        GroupDefinition(
            id=group["id"],
            name=group["name"],
            color=group.get("color", "#55d6c2"),
            description=group.get("description", ""),
            atom_keys=tuple(AtomKey(**key) for key in group.get("atom_keys", [])),
        )
        for group in data.get("groups", [])
    ]
    interfaces = [
        InterfaceDefinition(
            id=interface["id"],
            name=interface["name"],
            group_a=interface["group_a"],
            group_b=interface["group_b"],
        )
        for interface in data.get("interfaces", [])
    ]
    return Project(
        id=data["id"],
        name=data["name"],
        created_at=data["created_at"],
        modified_at=data["modified_at"],
        structure=structure,
        groups=groups,
        interfaces=interfaces,
        result_state=data.get("result_state"),
        view_state=data.get("view_state", {}),
        backend_settings=data.get("backend_settings", {}),
    )


def _structure_to_json(
    structure: StructureSource | None, project_directory: Path
) -> dict | None:
    if structure is None:
        return None
    source_path = structure.source_path
    try:
        stored_path = source_path.relative_to(project_directory).as_posix()
    except ValueError:
        stored_path = source_path.as_posix()
    return {
        "id": structure.id,
        "name": structure.name,
        "source_path": stored_path,
        "fingerprint": structure.fingerprint,
    }


def result_to_json(result: AnalysisResult) -> dict:
    """Serialize a solved result, including all mesh and scalar arrays."""
    return {
        "source": str(result.source) if result.source is not None else None,
        "name": result.name,
        "atoms": [asdict(atom) for atom in result.atoms],
        "bonds": [asdict(bond) for bond in result.bonds],
        "layers": [
            {
                "name": layer.name,
                "kind": layer.kind,
                "points": layer.points.tolist(),
                "lines": layer.lines.tolist() if layer.lines is not None else None,
                "faces": layer.faces.tolist() if layer.faces is not None else None,
                "source_path": str(layer.source_path) if layer.source_path else None,
                "color": layer.color,
                "opacity": layer.opacity,
                "visible": layer.visible,
                "cell_scalars": {
                    key: values.tolist() for key, values in layer.cell_scalars.items()
                },
                "color_scheme": layer.color_scheme,
                "color_map": layer.color_map,
                "scale_mode": layer.scale_mode,
                "interpretation": layer.interpretation,
            }
            for layer in result.layers
        ],
        "complete_cells": result.complete_cells,
        "surface_count": result.surface_count,
        "elapsed_seconds": result.elapsed_seconds,
        "info_sections": result.info_sections,
    }


def result_from_json(data: dict) -> AnalysisResult:
    layers = []
    for item in data.get("layers", []):
        source_path = item.get("source_path")
        layers.append(
            GeometryLayer(
                name=item["name"],
                kind=item["kind"],
                points=np.asarray(item.get("points", []), dtype=float),
                lines=(
                    np.asarray(item["lines"], dtype=np.int64)
                    if item.get("lines") is not None
                    else None
                ),
                faces=(
                    np.asarray(item["faces"], dtype=np.int64)
                    if item.get("faces") is not None
                    else None
                ),
                source_path=Path(source_path) if source_path else None,
                color=item.get("color", "#55a9d9"),
                opacity=float(item.get("opacity", 1.0)),
                visible=bool(item.get("visible", True)),
                cell_scalars={
                    key: np.asarray(values, dtype=float)
                    for key, values in item.get("cell_scalars", {}).items()
                },
                color_scheme=item.get("color_scheme", "solid"),
                color_map=item.get("color_map", "coolwarm"),
                scale_mode=item.get("scale_mode", "signed_log"),
                interpretation=item.get("interpretation", "magnitude"),
            )
        )
    atoms = [
        Atom(
            index=int(item["index"]),
            serial=int(item["serial"]),
            name=item["name"],
            element=item["element"],
            position=tuple(item["position"]),
            residue_name=item.get("residue_name", ""),
            residue_sequence=item.get("residue_sequence", ""),
            chain=item.get("chain", ""),
            radius=float(item.get("radius", 0.35)),
        )
        for item in data.get("atoms", [])
    ]
    return AnalysisResult(
        source=Path(data["source"]) if data.get("source") else None,
        name=data["name"],
        atoms=atoms,
        bonds=[Bond(int(item["atom_a"]), int(item["atom_b"])) for item in data.get("bonds", [])],
        layers=layers,
        complete_cells=int(data.get("complete_cells", 0)),
        surface_count=int(data.get("surface_count", 0)),
        elapsed_seconds=float(data.get("elapsed_seconds", 0.0)),
        info_sections={
            str(section): [(str(key), str(value)) for key, value in values]
            for section, values in data.get("info_sections", {}).items()
        },
    )