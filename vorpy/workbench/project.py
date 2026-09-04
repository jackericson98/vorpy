"""Versioned, Qt-free Workbench project state and JSON persistence."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

from vorpy.workbench.domain import Atom

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
    )


def _structure_to_json(
    structure: StructureSource | None, project_directory: Path
) -> dict | None:
    if structure is None:
        return None
    source_path = structure.source_path
    try:
        stored_path = str(source_path.relative_to(project_directory))
    except ValueError:
        stored_path = str(source_path)
    return {
        "id": structure.id,
        "name": structure.name,
        "source_path": stored_path,
        "fingerprint": structure.fingerprint,
    }
