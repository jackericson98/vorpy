import json

import pytest

from vorpy.workbench.domain import Atom
from vorpy.workbench.project import (
    AtomKey,
    GroupDefinition,
    InterfaceDefinition,
    Project,
    StructureSource,
    load_project,
    save_project,
)


def test_project_round_trip_uses_relative_structure_path(tmp_path):
    structure_path = tmp_path / "inputs" / "tiny.pdb"
    structure_path.parent.mkdir()
    structure_path.write_text("END\n", encoding="utf-8")
    atom = Atom(0, 1, "CA", "C", (0.0, 0.0, 0.0), "GLY", "7", "A")
    project = Project(
        name="Example",
        structure=StructureSource.from_path(structure_path),
        groups=[GroupDefinition.create("Backbone", [atom])],
        interfaces=[InterfaceDefinition("iface-id", "Contact", "Backbone", "Backbone")],
    )
    destination = tmp_path / "example.vpyworkbench.json"

    save_project(project, destination)
    payload = json.loads(destination.read_text(encoding="utf-8"))
    restored = load_project(destination)

    assert payload["schema_version"] == 1
    assert payload["project"]["structure"]["source_path"] == "inputs/tiny.pdb"
    assert restored.name == "Example"
    assert restored.structure.source_path == structure_path.resolve()
    assert restored.groups[0].atom_keys == (AtomKey.from_atom(atom),)
    assert restored.interfaces[0].name == "Contact"


def test_project_loader_rejects_unknown_schema(tmp_path):
    source = tmp_path / "future.vpyworkbench.json"
    source.write_text('{"schema_version": 99, "project": {}}', encoding="utf-8")

    with pytest.raises(ValueError, match="Unsupported project schema"):
        load_project(source)
