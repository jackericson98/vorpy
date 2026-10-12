import json
import sys
from dataclasses import replace
from types import SimpleNamespace

import pytest

from vorpy.src.command.vpy_cmnd import Command
from vorpy.src.group.shared_view import attach_group_view_results
from vorpy.src.io import load_network, save_network
from vorpy.src.network import Network
from vorpy.src.results import (
    GeometryLineageOrigin,
    ResultsDependency,
    SharedNetworkIdentity,
    Status,
    validate_results_dependency,
)
from vorpy.workbench.services.vorpy_backend import (
    VorPyBackend,
    VorPySolveSettings,
)


TINY_PDB = """\
ATOM      1  C1  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  C2  ALA A   2       2.000   0.000   0.000  1.00  0.00           C
ATOM      3  C3  ALA A   3       0.000   2.000   0.000  1.00  0.00           C
ATOM      4  C4  ALA A   4       0.000   0.000   2.000  1.00  0.00           C
END
"""


@pytest.mark.timeout(90)
def test_nexus_solve_results_archive_reload_cache_export(tmp_path, monkeypatch):
    """Exercise the shared checkpoint on a bounded synthetic shell solve."""
    monkeypatch.chdir(tmp_path)
    source = tmp_path / "tiny.pdb"
    source.write_text(TINY_PDB, encoding="utf-8")

    result = VorPyBackend(
        VorPySolveSettings(
            network_type="aw",
            boundary_mode="shell",
            max_vertices=5,
            typed_group_selections={"tiny-group": (0,)},
            source_network_id="producer:tiny",
        )
    ).solve(source, lambda *args: None, lambda: False)

    typed = result.typed_results[0]
    network = result.export_group.net
    context = network.get_solve_context_key()
    assert typed.status is Status.CERTIFIED
    assert typed.identity.source_network_id == "producer:tiny"
    assert typed.identity.source_identity_status == "AUTHORITATIVE"
    assert typed.identity.solve_universe_id == context.solve_universe_id
    assert typed.identity.selection_id.startswith("view:")
    assert typed.provenance.details["shared_network_identity"]["selection_id"] == (
        typed.identity.selection_id
    )
    assert typed.provenance.geometry_lineage is not None
    assert typed.provenance.geometry_lineage.origin is (
        GeometryLineageOrigin.MATERIALIZED_VIEW
    )
    assert typed.provenance.geometry_lineage.source_network_id == "producer:tiny"
    assert typed.provenance.geometry_lineage.source_view_id == typed.identity.selection_id
    assert typed.provenance.geometry_lineage.solve_universe_id == context.solve_universe_id

    identity = SharedNetworkIdentity(
        **typed.provenance.details["shared_network_identity"]
    )
    dependency = ResultsDependency(
        identity,
        source_id=typed.provenance.details["cache_dependencies"][0]["source_id"],
        source_digest=context.digest,
    )
    assert validate_results_dependency(dependency, dependency) is Status.CERTIFIED
    assert validate_results_dependency(
        dependency, replace(dependency, source_revision="changed")
    ) is Status.STALE
    assert validate_results_dependency(dependency, None) is Status.NOT_CALCULATED

    archive = save_network(network, tmp_path / "tiny.vpy")
    restored = load_network(archive)
    assert restored.source_network_id == "producer:tiny"
    assert restored.get_solve_context_key() == context

    restored_result = SimpleNamespace()
    restored_typed = attach_group_view_results(
        restored_result,
        restored,
        (SimpleNamespace(ball_ndxs=(0,), group_id="tiny-group"),),
        system="tiny",
        radius_configuration_id="vorpy-backend-radii-v1",
    )[0]
    assert restored_typed.status is Status.CERTIFIED
    assert restored_typed.identity.source_network_id == typed.identity.source_network_id
    assert restored_typed.identity.solve_universe_id == typed.identity.solve_universe_id
    assert restored_typed.identity.selection_id == typed.identity.selection_id
    assert restored_typed.provenance.details["dependency_digest"] == (
        typed.provenance.details["dependency_digest"]
    )

    unresolved_network = load_network(archive)
    unresolved_network.source_network_id = None
    unresolved = attach_group_view_results(
        SimpleNamespace(),
        unresolved_network,
        (SimpleNamespace(ball_ndxs=(0,), group_id="tiny-group"),),
        system="tiny",
        radius_configuration_id="vorpy-backend-radii-v1",
    )[0]
    assert unresolved.status is Status.UNRESOLVED
    assert unresolved.identity.source_identity_status == "PROVISIONAL"
    assert unresolved.provenance.geometry_lineage is not None
    # Identity is unresolved, while the materialized view remains validated.
    assert unresolved.provenance.geometry_lineage.origin is (
        GeometryLineageOrigin.MATERIALIZED_VIEW
    )
    assert unresolved.provenance.geometry_lineage.source_view_id == (
        unresolved.identity.selection_id
    )
    assert unresolved.provenance.geometry_lineage.source_network_id != "producer:tiny"

    monkeypatch.setattr(
        Network,
        "build",
        lambda *args, **kwargs: pytest.fail("cache-only export rebuilt the network"),
    )
    monkeypatch.setattr(
        sys,
        "argv",
        ["vorpy", "--profile-cache-export", str(archive), "-e", "only", "logs"],
    )
    command = Command()
    command.run()
    profile_path = tmp_path / "tiny_exports" / "cache_export_profile.json"
    payload = json.loads(profile_path.read_text(encoding="utf-8"))
    assert payload["format"] == "vorpy.cache_export_profile"
    assert payload["cache_state"]["interface_count"] == len(command.sys.ifaces)
    assert payload["export_phases"]["operations"]

    assert restored_typed.provenance.geometry_lineage is not None
    assert restored_typed.provenance.geometry_lineage == typed.provenance.geometry_lineage
