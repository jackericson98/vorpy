import json
import os
import subprocess
import sys
import zipfile
from dataclasses import replace
from pathlib import Path
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


def _grouped_pdb():
    atoms = (
        ('A', 1, 0.0, 0.0, 0.0), ('A', 2, 1.6, 0.0, 0.0),
        ('A', 3, 0.0, 1.6, 0.0), ('A', 4, 0.0, 0.0, 1.6),
        ('B', 5, 2.0, 1.0, 0.0), ('B', 6, 3.6, 1.0, 0.0),
        ('B', 7, 2.0, 2.6, 0.0), ('B', 8, 2.0, 1.0, 1.6),
        ('I', 9, 1.0, 0.5, 0.8), ('I', 10, 2.6, 0.5, 0.8),
        ('I', 11, 1.0, 2.1, 0.8), ('I', 12, 1.0, 0.5, 2.4),
    )
    return ''.join(
        f'ATOM  {serial:5d}  CA  ALA {chain}{residue:4d}'
        f'    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n'
        for serial, (chain, residue, x, y, z) in enumerate(atoms, 1)
    ) + 'END\n'


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

    network.sys.files['dir'] = str(tmp_path / 'geometry-checkpoint')
    checkpoint_command = Command(sys=network.sys)
    checkpoint = checkpoint_command._persist_geometry_checkpoint(
        network, scope='full_system_source'
    )
    restored_checkpoint = load_network(checkpoint)
    assert restored_checkpoint.archive_provenance['geometry_status'] == (
        'GEOMETRY_COMPLETE'
    )
    assert restored_checkpoint.archive_provenance['analysis_status'] == (
        'NOT_CALCULATED'
    )
    assert len(restored_checkpoint.verts) == len(network.verts)
    assert len(restored_checkpoint.edges) == len(network.edges)
    assert len(restored_checkpoint.surfs) == len(network.surfs)

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


@pytest.mark.timeout(120)
def test_real_cli_save_analyze_export_and_optional_failure(tmp_path):
    """Own the minimal CLI-to-filesystem integration contract."""
    root = Path(__file__).resolve().parents[2]
    source = tmp_path / 'grouped.pdb'
    source.write_text(_grouped_pdb(), encoding='ascii')
    output = tmp_path / 'pow-output'
    archive = tmp_path / 'pow.vpy'
    command = [
        sys.executable, '-m', 'vorpy', '--fail-if-exists', str(source),
        '-g', 'c', 'A', 'and', 'c', 'B', '-g', 'c', 'I', '-i', '-c', 's',
        '-s', 'nt', 'pow', 'and', 'sr', '0.2', 'and', 'mv', '30',
        '--boundary', 'shell', '--solvent-mode', 'MOLECULE_ONLY',
        '--profile', 'full', '-e', 'all', '-e', 'dir', str(output),
        '--save-network', str(archive),
    ]
    env = {**os.environ, 'PYTHONPATH': str(root)}
    solved = subprocess.run(
        command, cwd=root, capture_output=True, text=True, timeout=120, env=env,
    )
    assert solved.returncode == 0, solved.stdout + '\n' + solved.stderr
    assert solved.stdout.count('network built -') == 1, solved.stdout
    assert 'Interface Result status:' not in solved.stdout
    assert archive.is_file()

    files = [path for path in output.rglob('*') if path.is_file()]
    assert any(path.suffix.lower() == '.off' for path in files)
    assert any('edge' in path.stem.lower() for path in files)
    assert any('vert' in path.stem.lower() for path in files)
    off = next(path for path in files if path.suffix.lower() == '.off')
    assert off.read_text(encoding='utf-8').splitlines()[0] == 'OFF'
    assert any(path.suffix.lower() == '.csv' for path in files)
    with zipfile.ZipFile(archive) as bundle:
        members = set(bundle.namelist())
        assert 'checkpoint/manifest.json' in members
        assert any(name.startswith('checkpoint/payloads/') for name in members)

    source.unlink()
    reload_output = tmp_path / 'replayed-output'
    reload_code = (
        'import sys; from vorpy.src.command.vpy_cmnd import Command; '
        'sys.argv=["vorpy", sys.argv[1], "-e", "dir", sys.argv[2], "-e", "medium"]; '
        'Command().run()'
    )
    reloaded = subprocess.run(
        [sys.executable, '-c', reload_code, str(archive), str(reload_output)],
        cwd=root, capture_output=True, text=True, timeout=30, env=env,
    )
    assert reloaded.returncode == 0, reloaded.stdout + '\n' + reloaded.stderr
    assert (reload_output / source.stem / 'system' / 'system_info.txt').is_file()
    assert not (tmp_path / 'pow_exports').exists()

    shim = tmp_path / 'optional-failure-shim'
    shim.mkdir()
    source.write_text(_grouped_pdb(), encoding='ascii')
    (shim / 'sitecustomize.py').write_text(
        'from vorpy.src.command.vpy_cmnd import Command\n'
        'Command.run_aw_alpha_experimental = lambda self: (_ for _ in ()).throw(\n'
        "    SystemExit('synthetic GUDHI unavailable'))\n",
        encoding='ascii',
    )
    failure_output = tmp_path / 'optional-failure-output'
    failure_archive = tmp_path / 'optional-failure.vpy'
    failure_command = [
        sys.executable, '-m', 'vorpy', '--fail-if-exists', str(source),
        '-c', 's', '-s', 'nt', 'aw', 'and', 'sr', '0.2', 'and', 'mv', '30',
        '--boundary', 'shell', '--solvent-mode', 'MOLECULE_ONLY',
        '--profile', 'full', '--aw-alpha-experimental', '0', '-e', 'all',
        '-e', 'dir', str(failure_output), '--save-network', str(failure_archive),
    ]
    failure_env = {**os.environ, 'PYTHONPATH': f'{shim}{os.pathsep}{root}'}
    failed_analysis = subprocess.run(
        failure_command, cwd=root, capture_output=True, text=True,
        timeout=120, env=failure_env,
    )
    assert failed_analysis.returncode == 0, (
        failed_analysis.stdout + '\n' + failed_analysis.stderr
    )
    assert 'Optional analysis AW alpha analysis failed' in failed_analysis.stderr
    assert failure_archive.is_file()
    assert any(path.suffix.lower() == '.off' for path in failure_output.rglob('*'))
    restored_failure = Command()._load_canonical_checkpoint(failure_archive)
    failures = restored_failure.archive_provenance['optional_analysis_failures']
    assert failures[0]['status'] == 'UNRESOLVED'
    assert failures[0]['message'] == 'synthetic GUDHI unavailable'
