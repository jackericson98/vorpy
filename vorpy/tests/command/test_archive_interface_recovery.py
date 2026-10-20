"""Bounded recovery tests for an existing solved source checkpoint."""

import sys
from pathlib import Path

import pytest

from vorpy.src.command.vpy_cmnd import Command
from vorpy.src.group.group import Group
from vorpy.src.io import load_network, validate_network_source
from vorpy.src.network.network import Network
from vorpy.src.results.model import GeometryLineageOrigin, Status


CHECKPOINT = (
    Path(__file__).resolve().parents[3]
    / 'output' / 'benchmarks' / '1bna_DB1976' / 'pow_dry'
    / '1bna_DB1976_no_sol'
    / '1bna_DB1976_no_sol_full_system.geometry.vpy'
)


def test_checkpoint_source_validation_is_explicit_and_readable():
    if not CHECKPOINT.is_file():
        pytest.skip('verified 1BNA checkpoint is not present in this checkout')
    network = load_network(CHECKPOINT)
    validation = validate_network_source(
        network,
        network.archive_metadata['original_structure_path'],
        expected_scheme='pow',
        expected_solvent_mode='MOLECULE_ONLY',
    )
    assert validation['generator_count'] == 805
    assert validation['complete_cell_count'] == 525
    assert validation['source_sha256'] == network.archive_metadata['source_sha256']


def test_checkpoint_interface_recovery_does_not_rebuild_and_exports_partial_result(
    tmp_path, monkeypatch
):
    if not CHECKPOINT.is_file():
        pytest.skip('verified 1BNA checkpoint is not present in this checkout')

    def forbidden_build(*args, **kwargs):
        raise AssertionError('checkpoint recovery must not rebuild geometry')

    monkeypatch.setattr(Network, 'build', forbidden_build)
    monkeypatch.setattr(Group, 'build', forbidden_build)
    destination = tmp_path / 'recovered_pow_dry'
    repo = CHECKPOINT.parents[5]
    monkeypatch.setattr(
        sys,
        'argv',
        [
            'vorpy', '--load-network', str(CHECKPOINT),
            '-g', 'c A and c B', '-g', 'c C', '-i',
            '-e', 'dir', str(destination), '-e', 'interface',
        ],
    )

    command = Command()
    try:
        command.run()
    finally:
        # System output setup may change cwd; restore it before pytest cleans
        # its temporary directory on Windows.
        monkeypatch.chdir(repo)

    assert command.sys.archive_recovery_validation['generator_count'] == 805
    assert command.sys.archive_recovery_checkpoint == str(CHECKPOINT.resolve())
    assert len(command.sys.ifaces) == 1
    interface = command.sys.ifaces[0]
    view = interface.source_view
    assert view.selection_kind == 'interface'
    assert view.group_a_ids == tuple(sorted(interface.group1_indices))
    assert view.group_b_ids == tuple(sorted(interface.group2_indices))
    assert len(view.surfaces) == 191
    assert len(view.edges) == 512
    assert len(view.vertices) == 322
    assert view.completeness.status == 'INCOMPLETE'

    assert len(command.interface_results) == 1
    result = command.interface_results[0]
    assert result.status is Status.PARTIAL
    assert result.provenance.geometry_lineage.origin is GeometryLineageOrigin.MATERIALIZED_VIEW
    assert result.provenance.geometry_lineage.completeness['completeness_proven'] is False

    interface_dir = destination / '1bna_DB1976_no_sol' / interface.name
    for name in (
        'interface_surfs.off', 'interface_edges.off', 'interface_verts.off',
        'interface_surfaces.csv', 'interface_edges.csv', 'interface_vertices.csv',
        'group_a_atoms.pdb', 'group_b_atoms.pdb', 'interface_info.txt',
    ):
        assert (interface_dir / name).is_file(), name
    assert not list(destination.rglob('*.vpy'))
