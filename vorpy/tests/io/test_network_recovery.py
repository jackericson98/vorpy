"""Bounded tests for source-identity validation and checkpoint recovery."""

import hashlib
from pathlib import Path

import pandas as pd
import pytest

from vorpy.src.boundary import BoundaryConfig
from vorpy.src.boundary import BoundaryGenerator
from vorpy.src.io import (
    ArchiveError,
    group_from_network,
    load_network,
    recover_network_source,
    save_network,
)
from vorpy.src.interface import Interface
from vorpy.src.network.network import Network
from vorpy.src.system.system import System


CHECKPOINT = (
    Path(__file__).resolve().parents[3]
    / 'output' / 'benchmarks' / '1bna_DB1976' / 'pow_dry'
    / '1bna_DB1976_no_sol'
    / '1bna_DB1976_no_sol_full_system.geometry.vpy'
)


def test_actual_checkpoint_recovery_exposes_geometry_only_metadata():
    if not CHECKPOINT.is_file():
        pytest.skip('verified 1BNA checkpoint is not present in this checkout')

    network = load_network(CHECKPOINT)
    solve_universe_id = network.get_solve_context_key().solve_universe_id
    recovered = recover_network_source(
        network,
        expected_scheme='pow',
        expected_solvent_mode='MOLECULE_ONLY',
        expected_settings={'net_type': 'pow', 'solvent_mode': 'MOLECULE_ONLY'},
        expected_solve_universe_id=solve_universe_id,
    )
    identity = recovered.source_identity
    assert identity.generator_count == 805
    assert identity.generator_system_indices == tuple(range(805))
    assert identity.generator_stable_ids[-1] == 'system:804'
    assert len(identity.generators) == 805
    assert identity.generators[0].stable_id == 'system:0'
    assert identity.generators[0].coordinates == tuple(network.sys.balls.iloc[0]['loc'])
    assert identity.generators[0].radius == float(network.sys.balls.iloc[0]['rad'])
    assert identity.source_sha256 == network.archive_metadata['source_sha256']
    assert identity.geometry_status == 'GEOMETRY_COMPLETE'
    assert identity.analysis_status == 'NOT_CALCULATED'
    assert identity.solve_universe_id == solve_universe_id
    assert identity.source_network_id is None
    assert network.sys.ifaces == []
    assert identity.generator_policy['excluded_generator_count'] == 0


def test_recovery_rejects_a_mismatched_input_hash(tmp_path):
    if not CHECKPOINT.is_file():
        pytest.skip('verified 1BNA checkpoint is not present in this checkout')

    network = load_network(CHECKPOINT)
    source = tmp_path / 'mismatched.pdb'
    with open(network.archive_metadata['original_structure_path'], 'rb') as stream:
        original = stream.read()
    source.write_bytes(original + b'\nREMARK changed for recovery test\n')

    with pytest.raises(ArchiveError, match='source hash does not match'):
        recover_network_source(network, source)


def test_actual_checkpoint_view_preserves_order_and_partial_geometry_without_build(
    monkeypatch,
):
    if not CHECKPOINT.is_file():
        pytest.skip('verified 1BNA checkpoint is not present in this checkout')

    def forbidden_build(*args, **kwargs):
        raise AssertionError('source recovery must not rebuild geometry')

    monkeypatch.setattr(Network, 'build', forbidden_build)
    network = load_network(CHECKPOINT)
    recovered = recover_network_source(
        network,
        expected_scheme='pow',
        expected_solvent_mode='MOLECULE_ONLY',
    )
    group_a = [
        index for index, value in network.sys.balls['chain_name'].items()
        if value in {'A', 'B'}
    ]
    group_b = [
        index for index, value in network.sys.balls['chain_name'].items()
        if value == 'C'
    ]
    view = recovered.make_interface_view(group_a, group_b)

    assert view.group_a_ids == tuple(sorted(group_a))
    assert view.group_b_ids == tuple(sorted(group_b))
    assert len(view.surface_features) == 191
    assert len(view.edge_features) == 512
    assert len(view.vertex_features) == 322
    assert not view.completeness.proven


def _synthetic_source(path):
    lines = []
    for index, (x, y, z) in enumerate(
        ((0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 2.0, 0.0),
         (0.0, 0.0, 2.0), (4.0, 4.0, 4.0)),
        start=1,
    ):
        lines.append(
            f'ATOM  {index:5d}  C   GLY A{index:4d}    '
            f'{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C  '
        )
    path.write_text('\n'.join(lines) + '\nEND\n', encoding='ascii')
    return path


@pytest.mark.parametrize(
    ('scheme', 'solvent_mode'),
    [
        (scheme, solvent_mode)
        for scheme in ('aw', 'pow', 'prm')
        for solvent_mode in ('MOLECULE_ONLY', 'FULL_SOLVENT')
    ],
)
def test_synthetic_recovery_wraps_existing_network_and_keeps_order(
    tmp_path, monkeypatch, scheme, solvent_mode,
):
    source_path = _synthetic_source(tmp_path / 'synthetic.pdb')
    source_system = System(
        str(source_path), make_dir=False, print_actions=False,
        solvent_mode=solvent_mode,
    )

    locations = [tuple(row['loc']) for _, row in source_system.balls.iterrows()]
    balls = pd.DataFrame({
        'loc': locations,
        'rad': [float(row['rad']) for _, row in source_system.balls.iterrows()],
        'num': list(range(5)),
        'system_num': list(range(5)),
        'stable_id': [f'system:{index}' for index in range(5)],
        'is_boundary_generator': [False] * 5,
        'complete': [True] * 5,
    })
    verts = pd.DataFrame({
        'balls': [[0, 1, 4]] * 4,
        'complete': [True] * 4,
            'loc': [(0.2, 0.2, 0.2)] * 4,
        'rad': [0.2] * 4,
    })
    edges = pd.DataFrame({
        'balls': [[0, 1, 4]] * 4,
        'verts': [[0, 1, 2]] * 4,
        'points': [[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]] * 4,
        'complete': [True] * 4,
    })
    surfs = pd.DataFrame({
        'balls': [[0, 1, 4]],
        'verts': [[0, 1, 2]],
        'edges': [[0, 1, 2]],
        'complete': [True],
        'points': [[(0.2, 0.2, 0.2)] * 3],
        'tris': [[(0, 1, 2)]],
    })
    network = Network(
        locs=locations,
        rads=balls['rad'].tolist(),
        group=list(range(5)),
        settings={
            'net_type': scheme, 'surf_res': 0.2, 'max_vert': 5.0,
            'box_size': 1.25, 'solvent_mode': solvent_mode,
        },
        balls=balls, verts=verts, edges=edges, surfs=surfs,
        system=source_system,
        source_network_id=(
            f'producer:synthetic-{scheme}-{solvent_mode.lower()}'
        ),
        generator_system_indices=tuple(range(5)),
        generator_stable_ids=tuple(f'system:{index}' for index in range(5)),
        solvent_mode=solvent_mode,
    )
    network.boundary_config = BoundaryConfig(mode='none')
    network.boundary_mode = 'none'
    network.boundary_indices = ()
    network.boundary_generators = ()
    digest = hashlib.sha256(source_path.read_bytes()).hexdigest()
    network.archive_metadata = {
        'original_structure_path': str(source_path),
        'source_sha256': digest,
    }
    network.archive_provenance = {
        'geometry_status': 'GEOMETRY_COMPLETE',
        'analysis_status': 'NOT_CALCULATED',
    }

    def forbidden_build(*args, **kwargs):
        raise AssertionError('synthetic recovery must not build geometry')

    monkeypatch.setattr(Network, 'build', forbidden_build)
    producer_recovery = recover_network_source(
        network,
        expected_scheme=scheme,
        expected_solvent_mode=solvent_mode,
        expected_settings={
            'net_type': scheme,
            'surf_res': 0.2,
            'max_vert': 5.0,
            'box_size': 1.25,
            'solvent_mode': solvent_mode,
        },
    )
    expected_solve_universe_id = producer_recovery.source_identity.solve_universe_id
    expected_radii_digest = producer_recovery.source_identity.radii_digest
    with pytest.raises(ArchiveError, match='solve-universe identity'):
        recover_network_source(
            network,
            expected_scheme=scheme,
            expected_solvent_mode=solvent_mode,
            expected_solve_universe_id='solve:not-this-network',
        )
    checkpoint = save_network(network, tmp_path / 'synthetic_checkpoint.vpy')
    loaded = load_network(checkpoint)
    recovered = recover_network_source(
        loaded,
        expected_scheme=scheme,
        expected_solvent_mode=solvent_mode,
        expected_settings={'net_type': scheme, 'solvent_mode': solvent_mode},
        expected_solve_universe_id=expected_solve_universe_id,
        expected_radii_digest=expected_radii_digest,
    )
    view = recovered.make_interface_view([0], [1])

    assert recovered.source_identity.generator_count == 5
    assert recovered.source_identity.source_network_id == (
        f'producer:synthetic-{scheme}-{solvent_mode.lower()}'
    )
    assert recovered.source_identity.solve_universe_id == expected_solve_universe_id
    assert recovered.source_identity.geometry_checkpoint_valid is True
    assert recovered.source_identity.scientific_result_status == 'NOT_CALCULATED'
    assert recovered.source_identity.scientific_result_certified is False
    assert view.group_a_ids == (0,)
    assert view.group_b_ids == (1,)
    assert view.completeness.proven

    left = group_from_network(loaded, [0], 'left')
    right = group_from_network(loaded, [1], 'right')
    interface = Interface(loaded.sys, left, right)
    interface.net = loaded
    interface.source_network = loaded
    interface.source_view = view
    interface.group1_indices = {0}
    interface.group2_indices = {1}
    interface.dir = str(tmp_path / 'physical_interface')
    interface.export(surfs=True, edges=True, verts=True, info=True)
    for name in ('interface_surfs.off', 'interface_edges.off', 'interface_verts.off', 'interface_info.txt'):
        assert (tmp_path / 'physical_interface' / name).is_file()

    recovered_archive = save_network(loaded, tmp_path / 'recovered_interface.vpy')
    reloaded = load_network(recovered_archive)
    assert reloaded.archive_metadata['solve_universe_id'] == expected_solve_universe_id
    reloaded_recovery = recover_network_source(
        reloaded,
        expected_scheme=scheme,
        expected_solvent_mode=solvent_mode,
        expected_solve_universe_id=expected_solve_universe_id,
        expected_radii_digest=expected_radii_digest,
    )
    reloaded_view = reloaded_recovery.make_interface_view([0], [1])
    assert reloaded_view.selection_id == view.selection_id
    assert reloaded_view.completeness.proven
    assert reloaded_recovery.source_identity.geometry_checkpoint_valid is True
    assert reloaded_recovery.source_identity.scientific_result_status == 'NOT_CALCULATED'
    assert reloaded_recovery.source_identity.scientific_result_certified is False


def test_recovery_accepts_attested_boundary_rows_without_mixing_them_into_source(
    tmp_path,
):
    source_path = _synthetic_source(tmp_path / 'boundary-source.pdb')
    source_system = System(str(source_path), make_dir=False, print_actions=False)
    source_balls = source_system.balls.copy(deep=True)
    boundary = source_balls.iloc[[4]].copy()
    balls = pd.concat([source_balls, boundary], ignore_index=True)
    balls['system_num'] = [0, 1, 2, 3, 4, -1]
    balls['stable_id'] = [f'system:{index}' for index in range(5)] + ['boundary:5']
    balls['is_boundary_generator'] = [False] * 5 + [True]
    balls['complete'] = [True] * 6
    network = Network(
        locs=[tuple(value) for value in balls['loc']],
        rads=[float(value) for value in balls['rad']],
        group=list(range(5)),
        settings={
            'net_type': 'pow', 'surf_res': 0.2, 'max_vert': 5.0,
            'box_size': 1.25, 'solvent_mode': 'MOLECULE_ONLY',
        },
        balls=balls,
        verts=pd.DataFrame({'balls': [[0, 1, 4]], 'complete': [True]}),
        edges=pd.DataFrame({'balls': [[0, 1, 4]], 'verts': [[0]], 'complete': [True]}),
        surfs=pd.DataFrame({'balls': [[0, 1, 4]], 'verts': [[0]], 'edges': [[0]], 'complete': [True]}),
        system=source_system,
        generator_system_indices=tuple(range(5)) + (None,),
        generator_stable_ids=tuple(f'system:{index}' for index in range(5)) + ('boundary:5',),
        solvent_mode='MOLECULE_ONLY',
    )
    network.boundary_config = BoundaryConfig(mode='shell')
    network.boundary_mode = 'shell'
    network.boundary_indices = (5,)
    network.boundary_generators = (
        BoundaryGenerator(position=tuple(float(value) for value in balls.iloc[5]['loc'])),
    )
    network.archive_metadata = {
        'original_structure_path': str(source_path),
        'source_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
    }
    network.archive_provenance = {
        'geometry_status': 'GEOMETRY_COMPLETE',
        'analysis_status': 'NOT_CALCULATED',
    }

    recovered = recover_network_source(network)

    assert recovered.source_identity.generator_count == 5
    assert recovered.source_identity.generator_system_indices == tuple(range(5))
    assert recovered.source_identity.generator_stable_ids == tuple(
        f'system:{index}' for index in range(5)
    )
