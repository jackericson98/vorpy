"""Portable no-solve coverage for typed non-group InterfaceResult snapshots."""

import json
import hashlib
import zipfile

import pandas as pd
import pytest

from vorpy.src.group import Group
from vorpy.src.interface import Interface
from vorpy.src.io import ArchiveError, load_network, save_network
from vorpy.src.network import Network
from vorpy.src.results import (
    Environment,
    GeometryLineage,
    GeometryLineageOrigin,
    InterfaceResult,
    Provenance,
    Quantity,
    ResultIdentity,
    ResultsDependency,
    SharedNetworkIdentity,
    Status,
    Topology,
)
from vorpy.src.system import System


def _synthetic_network():
    system = System(
        spheres=pd.DataFrame({
            'loc': [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 2.0, 0.0)],
            'rad': [1.0, 1.0, 1.0],
            'mass': [1.0, 1.0, 1.0],
            'num': [0, 1, 2],
            'res_name': ['ALA', 'GLY', 'SER'],
            'name': ['CA', 'CA', 'CA'],
        }),
        name='typed-interface-archive',
        make_dir=False,
    )
    system.balls['res'] = [None, None, None]
    system.balls['chn'] = [None, None, None]
    balls = system.balls[['loc', 'rad', 'num']].copy()
    balls['complete'] = True
    balls['system_num'] = [0, 1, 2]
    empty = {
        'verts': pd.DataFrame(columns=['loc', 'rad', 'balls']),
        'edges': pd.DataFrame(columns=['balls', 'verts', 'points']),
        'surfs': pd.DataFrame(columns=['balls', 'verts', 'edges', 'points', 'tris']),
    }
    network = Network(
        locs=balls['loc'].tolist(),
        rads=balls['rad'].tolist(),
        group=[0, 1, 2],
        iface_grps=({0, 1}, {2}),
        group_name='typed-interface',
        settings={'net_type': 'pow'},
        balls=balls,
        **empty,
        sort_balls=False,
        system=system,
        source_network_id='producer:typed-interface-v1',
    )
    group_a = Group(system, name='A', atoms=[0, 1], settings={'net_type': 'pow'}, make_net=False)
    group_b = Group(system, name='B', atoms=[2], settings={'net_type': 'pow'}, make_net=False)
    group_a.group_id = 'A'
    group_b.group_id = 'B'
    interface = Interface(system, group_a, group_b, name='A_B', interface_id='A__B')
    interface.net = network
    interface.group1_indices = {0, 1}
    interface.group2_indices = {2}
    system.groups[:] = [group_a, group_b]
    system.ifaces[:] = [interface]
    return network, interface


def _typed_result(network, interface):
    archive_id = 'archive:typed-interface-v1'
    shared = SharedNetworkIdentity(
        source_network_id=network.source_network_id,
        source_identity_status='AUTHORITATIVE',
        solve_universe_id='solve:typed-interface-v1',
        frame_id='static',
        analysis_context_kind='GLOBAL_NETWORK',
        geometry_scheme='pow',
        partition_scheme='power',
        scope_mode='interface',
        solvent_mode='MOLECULE_ONLY',
        selection_id='view:typed-interface-v1',
        group_A_selection_id='view:A',
        group_B_selection_id='view:B',
        orientation='A -> B',
    )
    identity = ResultIdentity(
        result_kind='interface',
        system=network.sys.name,
        frame='static',
        network_id=archive_id,
        network_scope='interface',
        environment='dry',
        partition='power',
        representation='voronoi',
        radius_configuration_id='radii:typed-interface-v1',
        interface_id=interface.interface_id,
        group_A='A',
        group_B='B',
        side='shared',
        orientation='group1 -> group2',
        source_network_id=network.source_network_id,
        source_identity_status='AUTHORITATIVE',
        source_interface_id=interface.interface_id,
        solve_universe_id=shared.solve_universe_id,
        analysis_context_kind=shared.analysis_context_kind,
        selection_id=shared.selection_id,
        solvent_mode=shared.solvent_mode,
    )
    lineage = GeometryLineage(
        origin=GeometryLineageOrigin.DIRECT_SOLVE,
        solve_universe_id=shared.solve_universe_id,
        geometry_source_id='geometry:typed-interface-v1',
        geometry_version='pow-interface-v1',
        solver_settings_digest='settings:typed-interface-v1',
        solver_settings={'net_type': 'pow'},
        environment_digest='environment:dry-v1',
        environment={'network_environment': 'dry'},
        source_network_id=network.source_network_id,
        source_interface_id=interface.interface_id,
        selection_scope={'group_A_system_indices': (0, 1), 'group_B_system_indices': (2,)},
        orientation={'convention': 'physical_interface_A_to_B'},
        completeness={'completeness_proven': False},
    )
    dependency = ResultsDependency(
        shared,
        source_id=interface.interface_id,
        source_revision='producer:typed-interface-v1',
        status=Status.CERTIFIED,
    )
    provenance = Provenance(
        'synthetic typed producer',
        source_revision='producer:typed-interface-v1',
        details={
            'geometry_lineage': lineage.as_dict(),
            'dependency_digest': dependency.dependency_digest,
        },
    )
    result = InterfaceResult(
        identity,
        provenance,
        Status.PARTIAL,
        Environment('dry'),
        metrics={
            'area': Quantity(
                'area', 12.5, 'Å²', Status.CERTIFIED,
                'synthetic producer', 'interface', provenance=provenance,
            ),
        },
        topology=Topology(),
    )
    interface.typed_results = (result,)
    interface.results_dependency = dependency
    return result, dependency, lineage


def _archive_with_state(path, state):
    destination = path.with_name(path.stem + '-mutated.vpy')
    with zipfile.ZipFile(path) as source, zipfile.ZipFile(destination, 'w', zipfile.ZIP_DEFLATED) as target:
        for name in source.namelist():
            payload = json.dumps(state, separators=(',', ':')).encode() if name == 'state.json' else source.read(name)
            target.writestr(name, payload)
    return destination


def _interface_record(state):
    return next(record for record in state['objects'] if record['kind'] == 'interface')


def _encoded_map(encoded):
    assert set(encoded) == {'map'}
    return encoded['map']


def test_typed_interface_result_archive_round_trip_and_cache_export(tmp_path, monkeypatch):
    from vorpy.src.interface import Interface as InterfaceClass
    from vorpy.src.output import visualization

    network, interface = _synthetic_network()
    result, dependency, lineage = _typed_result(network, interface)
    assert not hasattr(network, 'archive_id')
    archive = save_network(network, tmp_path / 'typed.vpy')
    with zipfile.ZipFile(archive) as bundle:
        state = json.loads(bundle.read('state.json'))
    record = _interface_record(state)
    assert 'typed_results' in record['fields']
    assert 'results_dependency' in record['fields']
    group_record = next(item for item in state['objects'] if item['kind'] == 'group')
    assert 'typed_results' not in group_record['fields']
    assert 'results_dependency' not in group_record['fields']
    assert not any(item['kind'] == 'group_result' for item in state['objects'])

    monkeypatch.setattr(Network, 'build', lambda *args, **kwargs: pytest.fail('reload rebuilt topology'))
    monkeypatch.setattr(InterfaceClass, 'analyze_geometry', lambda *args, **kwargs: pytest.fail('reload rebuilt geometry'))
    restored = load_network(archive)
    restored_interface = restored.sys.ifaces[0]
    restored_result = restored_interface.typed_results[0]
    assert isinstance(restored_result, InterfaceResult)
    assert restored_result.identity.network_id == restored.archive_id
    assert restored_result.identity.source_network_id == network.source_network_id
    assert restored_result.identity.network_id != restored_result.identity.source_network_id
    assert restored_result.provenance.geometry_lineage == lineage
    assert restored_result.status is Status.PARTIAL
    assert restored_result.metrics['area'].status is Status.CERTIFIED
    assert restored_result.metrics['area'].value == pytest.approx(12.5)
    assert restored_interface.results_dependency.dependency_digest == dependency.dependency_digest

    output = tmp_path / 'cache-only' / 'logs.csv'
    visualization._compact_log(restored_interface, output)
    text = output.read_text(encoding='utf-8')
    assert '12.5' in text
    assert 'CERTIFIED' in text
    assert lineage.source_network_id in text


def test_interface_archive_metadata_records_source_hash_and_purpose(tmp_path):
    network, interface = _synthetic_network()
    _typed_result(network, interface)
    with pytest.raises(ValueError, match='Unknown archive purpose'):
        save_network(network, tmp_path / 'invalid-purpose.vpy', archive_purpose='other')
    source = tmp_path / 'source.pdb'
    source.write_bytes(b'ATOM synthetic source\n')
    network.sys.files['base_file'] = str(source)

    archive = save_network(
        network,
        tmp_path / 'checkpoint.vpy',
        archive_purpose='recovery_checkpoint',
    )
    with zipfile.ZipFile(archive) as bundle:
        metadata = json.loads(bundle.read('metadata.json'))
    expected_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    assert metadata['source_sha256'] == expected_hash
    assert metadata['archive_purpose'] == 'recovery_checkpoint'

    restored = load_network(archive)
    assert restored.archive_metadata['source_sha256'] == expected_hash
    assert restored.archive_metadata['archive_purpose'] == 'recovery_checkpoint'
    assert restored.sys.ifaces[0].typed_results[0].identity.source_network_id == (
        network.source_network_id
    )


def test_legacy_v3_archive_without_typed_fields_still_loads(tmp_path):
    network, interface = _synthetic_network()
    _typed_result(network, interface)
    archive = save_network(network, tmp_path / 'typed.vpy')
    with zipfile.ZipFile(archive) as bundle:
        state = json.loads(bundle.read('state.json'))
    record = _interface_record(state)
    record['fields'].pop('typed_results')
    record['fields'].pop('results_dependency')
    legacy = _archive_with_state(archive, state)
    restored = load_network(legacy)
    assert restored.sys.ifaces[0].typed_results == ()
    assert restored.sys.ifaces[0].results_dependency is None


@pytest.mark.parametrize('mutation,match', [
    ('missing_payload', 'snapshot payload'),
    ('bad_schema', 'snapshot schema version'),
])
def test_typed_snapshot_corruption_is_rejected(tmp_path, mutation, match):
    network, interface = _synthetic_network()
    _typed_result(network, interface)
    archive = save_network(network, tmp_path / 'typed.vpy')
    with zipfile.ZipFile(archive) as bundle:
        state = json.loads(bundle.read('state.json'))
    encoded_results = _interface_record(state)['fields']['typed_results']
    encoded_snapshot = encoded_results['tuple'][0]
    snapshot_map = _encoded_map(encoded_snapshot)
    if mutation == 'missing_payload':
        snapshot_map[:] = [item for item in snapshot_map if item[0] != 'payload']
    else:
        for item in snapshot_map:
            if item[0] == 'schema_version':
                item[1] = 99
    corrupted = _archive_with_state(archive, state)
    with pytest.raises(ArchiveError, match=match):
        load_network(corrupted)


def test_typed_result_owner_mismatch_is_rejected_before_save(tmp_path):
    network, interface = _synthetic_network()
    result, _, _ = _typed_result(network, interface)
    network.archive_id = result.identity.network_id
    wrong = object.__new__(InterfaceResult)
    object.__setattr__(wrong, 'identity', ResultIdentity(
        **{**result.identity.__dict__, 'network_id': 'archive:other'}
    ))
    for field in ('provenance', 'status', 'environment', 'metrics', 'topology',
                  'mean_curvature', 'gaussian_curvature', 'contacts', 'selection',
                  'extensions', 'schema_version'):
        object.__setattr__(wrong, field, getattr(result, field))
    interface.typed_results = (wrong,)
    with pytest.raises(ArchiveError, match='network identity'):
        save_network(network, tmp_path / 'wrong-owner.vpy')


def test_fail_if_exists_remains_atomic_for_typed_archive(tmp_path):
    network, interface = _synthetic_network()
    _typed_result(network, interface)
    archive = save_network(network, tmp_path / 'typed.vpy')
    original = archive.read_bytes()
    with pytest.raises(FileExistsError):
        save_network(network, archive, fail_if_exists=True)
    assert archive.read_bytes() == original
