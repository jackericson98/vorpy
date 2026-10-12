"""Real solved-network persistence and hostile/malformed archive regression tests."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import zipfile

import numpy as np
import pytest

from vorpy.src.io import ArchiveError, group_from_network, load_network, save_network
from vorpy.workbench.services.vorpy_backend import VorPyBackend, VorPySolveSettings


@pytest.fixture(scope='module')
def solved(tmp_path_factory):
    directory = tmp_path_factory.mktemp('archive-solve')
    source = directory / 'molecule.pdb'
    shutil.copy(Path(__file__).resolve().parents[2] / 'data' / 'EDTA.pdb', source)
    cwd = Path.cwd()
    try:
        result = VorPyBackend(VorPySolveSettings(network_type='pow')).solve(source, lambda *a: None, lambda: False)
        yield result
    finally:
        os.chdir(cwd)


@pytest.fixture
def archive(solved, tmp_path):
    path = tmp_path / 'network.vpy'
    save_network(solved.export_group.net, path)
    return path


def test_preset_archive_can_be_loaded_without_rebuilding(solved, tmp_path, monkeypatch):
    from vorpy.src.output.output import export_micro
    from vorpy.src.network import Network
    system = solved.export_group.sys
    monkeypatch.setitem(system.files, 'dir', str(tmp_path))
    for group in system.groups:
        monkeypatch.setattr(group, 'dir', group.dir)
    cwd = Path.cwd()
    try:
        export_micro(system)
        destination = tmp_path / (system.name + '.vpy')
        assert zipfile.is_zipfile(destination)
        with zipfile.ZipFile(destination) as archive:
            assert {'state.json', 'metadata.json'} <= set(archive.namelist())
            units = json.loads(archive.read('metadata.json'))['units']
            assert units['mean_curvature'] == 'angstrom^-1'
            assert units['gaussian_curvature'] == 'angstrom^-2'
            assert units['integrated_mean_curvature'] == 'angstrom'
            assert units['integrated_gaussian_curvature'] == '1'
        monkeypatch.setattr(Network, 'build', lambda *a, **kw: pytest.fail('Archive loader rebuilt topology'))
        loaded = load_network(destination)
        assert len(loaded.surfs) == len(solved.export_group.net.surfs)
        assert len(loaded.verts) == len(solved.export_group.net.verts)
    finally:
        os.chdir(cwd)


def same(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a:
            same(a[key], b[key])
    elif isinstance(a, (list, tuple, set, np.ndarray)):
        if isinstance(a, set):
            assert a == b
        else:
            assert len(a) == len(b)
            for left, right in zip(a, b):
                same(left, right)
    elif isinstance(a, (float, np.floating)):
        assert float(a) == pytest.approx(float(b), nan_ok=True)
    else:
        assert a == b


def _write_legacy_archive(source, destination, version):
    """Derive a schema-valid v1/v2 fixture from a current solved archive."""
    from vorpy.src.io.network_archive import LINKS, V1_FIELDS, V2_FIELDS

    allowed_fields = V1_FIELDS if version == 1 else V2_FIELDS
    with zipfile.ZipFile(source, 'r') as archive:
        members = {
            name: archive.read(name)
            for name in archive.namelist()
        }

    metadata = json.loads(members['metadata.json'])
    metadata['format_version'] = version
    state = json.loads(members['state.json'])
    records = [
        record for record in state['objects']
        if record['kind'] in allowed_fields
    ]
    old_to_new = {record['id']: index for index, record in enumerate(records)}

    def remap_references(value):
        if isinstance(value, list):
            return [remap_references(item) for item in value]
        if not isinstance(value, dict):
            return value
        if set(value) == {'ref'}:
            return {'ref': old_to_new[value['ref']]}
        return {name: remap_references(item) for name, item in value.items()}

    for index, record in enumerate(records):
        allowed = set(allowed_fields[record['kind']]) | set(LINKS[record['kind']])
        record['fields'] = {
            name: remap_references(value)
            for name, value in record['fields'].items()
            if name in allowed
        }
        record['id'] = index
        if version == 1:
            record.pop('schema_version', None)
    metadata['root_network'] = old_to_new[metadata['root_network']]
    state['objects'] = records

    # Object references also occur in table cells; remap every encoded table
    # member after dropping newer object kinds from the legacy projection.
    for name, content in tuple(members.items()):
        if name.startswith('tables/') and name.endswith('.json'):
            members[name] = json.dumps(
                remap_references(json.loads(content)), allow_nan=False
            ).encode('utf-8')

    members['metadata.json'] = json.dumps(metadata, allow_nan=False).encode('utf-8')
    members['state.json'] = json.dumps(state, allow_nan=False).encode('utf-8')
    with zipfile.ZipFile(destination, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
        for name, content in members.items():
            archive.writestr(name, content)


@pytest.mark.parametrize('version', [1, 2])
def test_legacy_archive_versions_load_without_rebuilding(solved, archive, tmp_path, monkeypatch, version):
    from vorpy.src.network import Network

    legacy = tmp_path / f'network-v{version}.vpy'
    _write_legacy_archive(archive, legacy, version)
    monkeypatch.setattr(Network, 'build', lambda *args, **kwargs: pytest.fail('Archive loader rebuilt topology'))

    loaded = load_network(legacy)

    original = solved.export_group.net
    assert len(loaded.balls) == len(original.balls)
    assert len(loaded.surfs) == len(original.surfs)
    assert len(loaded.edges) == len(original.edges)
    assert len(loaded.verts) == len(original.verts)
    assert loaded.settings == original.settings


def test_solve_context_key_is_json_safe_and_detects_source_changes(solved):
    from copy import deepcopy
    from dataclasses import asdict, fields
    from vorpy.src.network.aw_geometry_reuse import SolveUniverseKey

    network = solved.export_group.net
    options = {
        'frame_id': 'frame-v6',
        'environment': {'mode': 'dry'},
        'partition_scheme': 'power',
    }
    key = network.get_solve_context_key(**options)
    payload = asdict(key)
    assert set(payload) == {field.name for field in fields(SolveUniverseKey)}

    restored_payload = json.loads(json.dumps(payload, sort_keys=True, allow_nan=False))
    for name in ('physical_generator_ids', 'boundary_generator_ids'):
        restored_payload[name] = tuple(restored_payload[name])
    restored_payload['local_to_system'] = tuple(
        tuple(pair) for pair in restored_payload['local_to_system']
    )
    assert SolveUniverseKey(**restored_payload) == key
    assert key.solve_universe_id == f'solve:{key.digest}'
    assert key.context_id == key.digest

    geometry_changed = deepcopy(network)
    first_index = geometry_changed.balls.index[0]
    geometry_changed.balls.loc[first_index, 'rad'] = (
        float(geometry_changed.balls.loc[first_index, 'rad']) + 0.01
    )
    selection_changed = deepcopy(network)
    selected = tuple(selection_changed.group)
    assert len(selected) > 1
    selection_changed.group = selected[:-1]

    assert network.get_solve_context_key(**options).digest == key.digest
    assert geometry_changed.get_solve_context_key(**options).digest != key.digest
    assert selection_changed.get_solve_context_key(**options).digest != key.digest


def test_archive_round_trip_preserves_network_identity_and_context_mapping(
    solved, archive, monkeypatch
):
    from vorpy.src.network import Network

    original = solved.export_group.net
    before = original.get_solve_context_key()
    with zipfile.ZipFile(archive) as bundle:
        members = set(bundle.namelist())
        metadata = json.loads(bundle.read('metadata.json'))
        state = json.loads(bundle.read('state.json'))
    assert metadata['format_version'] == 3
    assert {'metadata.json', 'state.json'} <= members
    assert any(record['kind'] == 'network' for record in state['objects'])

    monkeypatch.setattr(Network, 'build', lambda *a, **kw: pytest.fail('Archive loader rebuilt topology'))
    loaded = load_network(archive)
    after = loaded.get_solve_context_key()

    assert loaded.archive_id == original.archive_id
    assert after == before
    assert after.local_to_system == before.local_to_system
    assert after.physical_generator_ids == before.physical_generator_ids


def test_v15_real_default_edta_source_network_identity_and_results_round_trip(
    tmp_path, monkeypatch
):
    from types import SimpleNamespace
    from vorpy.src.group.shared_view import attach_group_view_results, make_group_view
    from vorpy.src.io.network_archive import FIELDS, V1_FIELDS, V2_FIELDS
    from vorpy.src.io.scientific_adapters import stable_id
    from vorpy.src.boundary import BoundaryConfig
    from vorpy.src.group import Group
    from vorpy.src.network import Network
    from vorpy.src.network.views import NetworkViewError
    from vorpy.src.results import Status
    from vorpy.src.system import System

    edta = Path(__file__).resolve().parents[2] / 'data' / 'EDTA.pdb'
    system = System(
        file=str(edta), make_dir=False,
        boundary_config=BoundaryConfig('shell', shell_spacing=4.0),
    )
    group = Group(
        system, name='edta-default', atoms=list(range(32)), max_vert=5,
        make_net=True, build_net=False,
    )
    system.groups.append(group)
    builds = []
    original_build = Network.build

    def counted_build(*args, **kwargs):
        builds.append(True)
        return original_build(*args, **kwargs)

    monkeypatch.setattr(Network, 'build', counted_build)
    group.build()
    group.get_layers(max_layers=50)
    network = group.net
    assert len(builds) == 1
    source_network_id = 'producer:edta-default'
    monkeypatch.setattr(network, 'source_network_id', source_network_id)
    assert 'source_network_id' in FIELDS['network']
    assert 'source_network_id' not in V1_FIELDS['network']
    assert 'source_network_id' not in V2_FIELDS['network']

    before_context = network.get_solve_context_key()
    assert network.solvent_mode == 'MOLECULE_ONLY'
    before_solvent = network.solvent_provenance
    before_retention = network.retained_generator_provenance
    archive = save_network(network, tmp_path / 'edta-default-source-id.vpy')
    archive_id = network.archive_id
    with zipfile.ZipFile(archive) as bundle:
        state = json.loads(bundle.read('state.json'))
    record = next(item for item in state['objects'] if item['kind'] == 'network')
    assert record['fields']['source_network_id'] == source_network_id
    assert archive_id != source_network_id

    monkeypatch.setattr(
        Network, 'build', lambda *args, **kwargs: pytest.fail('Archive path rebuilt Network')
    )
    context = network.get_solve_context_key()
    groups = []
    for _, system_id in context.local_to_system:
        if system_id is None:
            continue
        group = SimpleNamespace(
            ball_ndxs=(system_id,), group_id=f'edta-default-{system_id}'
        )
        try:
            make_group_view(group, network)
        except NetworkViewError:
            continue
        groups.append(group)
        if len(groups) == 2:
            break
    assert len(groups) == 2
    groups = tuple(groups)
    before_results = attach_group_view_results(
        SimpleNamespace(), network, groups,
        system='EDTA', radius_configuration_id='edta-default-radii-v1',
        source_revision=archive_id,
    )

    restored = load_network(archive)
    after_results = attach_group_view_results(
        SimpleNamespace(), restored, groups,
        system='EDTA', radius_configuration_id='edta-default-radii-v1',
        source_revision=restored.archive_id,
    )

    assert restored.source_network_id == source_network_id
    assert restored.archive_id == archive_id and restored.archive_id != source_network_id
    after_context = restored.get_solve_context_key()
    assert after_context == before_context
    assert after_context.boundary_generator_digest == before_context.boundary_generator_digest
    assert after_context.environment_digest == before_context.environment_digest
    assert restored.solvent_provenance == before_solvent
    assert restored.retained_generator_provenance == before_retention
    assert restored.generator_system_indices == network.generator_system_indices
    assert restored.generator_stable_ids == network.generator_stable_ids
    for before, after in zip(before_results, after_results, strict=True):
        assert before.status is after.status is Status.CERTIFIED
        assert before.identity.source_network_id == after.identity.source_network_id == source_network_id
        assert before.identity.solve_universe_id == after.identity.solve_universe_id
        assert before.identity.selection_id == after.identity.selection_id
        assert before.provenance.details['dependency_digest'] == after.provenance.details['dependency_digest']
        assert after.provenance.details['completeness']['proven']
        assert after.metrics['area'].value == pytest.approx(before.metrics['area'].value)
        assert after.metrics['volume'].value == pytest.approx(before.metrics['volume'].value)
    assert before_results[0].identity.selection_id != before_results[1].identity.selection_id
    assert len(builds) == 1


def test_v15_legacy_archive_without_source_network_id_stays_absent(archive):
    def omit_source_network_id(state):
        for record in state['objects']:
            if record['kind'] == 'network':
                record['fields'].pop('source_network_id', None)
        return state

    rewrite(archive, 'state.json', omit_source_network_id)
    restored = load_network(archive)
    assert restored.source_network_id is None


def test_v15_missing_authoritative_source_id_does_not_certify_results(
    solved, tmp_path, monkeypatch
):
    from vorpy.src.network import Network
    from vorpy.src.network.views import _resolve_source_network_id
    from vorpy.src.results import (
        ResultsDependency,
        SharedNetworkIdentity,
        Status,
        validate_results_dependency,
    )

    network = solved.export_group.net
    monkeypatch.setattr(network, 'source_network_id', None)
    archive = save_network(network, tmp_path / 'edta-missing-source-id.vpy')
    monkeypatch.setattr(
        Network, 'build', lambda *args, **kwargs: pytest.fail('Archive path rebuilt Network')
    )
    restored = load_network(archive)
    assert restored.source_network_id is None

    context = restored.get_solve_context_key()
    resolved_source_id = _resolve_source_network_id(restored, context)
    identity = SharedNetworkIdentity(
        source_network_id=resolved_source_id,
        source_identity_status="PROVISIONAL",
        solve_universe_id=context.solve_universe_id,
        frame_id=context.frame_id,
        analysis_context_kind='GLOBAL_NETWORK',
        geometry_scheme=context.scheme,
        partition_scheme=context.partition_scheme,
        scope_mode=context.scope_mode,
        selection_id='view:missing-source-regression',
    )
    dependency = ResultsDependency(
        identity,
        'missing-source-regression',
        source_revision=restored.archive_id,
    )
    assert validate_results_dependency(dependency, dependency) is Status.UNRESOLVED

    archive_identity = SharedNetworkIdentity(
        source_network_id=restored.archive_id,
        solve_universe_id=context.solve_universe_id,
        frame_id=context.frame_id,
        analysis_context_kind='GLOBAL_NETWORK',
        geometry_scheme=context.scheme,
        partition_scheme=context.partition_scheme,
        scope_mode=context.scope_mode,
        selection_id='view:archive-id-only-regression',
    )
    archive_dependency = ResultsDependency(
        archive_identity,
        'archive-id-only-regression',
        source_revision=restored.archive_id,
    )
    assert archive_identity.source_identity_status == 'UNRESOLVED'
    assert validate_results_dependency(
        archive_dependency, archive_dependency
    ) is Status.UNRESOLVED


def test_v15_archive_rejects_malformed_source_network_id(solved, tmp_path, monkeypatch):
    network = solved.export_group.net
    monkeypatch.setattr(network, 'source_network_id', 'producer:valid-id')
    archive = save_network(network, tmp_path / 'edta-source-id-valid.vpy')

    monkeypatch.setattr(network, 'source_network_id', '  ')
    with pytest.raises(ArchiveError, match='network.source_network_id'):
        save_network(network, tmp_path / 'edta-source-id-invalid.vpy')

    rewrite(
        archive,
        'state.json',
        lambda state: next(
            record['fields'].__setitem__('source_network_id', '') or state
            for record in state['objects']
            if record['kind'] == 'network'
        ),
    )
    with pytest.raises(ArchiveError, match='network.source_network_id'):
        load_network(archive)


def test_real_edta_partial_view_fails_closed_without_certification(solved, monkeypatch):
    """A partial real solve must reject an interface view, not certify it."""
    from vorpy.src.network import Network
    from vorpy.src.network.views import NetworkViewError

    original = solved.export_group.net
    context = original.get_solve_context_key()
    local_to_system = dict(context.local_to_system)
    source_balls = original.surfs.iloc[0].get('balls', ())
    physical_ids = [
        local_to_system[int(local)]
        for local in source_balls
        if local_to_system[int(local)] is not None
    ]
    assert len(physical_ids) >= 2
    group_a, group_b = physical_ids[:1], physical_ids[1:2]
    monkeypatch.setattr(Network, 'build', lambda *a, **kw: pytest.fail('Archive loader rebuilt topology'))
    with pytest.raises(NetworkViewError, match='Source cell completeness is not proven'):
        original.make_view(iface_grps=(group_a, group_b))


def test_v8_archive_to_group_view_round_trip_and_fail_closed_identity(
    solved, tmp_path, monkeypatch
):
    from copy import copy
    from dataclasses import replace
    from vorpy.src.network import Network
    from vorpy.src.results import (
        ResultsDependency,
        SharedNetworkIdentity,
        Status,
        validate_results_dependency,
    )

    from vorpy.tests.network.test_views import _network as complete_solved_network

    original = complete_solved_network(scheme='prm')
    original.sys = copy(solved.export_group.sys)
    original.sys.groups = []
    original.sys.ifaces = []
    original.sys.data = []
    original.sys.interface_geometry_cache = {'verts': {}, 'edges': {}, 'surfs': {}}
    original.sys._interface_phase1_dry_cache = {}
    original.sys.residues = []
    original.sys.chains = []
    original.sys.sol = None
    original.sys.balls = original.sys.balls[
        original.sys.balls['num'].isin((10, 11, 12, 13))
    ].copy(deep=True)
    for column in ('sys', 'res', 'chn', 'chain', 'residue'):
        if column in original.sys.balls:
            original.sys.balls[column] = None
    original.group = [0, 1, 2, 3]
    original.boundary_generators = ()
    monkeypatch.setattr(original, 'archive_id', 'vector-v8-source-network', raising=False)
    archive = save_network(original, tmp_path / 'v8.vpy')
    context = original.get_solve_context_key()
    group_id = next(system_id for _, system_id in context.local_to_system if system_id == 10)
    view = original.make_view(group_ids=[group_id])
    assert view.selection_kind == 'group'
    assert view.group_a_ids == view.group_b_ids == ()
    assert view.completeness.proven
    assert view.completeness.status == 'COMPLETE'

    with zipfile.ZipFile(archive) as bundle:
        metadata = json.loads(bundle.read('metadata.json'))
        state = json.loads(bundle.read('state.json'))
    assert metadata['format_version'] == 3
    assert 'network_view' not in {record['kind'] for record in state['objects']}

    monkeypatch.setattr(Network, 'build', lambda *a, **kw: pytest.fail('Archive loader rebuilt topology'))
    loaded = load_network(archive)
    restored = loaded.make_view(group_ids=[group_id])

    assert restored.source_network_id == original.source_network_id
    assert loaded.archive_id == original.archive_id
    assert restored.source_network_id != loaded.archive_id
    assert restored.solve_universe_id == view.solve_universe_id
    assert restored.frame_id == view.frame_id == context.frame_id
    assert restored.selection_id == view.selection_id
    assert restored.local_to_system == view.local_to_system
    assert restored.completeness == view.completeness
    assert tuple(
        (feature.generator_local_indices, feature.generator_ids,
         feature.generator_system_ids)
        for feature in restored.cells
    ) == tuple(
        (feature.generator_local_indices, feature.generator_ids,
         feature.generator_system_ids)
        for feature in view.cells
    )

    identity = SharedNetworkIdentity(
        source_network_id=restored.source_network_id,
        source_identity_status='AUTHORITATIVE',
        solve_universe_id=restored.solve_universe_id,
        frame_id=restored.frame_id,
        analysis_context_kind='GLOBAL_NETWORK',
        geometry_scheme=restored.scheme,
        partition_scheme=restored.solve_context.partition_scheme,
        scope_mode=restored.solve_context.scope_mode,
        selection_id=restored.selection_id,
    )
    expected = ResultsDependency(identity, 'group-cache', source_revision='v8')
    assert validate_results_dependency(expected, expected) is Status.CERTIFIED
    assert validate_results_dependency(expected, None) is Status.NOT_CALCULATED
    assert validate_results_dependency(
        expected,
        ResultsDependency(replace(identity, selection_id='view:stale-selection'),
                          'group-cache', source_revision='v8'),
    ) is Status.STALE
    assert validate_results_dependency(
        expected,
        ResultsDependency(
            replace(identity, source_network_id=f'network:{identity.solve_universe_id.removeprefix("solve:")}',
                    source_identity_status='PROVISIONAL'),
            'group-cache', source_revision='v8',
        ),
    ) is Status.UNRESOLVED
    with pytest.raises(ValueError, match='canonical view'):
        replace(identity, selection_id='invalid-selection')


def test_v9_boundary_generators_round_trip_and_legacy_absence_stays_unresolved(
    tmp_path, monkeypatch
):
    from vorpy.src.boundary import BoundaryConfig, BoundaryGenerator
    from vorpy.src.io.network_archive import FIELDS, V1_FIELDS, V2_FIELDS
    from vorpy.src.network import Network
    from vorpy.src.system import System
    from vorpy.tests.network.test_views import _network as complete_network

    original = complete_network(scheme='prm')
    system = System(simple=True, make_dir=False)
    system.frame_index = 7
    system.balls = original.balls.iloc[:4][['loc', 'rad', 'num']].copy(deep=True)
    system.balls['res'] = None
    system.balls['chn'] = None
    system.groups = []
    system.ifaces = []
    system.boundary_config = BoundaryConfig(mode='shell')
    system.boundary_mode = 'shell'
    system.boundary_generators = ()
    original.sys = system
    original.group = [0, 1, 2, 3]
    original.boundary_config = BoundaryConfig(mode='shell')
    original.boundary_mode = 'shell'
    original.boundary_indices = (4,)
    original.boundary_generators = (
        BoundaryGenerator(position=(4.0, 4.0, 4.0), radius=1.0),
    )
    monkeypatch.setattr(original, 'archive_id', 'vector-v9-boundary-network', raising=False)
    assert 'boundary_generators' in FIELDS['network']
    assert 'boundary_generators' not in V1_FIELDS['network']
    assert 'boundary_generators' not in V2_FIELDS['network']

    archive = save_network(original, tmp_path / 'v9-boundary.vpy')
    with zipfile.ZipFile(archive, 'r') as bundle:
        metadata = json.loads(bundle.read('metadata.json'))
        records = json.loads(bundle.read('state.json'))['objects']
    network_record = next(record for record in records if record['kind'] == 'network')
    assert metadata['format_version'] == 3
    assert 'boundary_generators' in network_record['fields']
    assert any(record['kind'] == 'boundary_generator' for record in records)
    before = original.get_solve_context_key()
    view = original.make_view(group_ids=[10])

    monkeypatch.setattr(Network, 'build', lambda *a, **kw: pytest.fail('Archive loader rebuilt topology'))
    loaded = load_network(archive)
    after = loaded.get_solve_context_key()
    restored = loaded.make_view(group_ids=[10])

    assert after == before
    assert loaded.boundary_generators == original.boundary_generators
    assert restored.source_network_id == original.source_network_id
    assert loaded.archive_id == original.archive_id
    assert restored.source_network_id != loaded.archive_id
    assert restored.solve_universe_id == view.solve_universe_id
    assert restored.frame_id == view.frame_id
    assert restored.selection_id == view.selection_id
    assert restored.local_to_system == view.local_to_system
    assert restored.boundary_generator_ids == view.boundary_generator_ids == ('boundary:4',)
    assert restored.completeness == view.completeness

    legacy = tmp_path / 'v3-missing-boundary.vpy'
    with zipfile.ZipFile(archive, 'r') as source:
        members = {name: source.read(name) for name in source.namelist()}
    state = json.loads(members['state.json'])
    for record in state['objects']:
        if record['kind'] == 'network':
            record['fields'].pop('boundary_generators', None)
    members['state.json'] = json.dumps(state, allow_nan=False).encode('utf-8')
    with zipfile.ZipFile(legacy, 'w', compression=zipfile.ZIP_DEFLATED) as target:
        for name, content in members.items():
            target.writestr(name, content)

    legacy_loaded = load_network(legacy)
    assert not hasattr(legacy_loaded, 'boundary_generators')
    legacy_context = legacy_loaded.get_solve_context_key()
    assert legacy_context.boundary_generator_ids == ('boundary:4',)
    assert legacy_context.boundary_generator_digest != before.boundary_generator_digest


def test_v10_solvent_provenance_round_trip_and_legacy_absence_stays_unresolved(
    tmp_path, monkeypatch
):
    from vorpy.src.boundary import BoundaryConfig, SolventMode
    from vorpy.src.io.network_archive import FIELDS, V1_FIELDS, V2_FIELDS
    from vorpy.src.network import Network
    from vorpy.src.system import System
    from vorpy.tests.network.test_views import _network as complete_network

    original = complete_network(scheme='prm')
    system = System(simple=True, make_dir=False, solvent_mode=SolventMode.SOLVENT_CONTEXT)
    system.frame_index = 8
    system.balls = original.balls.iloc[:4][['loc', 'rad', 'num']].copy(deep=True)
    system.balls['res'] = None
    system.balls['chn'] = None
    system.groups = []
    system.ifaces = []
    system.boundary_config = BoundaryConfig(mode='shell')
    system.boundary_mode = 'shell'
    system.boundary_generators = ()
    original.sys = system
    original.group = [0, 1, 2, 3]
    original.generator_system_indices = (0, 1, 2, 3, None)
    original.generator_stable_ids = ('atom-A', 'atom-B', 'atom-C', 'atom-D', 'boundary:4')
    original.balls['system_num'] = [0, 1, 2, 3, -1]
    original.solvent_local_indices = (1,)
    original.solvent_generator_ids = ('atom-B',)
    original.solvent_mode = SolventMode.SOLVENT_CONTEXT.value
    original.balls['is_solvent_generator'] = [False, True, False, False, False]
    original.refresh_solvent_provenance()
    monkeypatch.setattr(original, 'archive_id', 'vector-v10-solvent-network', raising=False)

    assert 'solvent_mode' in FIELDS['system']
    assert 'solvent_generator_digest' in FIELDS['network']
    assert 'solvent_mode' not in V1_FIELDS['network']
    assert 'solvent_mode' not in V2_FIELDS['network']
    before_context = original.get_solve_context_key()
    before_provenance = original.solvent_provenance
    archive = save_network(original, tmp_path / 'v10-solvent.vpy')
    with zipfile.ZipFile(archive, 'r') as bundle:
        records = json.loads(bundle.read('state.json'))['objects']
    network_record = next(record for record in records if record['kind'] == 'network')
    system_record = next(record for record in records if record['kind'] == 'system')
    assert network_record['fields']['solvent_mode'] == 'SOLVENT_CONTEXT'
    assert system_record['fields']['solvent_mode'] == 'SOLVENT_CONTEXT'
    for field in (
        'generator_system_indices', 'generator_stable_ids',
        'solvent_local_indices', 'solvent_generator_ids',
        'solvent_cells_complete', 'solvent_generator_digest',
    ):
        assert field in network_record['fields']

    monkeypatch.setattr(Network, 'build', lambda *a, **kw: pytest.fail('Archive loader rebuilt topology'))
    loaded = load_network(archive)
    assert loaded.sys.solvent_mode is SolventMode.SOLVENT_CONTEXT
    assert loaded.solvent_provenance == before_provenance
    assert loaded.generator_system_indices == original.generator_system_indices
    assert loaded.generator_stable_ids == original.generator_stable_ids
    assert loaded.get_solve_context_key() == before_context

    legacy = tmp_path / 'v3-missing-solvent.vpy'
    with zipfile.ZipFile(archive, 'r') as source:
        members = {name: source.read(name) for name in source.namelist()}
    state = json.loads(members['state.json'])
    for record in state['objects']:
        if record['kind'] == 'network':
            for field in (
                'generator_system_indices', 'generator_stable_ids',
                'solvent_local_indices', 'solvent_generator_ids', 'solvent_mode',
                'solvent_cells_complete', 'solvent_generator_digest',
            ):
                record['fields'].pop(field, None)
        elif record['kind'] == 'system':
            record['fields'].pop('solvent_mode', None)
    members['state.json'] = json.dumps(state, allow_nan=False).encode('utf-8')
    with zipfile.ZipFile(legacy, 'w', compression=zipfile.ZIP_DEFLATED) as target:
        for name, content in members.items():
            target.writestr(name, content)

    legacy_loaded = load_network(legacy)
    assert legacy_loaded.solvent_provenance is None
    assert legacy_loaded.solvent_mode is None
    assert legacy_loaded.sys.solvent_mode is None


def test_real_f10_solvent_provenance_round_trip_without_rebuild(solved, tmp_path, monkeypatch):
    from vorpy.src.network import Network

    network = solved.export_group.net
    assert network.solvent_mode == 'MOLECULE_ONLY'
    assert network.solvent_provenance is not None
    assert len(network.generator_system_indices) == len(network.balls)
    before_context = network.get_solve_context_key()
    before = network.solvent_provenance

    archive = save_network(network, tmp_path / 'real-f10-solvent.vpy')
    monkeypatch.setattr(Network, 'build', lambda *a, **kw: pytest.fail('Archive loader rebuilt topology'))
    loaded = load_network(archive)

    assert loaded.solvent_provenance == before
    assert loaded.sys.solvent_mode.value == 'MOLECULE_ONLY'
    assert loaded.get_solve_context_key() == before_context
    assert network.get_solve_context_key() == before_context
    restored_context = loaded.get_solve_context_key()
    assert restored_context.boundary_generator_digest == before_context.boundary_generator_digest
    assert restored_context.environment_digest == before_context.environment_digest
    assert loaded.archive_id == network.archive_id
    assert loaded.generator_system_indices == network.generator_system_indices
    assert loaded.generator_stable_ids == network.generator_stable_ids
    assert loaded.solvent_cells_complete == network.solvent_cells_complete


def test_archive_rejects_malformed_physical_generator_id(tmp_path):
    from vorpy.src.io.network_archive import ArchiveError
    from vorpy.src.system import System
    from vorpy.tests.network.test_views import _network as complete_network

    network = complete_network(scheme='prm')
    network.sys = System(simple=True, make_dir=False)
    network.sys.balls = network.balls.iloc[:4][['loc', 'rad', 'num']].copy(deep=True)
    network.sys.balls['res'] = None
    network.sys.balls['chn'] = None
    network.sys.groups = []
    network.sys.ifaces = []
    network.group = [0, 1, 2, 3]
    network.generator_system_indices = (-1, 1, 2, 3, None)
    network.generator_stable_ids = ('atom-A', 'atom-B', 'atom-C', 'atom-D', 'boundary:4')
    network.balls['system_num'] = [-1, 1, 2, 3, -1]

    with pytest.raises(ArchiveError, match='invalid physical generator ID'):
        save_network(network, tmp_path / 'malformed-generator.vpy')


def test_archive_rejects_invalid_boundary_metadata(tmp_path):
    from vorpy.src.io.network_archive import ArchiveError
    from vorpy.src.system import System
    from vorpy.tests.network.test_views import _network as complete_network

    network = complete_network(scheme='prm')
    network.sys = System(simple=True, make_dir=False)
    network.sys.balls = network.balls.iloc[:4][['loc', 'rad', 'num']].copy(deep=True)
    network.sys.balls['res'] = None
    network.sys.balls['chn'] = None
    network.sys.groups = []
    network.sys.ifaces = []
    network.group = [0, 1, 2, 3]
    network.boundary_indices = (len(network.balls),)

    with pytest.raises(ArchiveError, match='boundary_indices.*invalid/dangling ID'):
        save_network(network, tmp_path / 'invalid-boundary.vpy')


def test_archive_preserves_boundary_layer_atom_sentinel(solved, tmp_path):
    complete_system_id = next(
        int(system_id)
        for system_id, complete in zip(
            solved.export_group.net.balls['system_num'],
            solved.export_group.net.balls['complete'],
        )
        if int(system_id) >= 0 and bool(complete)
    )
    selection = group_from_network(
        solved.export_group.net, [complete_system_id], 'boundary-layer-sentinel'
    )
    archive = save_network(
        solved.export_group.net, tmp_path / 'layer-sentinel-legacy-source.vpy'
    )
    legacy = tmp_path / 'layer-sentinel-legacy.vpy'
    shutil.copyfile(archive, legacy)

    def add_historical_sentinel(state):
        record = next(
            record for record in state['objects']
            if record['kind'] == 'group'
            and record['fields'].get('name') == 'boundary-layer-sentinel'
        )
        record['fields']['layer_atoms'] = {'list': [{'list': [complete_system_id, -1]}]}
        return state

    rewrite(legacy, 'state.json', add_historical_sentinel)
    rewrite(
        legacy, 'metadata.json',
        lambda metadata: {
            key: value for key, value in metadata.items()
            if key != 'producer_contracts'
        },
    )
    loaded = load_network(legacy)
    restored = next(
        item for item in loaded.sys.groups if item.name == 'boundary-layer-sentinel'
    )
    assert restored.layer_atoms == [[complete_system_id, -1]]
    assert loaded.archive_metadata['format_version'] == 3
    assert 'producer_contracts' not in loaded.archive_metadata


def test_new_archives_reject_boundary_layer_sentinel(solved, tmp_path):
    from vorpy.src.io.network_archive import GROUP_LAYER_ATOMS_CONTRACT

    network = solved.export_group.net
    archive = save_network(network, tmp_path / 'marked-valid.vpy')
    with zipfile.ZipFile(archive) as bundle:
        metadata = json.loads(bundle.read('metadata.json'))
    assert metadata['producer_contracts']['group.layer_atoms'] == GROUP_LAYER_ATOMS_CONTRACT

    complete_system_id = next(
        int(system_id)
        for system_id, complete in zip(network.balls['system_num'], network.balls['complete'])
        if int(system_id) >= 0 and bool(complete)
    )
    selection = group_from_network(network, [complete_system_id], 'invalid-new-layer')
    try:
        selection.layer_atoms = [[complete_system_id, -1]]
        with pytest.raises(ArchiveError, match='group.layer_atoms.*-1'):
            save_network(network, tmp_path / 'invalid-new-layer.vpy')
        selection.layer_atoms = [[complete_system_id, True]]
        with pytest.raises(ArchiveError, match='group.layer_atoms.*True'):
            save_network(network, tmp_path / 'boolean-new-layer.vpy')
    finally:
        network.sys.groups.remove(selection)

    def corrupt_current_archive(state):
        group = next(record for record in state['objects'] if record['kind'] == 'group')
        group['fields']['layer_atoms'] = {'list': [{'list': [-1]}]}
        return state

    rewrite(archive, 'state.json', corrupt_current_archive)
    with pytest.raises(ArchiveError, match='group.layer_atoms.*-1'):
        load_network(archive)


def test_unknown_group_layer_producer_contract_fails_closed(archive):
    rewrite(
        archive, 'metadata.json',
        lambda metadata: {
            **metadata,
            'producer_contracts': {'group.layer_atoms': 'future_unknown_contract'},
        },
    )
    with pytest.raises(ArchiveError, match='Unsupported group.layer_atoms producer contract'):
        load_network(archive)


def test_incomplete_group_layer_producer_contract_fails_closed(archive):
    rewrite(
        archive, 'metadata.json',
        lambda metadata: {**metadata, 'producer_contracts': {}},
    )
    with pytest.raises(ArchiveError, match='Missing group.layer_atoms producer contract'):
        load_network(archive)


def test_real_f13_retained_edta_archive_reconstructs_two_results(
    tmp_path, monkeypatch
):
    from dataclasses import FrozenInstanceError
    from types import SimpleNamespace
    from vorpy.src.boundary import BoundaryConfig
    from vorpy.src.group import Group
    from vorpy.src.group.shared_view import attach_group_view_results, make_group_view
    from vorpy.src.io.scientific_adapters import stable_id
    from vorpy.src.network import Network
    from vorpy.src.network.views import NetworkViewError
    from vorpy.src.results import Status
    from vorpy.src.system import System

    edta = Path(__file__).resolve().parents[2] / 'data' / 'EDTA.pdb'
    system = System(
        file=str(edta), make_dir=False, retained_system_indices=(32,),
        boundary_config=BoundaryConfig('shell', shell_spacing=4.0),
    )
    group = Group(
        system, name='edta-retained-mg', atoms=list(range(32)), max_vert=5,
        make_net=True, build_net=False,
    )
    system.groups.append(group)
    builds = []
    original_build = Network.build

    def counted_build(*args, **kwargs):
        builds.append(True)
        return original_build(*args, **kwargs)

    monkeypatch.setattr(Network, 'build', counted_build)
    group.build()
    group.get_layers(max_layers=50)
    network = group.net
    assert len(builds) == 1
    producer_source_id = 'producer:edta-retained-mg'
    monkeypatch.setattr(network, 'source_network_id', producer_source_id)
    assert network.solvent_mode == 'MOLECULE_ONLY'
    assert network.retained_system_indices == (32,)
    assert network.retained_generator_ids == ('system:32',)
    assert network.retained_generator_digest

    boundary_local = set(network.boundary_indices)
    boundary_topology = {
        int(network.balls.iloc[index]['num']) for index in boundary_local
    }
    layer_system_ids = [
        int(value) for layer in group.layer_atoms for value in layer
    ]
    layer_topology_ids = {
        int(value) for layer in group.layer_net_atoms for value in layer
    }
    assert layer_system_ids and all(0 <= value < len(system.balls) for value in layer_system_ids)
    assert -1 not in layer_system_ids
    assert boundary_topology.intersection(layer_topology_ids)
    assert all(network.generator_system_indices[index] is None for index in boundary_local)
    assert all(bool(network.balls.iloc[index]['is_boundary_generator']) for index in boundary_local)

    archive_id = stable_id(network)
    assert archive_id != producer_source_id
    before_context = network.get_solve_context_key()
    before_solvent = network.solvent_provenance
    before_retention = network.retained_generator_provenance
    selected_system_ids = []
    before_views = []
    incomplete_mg_rejections = 0
    for system_id in range(32):
        candidate = SimpleNamespace(
            ball_ndxs=[system_id], group_id=f'edta-{system_id}'
        )
        try:
            view = make_group_view(candidate, network)
        except NetworkViewError as error:
            if str(error) != 'Source cell completeness is not proven for generator 32.':
                raise
            incomplete_mg_rejections += 1
            continue
        selected_system_ids.append(system_id)
        before_views.append(view)
        if len(selected_system_ids) == 2:
            break
    assert len(selected_system_ids) == 2
    assert incomplete_mg_rejections > 0
    selections = tuple(
        SimpleNamespace(ball_ndxs=[system_id], group_id=f'edta-{system_id}')
        for system_id in selected_system_ids
    )
    before_views = tuple(before_views)
    assert before_views[0].source_network_id == before_views[1].source_network_id == producer_source_id
    assert before_views[0].solve_universe_id == before_views[1].solve_universe_id
    assert before_views[0].selection_id != before_views[1].selection_id
    assert all(view.completeness.proven for view in before_views)
    with pytest.raises(FrozenInstanceError):
        before_views[0].selection_id = 'mutated'

    before_results = attach_group_view_results(
        SimpleNamespace(), network, selections,
        system='EDTA', radius_configuration_id='edta-retained-radii-v1',
        source_revision=archive_id,
    )
    assert all(item.status is Status.CERTIFIED for item in before_results)
    assert all(
        item.provenance.details['source_retained_generator_provenance']
        == before_retention
        for item in before_results
    )

    archive = save_network(network, tmp_path / 'f13-edta-retained.vpy')
    with zipfile.ZipFile(archive) as bundle:
        metadata = json.loads(bundle.read('metadata.json'))
        state = json.loads(bundle.read('state.json'))
    assert metadata['producer_contracts']['group.layer_atoms'] == 'system_global_nonnegative_v1'
    network_record = next(record for record in state['objects'] if record['kind'] == 'network')
    assert network_record['fields']['source_network_id'] == producer_source_id

    restored = load_network(archive)
    restored_group = next(
        item for item in restored.sys.groups if item.name == 'edta-retained-mg'
    )
    after_context = restored.get_solve_context_key()
    after_views = tuple(make_group_view(item, restored) for item in selections)
    after_results = attach_group_view_results(
        SimpleNamespace(), restored, selections,
        system='EDTA', radius_configuration_id='edta-retained-radii-v1',
        source_revision=restored.archive_id,
    )

    assert restored.archive_id == archive_id
    assert restored.source_network_id == producer_source_id
    assert after_context == before_context
    assert after_context.boundary_generator_digest == before_context.boundary_generator_digest
    assert after_context.environment_digest == before_context.environment_digest
    assert restored.solvent_provenance == before_solvent
    assert restored.retained_generator_provenance == before_retention
    assert restored.generator_system_indices == network.generator_system_indices
    assert restored.generator_stable_ids == network.generator_stable_ids
    assert restored.solvent_cells_complete == network.solvent_cells_complete
    assert restored_group.layer_atoms == group.layer_atoms
    assert restored_group.layer_net_atoms == group.layer_net_atoms
    assert all(view.completeness.proven for view in after_views)
    assert tuple(view.selection_id for view in after_views) == tuple(
        view.selection_id for view in before_views
    )
    for before, after in zip(before_results, after_results):
        assert after.status is Status.CERTIFIED
        assert before.identity.source_network_id == producer_source_id
        assert before.identity.source_network_id != archive_id
        assert after.identity.source_network_id == before.identity.source_network_id
        assert after.identity.solve_universe_id == before.identity.solve_universe_id
        assert after.identity.selection_id == before.identity.selection_id
        assert before.provenance.details['dependency_digest'] == after.provenance.details['dependency_digest']
        assert after.provenance.details['source_retained_generator_provenance'] == (
            restored.retained_generator_provenance
        )
        assert after.metrics['area'].value == pytest.approx(
            before.metrics['area'].value, abs=1e-12
        )
        assert after.metrics['volume'].value == pytest.approx(
            before.metrics['volume'].value, abs=1e-12
        )
    assert len(builds) == 1


def test_archive_writes_versioned_timing_sidecar(archive):
    sidecar = archive.with_name(archive.name + '.timing.json')

    manifest = json.loads(sidecar.read_text(encoding='utf-8'))

    assert manifest['format'] == 'vorpy.archive_timing'
    assert manifest['manifest_version'] == 1
    assert manifest['archive'] == archive.name
    timing = manifest['timing']
    assert timing['total'] >= 0.0
    assert timing['member_count'] > 0
    assert timing['archive_bytes'] == archive.stat().st_size
    assert timing['array_count'] >= 0
    assert timing['table_count'] >= 0
    assert timing['uncompressed_member_bytes'] >= timing['compressed_member_bytes']


def test_round_trip_all_scientific_columns(solved, archive, monkeypatch):
    from vorpy.src.network import Network
    from vorpy.src.group import Group
    def forbidden(*a, **kw):
        pytest.fail('Load attempted to build or analyze')
    for name in ('build', 'build_edges', 'build_surfaces', 'find_verts'):
        monkeypatch.setattr(Network, name, forbidden)
    monkeypatch.setattr(Group, 'get_info', forbidden)
    loaded = load_network(archive)
    original = solved.export_group.net
    for key in ('balls', 'surfs', 'edges', 'verts'):
        left, right = getattr(original, key), getattr(loaded, key)
        assert list(left.columns) == list(right.columns)
        assert left.index.tolist() == right.index.tolist()
        for column in left.columns:
            same(left[column].tolist(), right[column].tolist())
    same(original.settings, loaded.settings)
    same(original.box, loaded.box)
    same(original.sys.balls['loc'].tolist(), loaded.sys.balls['loc'].tolist())
    same(original.sys.balls['rad'].tolist(), loaded.sys.balls['rad'].tolist())
    assert loaded.sys.groups[0].net is loaded
    assert loaded.sys.balls.iloc[0]['res'] is loaded.sys.residues[0]
    assert loaded.sys.residues[0].sys is loaded.sys
    assert loaded.sys.balls.iloc[0]['chn'] in loaded.sys.chains
    assert loaded.sys.files['base_file'] is None


def test_export_equivalence_and_no_source(solved, archive, tmp_path):
    from vorpy.src.output.surfs import write_surfs
    from vorpy.src.output.edges import write_edges
    from vorpy.src.output.verts import write_off_verts
    from vorpy.src.output.pdb import write_pdb
    original = solved.export_group.net
    loaded = load_network(archive)
    old_source = solved.source
    hidden_source = old_source.with_suffix('.hidden')
    old_source.rename(hidden_source)
    cwd = Path.cwd()
    try:
        # Export all supported mesh formats from the same scientific values.
        for fmt in ('off', 'ply', 'vtp'):
            dirs = []
            for index, network in enumerate((original, loaded)):
                directory = tmp_path / f'{fmt}-{index}'
                directory.mkdir()
                dirs.append(directory)
                for name, function in (('surfs', write_surfs), ('edges', write_edges), ('verts', write_off_verts)):
                    function(network, [0, 1], name, directory=str(directory), file_type=fmt)
            for file in dirs[0].iterdir():
                assert file.read_bytes() == (dirs[1] / file.name).read_bytes()
        group = loaded.sys.groups[0]
        loaded.sys.files['dir'] = str(tmp_path)
        group.dir = str(tmp_path / 'group-export')
        Path(group.dir).mkdir()
        group.exports(atoms=True, info=True, logs=True, shell_surfs=True, shell_edges=True,
                      shell_verts=True, surr_atoms=True, surr_resids=True, atom_surfs=True,
                      atom_edges=True, atom_verts=True, sep_verts=True, file_type='ply')
        assert (Path(group.dir) / 'group_atoms.pdb').stat().st_size > 0
        assert (Path(group.dir) / 'info.txt').is_file()
        assert list(Path(group.dir).glob('*logs.csv'))
        assert list((Path(group.dir) / 'atoms').glob('*.ply'))
        assert list((Path(group.dir) / 'verts').glob('*.ply'))
        # A new complete-cell selection reuses existing geometry.
        selection = group_from_network(loaded, [0, 1], 'new_selection')
        assert selection.net.surfs is loaded.surfs
        selection.exports(atoms=True, shell_surfs=True, file_type='off')
        assert (Path(selection.dir) / 'shell_surfs.off').is_file()
    finally:
        os.chdir(cwd)
        hidden_source.rename(old_source)


def rewrite(path, member, transform):
    with zipfile.ZipFile(path) as z:
        files = {name: z.read(name) for name in z.namelist()}
    files[member] = json.dumps(transform(json.loads(files[member]))).encode()
    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as z:
        for name, data in files.items():
            z.writestr(name, data)


@pytest.mark.parametrize('field,value,match', [('format','wrong','Not a VorPy'),
                                               ('format_version',999,'Unsupported'),
                                               ('root_network',999,'root network')])
def test_invalid_metadata(archive, field, value, match):
    rewrite(archive, 'metadata.json', lambda data: dict(data, **{field:value}))
    with pytest.raises(ArchiveError, match=match):
        load_network(archive)


def test_count_mismatch(archive):
    rewrite(archive, 'metadata.json', lambda data: dict(data, counts=dict(data['counts'], verts=99)))
    with pytest.raises(ArchiveError, match='count mismatch'):
        load_network(archive)


def test_dangling_object_reference(archive):
    def corrupt(state):
        state['objects'][0]['fields']['sys'] = {'ref': 999999}
        return state
    rewrite(archive, 'state.json', corrupt)
    with pytest.raises(ArchiveError, match='Dangling object'):
        load_network(archive)


def test_missing_array(archive):
    with zipfile.ZipFile(archive) as z:
        files = {name:z.read(name) for name in z.namelist()}
    del files[next(name for name in files if name.startswith('arrays/'))]
    with zipfile.ZipFile(archive,'w') as z:
        for name,data in files.items():
            z.writestr(name,data)
    with pytest.raises(ArchiveError, match='array'):
        load_network(archive)


def test_process_exit_load_export(archive, tmp_path):
    import sys
    project_root = Path(__file__).resolve().parents[3]
    environment = os.environ.copy()
    environment['PYTHONPATH'] = os.pathsep.join(
        value for value in (str(project_root), environment.get('PYTHONPATH')) if value
    )
    code = '''from vorpy.src.io import load_network
from pathlib import Path
import sys
n=load_network(sys.argv[1]); g=n.sys.groups[0]
g.exports(atoms=True,shell_surfs=True,verts=True,file_type='off')
assert (Path(g.dir)/'group_atoms.pdb').is_file()
'''
    completed = subprocess.run(
        [sys.executable, '-c', code, str(archive)],
        cwd=project_root,
        env=environment,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert completed.returncode == 0, completed.stderr


@pytest.mark.timeout(120)
def test_aw_curvature_and_interface_records_round_trip(tmp_path):
    from vorpy.src.interface import Interface
    source = Path(__file__).resolve().parents[2] / 'data' / 'EDTA.pdb'
    cwd = Path.cwd()
    try:
        result = VorPyBackend().solve(source, lambda *a:None, lambda:False)
        group = result.export_group
        group.get_info()
        # Definitions and contact relationships on the existing solved topology.
        complete_system_ids = [
            int(system_id)
            for system_id, complete in zip(
                group.net.balls['system_num'], group.net.balls['complete']
            )
            if int(system_id) >= 0 and bool(complete)
        ]
        assert len(complete_system_ids) >= 2
        left_ids = {complete_system_ids[0]}
        right_ids = {complete_system_ids[1]}
        left = group_from_network(group.net, left_ids, 'left')
        right = group_from_network(group.net, right_ids, 'right')
        interface = Interface(group.sys, left, right)
        interface.net = group.net
        interface.group1_indices = set(left_ids)
        interface.group2_indices = set(right_ids)
        interface.water_topology = {'waters': {}, 'test_reference': group.sys.residues[0]}
        group.sys.ifaces.append(interface)
        archive = save_network(group.net, tmp_path / 'aw.vpy')
        loaded = load_network(archive)
        for table_name in ('balls','surfs','edges','verts'):
            original = getattr(group.net, table_name)
            restored = getattr(loaded, table_name)
            for column in original.columns:
                if 'curv' in column or 'energy' in column or 'angle' in column:
                    same(original[column].tolist(), restored[column].tolist())
        loaded_interface = loaded.sys.ifaces[0]
        assert loaded_interface.net is loaded
        assert loaded_interface.group1.interfaces[interface.interface_id] is loaded_interface
        assert loaded_interface.water_topology['test_reference'] is loaded.sys.residues[0]
        same(group.vol, loaded.sys.groups[0].vol)
        same(group.int_mean_curv, loaded.sys.groups[0].int_mean_curv)
        loaded.sys.set_output_directory(str(tmp_path / 'interface_output'))
        loaded_interface.export(atoms=True, surfs=True, edges=True, verts=True, info=True)
        assert (Path(loaded_interface.dir) / 'interface_atoms.pdb').is_file()
    finally:
        os.chdir(cwd)


def test_cli_archive_input_without_build(archive, monkeypatch, tmp_path):
    import sys
    from vorpy.src.command.vpy_cmnd import Command
    from vorpy.src.network import Network
    monkeypatch.setattr(Network, 'build', lambda *a, **kw: pytest.fail('CLI tried to rebuild'))
    monkeypatch.setattr(sys, 'argv', ['vorpy', '--load-network', str(archive), '-e', 'only', 'logs',
                                    '--save-network', str(tmp_path / 'resaved.vpy')])
    cwd = Path.cwd()
    try:
        command = Command()
        command.run()
        assert (tmp_path / 'resaved.vpy').is_file()
        assert list(Path(command.sys.files['dir']).rglob('*logs.csv'))
    finally:
        os.chdir(cwd)


def test_cached_alpha_dual_archive_round_trip_writes_logs_without_rebuild(
    solved, tmp_path, monkeypatch
):
    import csv
    import json
    from vorpy.src.interface import Interface
    from vorpy.src.interface.geometry_analysis import InterfaceGeometryAnalysis
    from vorpy.src.geometry.duals import build_dual
    from vorpy.src.io.scientific_adapters import stable_id
    from vorpy.src.network import Network
    from vorpy.src.output import visualization
    from vorpy.src.results import Status

    group = solved.export_group
    complete_system_ids = [
        int(system_id)
        for system_id, complete in zip(
            group.net.balls['system_num'], group.net.balls['complete']
        )
        if int(system_id) >= 0 and bool(complete)
    ]
    assert len(complete_system_ids) >= 2
    left_ids = {complete_system_ids[0]}
    right_ids = {complete_system_ids[1]}
    left = group_from_network(group.net, left_ids, 'cached_left')
    right = group_from_network(group.net, right_ids, 'cached_right')
    interface = Interface(
        group.sys, left, right, name='cached_alpha_dual_regression',
        interface_id='cached_alpha_dual_regression',
    )
    interface.net = group.net
    interface.group1_indices = set(left_ids)
    interface.group2_indices = set(right_ids)
    group.sys.ifaces.append(interface)
    analysis = InterfaceGeometryAnalysis(
        **json.loads(
            (Path(__file__).parent.parent / 'results' / 'fixtures' / '2kai_aw_cached_analysis.json')
            .read_text(encoding='utf-8')
        )
    )
    analysis.metadata.update(
        scheme='pow', alpha_convention='power_distance', alpha_units='A^2'
    )
    analysis.archive_provenance = {
        'fixture': '2kai_aw_cached_analysis',
        'network_id': stable_id(group.net),
    }
    interface.geometry_analysis = analysis
    interface._geometry_dual = build_dual(group.net)
    interface._geometry_dual.record_explicit_group_membership(
        interface.group1_indices, interface.group2_indices
    )

    archive = save_network(group.net, tmp_path / 'cached.vpy')
    with zipfile.ZipFile(archive) as bundle:
        state = json.loads(bundle.read('state.json'))
    kinds = {record['kind'] for record in state['objects']}
    assert {'interface_analysis', 'dual_complex', 'dual_simplex'} <= kinds
    assert not {'interface_result', 'dual_result', 'alpha_selection_result'} & kinds
    monkeypatch.setattr(Network, 'build', lambda *a, **kw: pytest.fail('Archive loader rebuilt topology'))
    monkeypatch.setattr(Interface, 'analyze_geometry', lambda *a, **kw: pytest.fail('Archive log rebuilt geometry'))

    loaded = load_network(archive)
    restored = next(
        item for item in loaded.sys.ifaces
        if item.interface_id == interface.interface_id
    )
    assert restored.geometry_analysis.alpha_selection == analysis.alpha_selection
    assert restored.geometry_analysis.archive_provenance['fixture'] == '2kai_aw_cached_analysis'
    assert restored._geometry_dual.network is loaded
    assert restored.group1_indices == interface.group1_indices
    assert restored.group2_indices == interface.group2_indices
    assert restored._geometry_dual.metadata['explicit_group_membership'] == (
        interface._geometry_dual.metadata['explicit_group_membership']
    )
    alpha_result = visualization._cached_alpha_selection_result(restored)
    assert alpha_result.status is Status.PARTIAL
    assert alpha_result.provenance.details['fixture'] == '2kai_aw_cached_analysis'

    first = tmp_path / 'restored-1' / 'logs.csv'
    second = tmp_path / 'restored-2' / 'logs.csv'
    visualization._compact_log(restored, first)
    visualization._compact_log(restored, second)
    assert first.read_bytes() == second.read_bytes()
    with first.open(newline='', encoding='utf-8') as stream:
        rows = list(csv.DictReader(stream))
    dual_rows = [row for row in rows if row['representation'].startswith('dual_')]
    assert {
        (row['representation'], row['side']) for row in dual_rows
    } >= {
        ('dual_contact_complex', 'shared'),
        ('dual_generator_complex', 'A'),
        ('dual_generator_complex', 'B'),
        ('dual_generator_complex', 'AB/contact'),
    }
    for row in dual_rows:
        details = json.loads(row['provenance'])['details']
        assert row['interface_id'] == restored.archive_id
        assert details['network_id'] == loaded.archive_id
        assert details['partition'] == 'power'
