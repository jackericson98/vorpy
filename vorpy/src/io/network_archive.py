"""Versioned, self-contained solved-network archives (ZIP + JSON + NumPy).

Public API accepts a Network, Group, or System on save and returns the selected
Network on load. Its .sys retains all archived groups/interfaces and networks.
No solve, triangulation, or scientific analysis is invoked by the loader.
"""
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import json
import os
import tempfile
import zipfile
import hashlib
from time import perf_counter

import numpy as np
import pandas as pd

from vorpy.src.io.archive_codec import ArchiveError, Reader, Writer
from vorpy.src.io.archive_fields import (
    GROUP_FIELDS, INTERFACE_FIELDS, NETWORK_FIELDS, SYSTEM_FIELDS, RESIDUE_FIELDS, CHAIN_FIELDS,
)
from vorpy.src.version import __version__
from .session import Session
from .scientific_adapters import (
    registry, stable_id, scope, canonical_environment, canonical_network_scope,
    normalize_vocabulary, restore_runtime_keys,
)


def _process_memory():
    """Return a process memory reading without making psutil mandatory."""
    try:
        import resource
    except ImportError:
        resource = None
    if resource is not None:
        value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return int(value * (1024 if os.name != 'darwin' else 1)), 'resource_ru_maxrss'
    try:
        import psutil
    except ImportError:
        return None, None
    return int(psutil.Process().memory_info().rss), 'psutil_rss_current'


def _progress(progress, stage, completed, total=None):
    if progress is not None:
        progress(stage, completed, total)


def _timing_manifest_path(path):
    path = Path(path)
    return path.with_name(path.name + '.timing.json')


def _write_timing_manifest(path, profile):
    """Atomically persist diagnostics beside the archive, not in its schema."""
    destination = _timing_manifest_path(path)
    manifest = {
        'format': 'vorpy.archive_timing',
        'manifest_version': 1,
        'archive': Path(path).name,
        'timing': profile,
    }
    fd, temp_name = tempfile.mkstemp(prefix='.vorpy-timing-', suffix='.tmp', dir=destination.parent)
    try:
        with os.fdopen(fd, 'w', encoding='utf-8', newline='') as stream:
            json.dump(manifest, stream, allow_nan=False, separators=(',', ':'))
            stream.write('\n')
        os.replace(temp_name, destination)
    finally:
        if os.path.exists(temp_name):
            os.unlink(temp_name)
    return destination

FORMAT = 'VorPy Network Archive'
FORMAT_VERSION = 3
FIELDS = {'system': SYSTEM_FIELDS, 'network': NETWORK_FIELDS, 'group': GROUP_FIELDS,
          'interface': INTERFACE_FIELDS, 'residue': RESIDUE_FIELDS,
          'chain': CHAIN_FIELDS, 'solvent': CHAIN_FIELDS}
V1_FIELDS = dict(FIELDS)
LINKS = {
    'system': ('groups', 'ifaces', 'residues', 'chains', 'sol'),
    'network': ('sys',),
    'group': ('sys', 'net', 'interfaces', 'parent_interface'),
    'interface': ('sys', 'net', 'group1', 'group2', 'partial_group1', 'partial_group2', 'water_groups', 'buried_water_groups'),
    'residue': ('sys', 'chain', 'mol'),
    'chain': ('residues',), 'solvent': ('residues',),
}
ADAPTERS = registry()
OPTIONAL_ADAPTER_FIELDS = {
    # These fields were added after the first molecular-contact archives.
    # Missing provenance remains unavailable; it is never inferred on load.
    'molecular_contact_surface': frozenset({'contact_selection_source', 'curvature'}),
}
FIELDS.update({kind: names + ('archive_id', 'archive_provenance') for kind, (_, names) in ADAPTERS.items()})
FIELDS['session'] = ('archive_id', 'active_network_id', 'workbench')
LINKS['session'] = ('systems', 'networks')
for kind in ADAPTERS:
    LINKS[kind] = ()
for kind in ('system', 'network', 'group', 'interface'):
    FIELDS[kind] += ('archive_id',)
FIELDS['network'] += ('network_scope', 'generator_provenance', 'boundary_indices', 'boundary_config', 'boundary_mode')
FIELDS['system'] += ('boundary_config', 'boundary_generators', 'boundary_mode',
    'interface_geometry_cache', '_interface_phase1_dry_cache')
FIELDS['interface'] += ('geometry_analysis', '_geometry_dual', '_geometry_filtration',
    '_geometry_representation', 'analysis_provenance', 'archive_cache_state')
V2_FIELDS = dict(FIELDS)
FIELDS['interface'] += ('representation_caches', 'representation_links')

_REPRESENTATION_FIELDS = frozenset({
    'geometry_analysis', '_geometry_dual', '_geometry_filtration',
    '_geometry_representation', 'representation_caches', 'representation_links',
    'analysis_provenance', 'archive_cache_state',
})
_OBJECT_FIELDS = {kind: FIELDS[kind] + LINKS[kind] for kind in FIELDS}


def _water_topology_snapshot(value):
    # Discovery's local lookup closure is executable runtime state, not
    # scientific data. Preserve the lookup maps and all water results.
    if isinstance(value, dict):
        value = dict(value)
        context = value.get('context')
        if isinstance(context, dict):
            value['context'] = {key: item for key, item in context.items()
                                if key != 'topology_to_parent'}
    return value


def _archive_field_value(kind, name, value):
    if kind == 'interface' and name == 'water_topology':
        return _water_topology_snapshot(value)
    if kind == 'interface' and name == 'water_geometries' and value is not None:
        return [dict(geometry, topology_analysis=_water_topology_snapshot(geometry['topology_analysis']))
                if isinstance(geometry, dict) and 'topology_analysis' in geometry else geometry
                for geometry in value]
    return value


def _archive_field(obj, kind, name):
    return _archive_field_value(kind, name, getattr(obj, name))


def _classes():
    from vorpy.src.system import System
    from vorpy.src.network import Network
    from vorpy.src.group import Group
    from vorpy.src.interface import Interface
    from vorpy.src.objects.residue import Residue
    from vorpy.src.objects.chain import Chain, Sol
    return {'system': System, 'network': Network, 'group': Group, 'interface': Interface,
            'residue': Residue, 'chain': Chain, 'solvent': Sol, 'session': Session,
            **{kind: cls for kind, (cls, _) in ADAPTERS.items()}}


def _root(value):
    classes = _classes()
    if isinstance(value, classes['network']):
        return value
    if isinstance(value, classes['group']):
        return value.net
    if isinstance(value, classes['system']):
        return next((g.net for g in (value.groups or []) + (value.ifaces or []) if g.net is not None), None)
    raise ArchiveError('Expected a VorPy Network, Group, or System')


def _collect(root):
    types = {cls: kind for kind, cls in _classes().items()}
    objects, ids = [], {}

    def visit(value):
        if type(value) in types:
            value_id = id(value)
            if value_id in ids:
                return
            ids[value_id] = len(objects)
            objects.append(value)
            kind = types[type(value)]
            for name in _OBJECT_FIELDS[kind]:
                visit(getattr(value, name, None))
        elif isinstance(value, dict):
            for item in value.values():
                visit(item)
        elif isinstance(value, (list, tuple, set, frozenset)):
            for item in value:
                visit(item)
    visit(root)
    return objects, ids, types


def save_network(network, path, *, progress=None):
    """Atomically save the root network and its molecular/group/interface context."""
    root = _root(network)
    if root is None or root.sys is None:
        raise ArchiveError('A solved Network with a parent System is required')
    validate_network(root)
    return _save_document(root, path, root, 'network', progress=progress)


def save_session(session, path, *, progress=None):
    """Save a complete session in the same .vpy object graph format."""
    if not isinstance(session, Session):
        session = Session.from_network(_root(session))
    # Normalize document membership using existing references only. This also
    # supports a session built from systems without a separately supplied list
    # of their reachable networks.
    for network in session.networks:
        if all(network.sys is not system for system in session.systems):
            session.systems.append(network.sys)
    for system in session.systems:
        for owner in (getattr(system, 'groups', None) or []) + (getattr(system, 'ifaces', None) or []):
            network = getattr(owner, 'net', None)
            if network is not None and all(network is not item for item in session.networks):
                session.networks.append(network)
    for network in session.networks:
        if all(getattr(network, name, None) is not None for name in ('balls', 'verts', 'edges', 'surfs')):
            validate_network(network)
    return _save_document(session, path, session.active_network, 'session', progress=progress)


def _save_document(document, path, root, document_kind, *, progress=None):
    profile = {name: 0.0 for name in (
        'archive_preparation', 'object_collection', 'object_table_traversal',
        'scientific_object_preparation', 'scientific_object_conversion',
        'representation_cache_serialization', 'array_conversion',
        'array_serialization', 'array_serialization_and_write',
        'table_traversal', 'json_serialization', 'metadata_serialization',
        'archive_member_write', 'compression', 'zip_archive_write',
        'member_compression_write', 'zip_finalization', 'filesystem_replace', 'total',
    )}
    profile['object_type_timings'] = {}
    profile['table_diagnostics'] = {}
    profile['size_histograms'] = {'members': {}, 'arrays': {}}
    profile.update({
        'member_count': 0,
        'uncompressed_member_bytes': 0,
        'compressed_member_bytes': 0,
        'array_uncompressed_bytes': 0,
        'array_compressed_bytes': 0,
        'table_uncompressed_bytes': 0,
        'table_compressed_bytes': 0,
        'json_member_sizes': {},
    })
    profile['timing_semantics'] = {
        'archive_preparation': 'inclusive document normalization and metadata preparation',
        'object_collection': 'exclusive reachable-object collection',
        'scientific_object_preparation': 'exclusive stable identity and provenance preparation',
        'object_table_traversal': 'inclusive object field encoding and record construction',
        'scientific_object_conversion': 'exclusive field encoding for scientific objects',
        'representation_cache_serialization': 'exclusive representation/cache field encoding',
        'array_conversion': 'exclusive NumPy coercion before NPY serialization',
        'table_traversal': 'inclusive table column conversion and table member write',
        'json_serialization': 'exclusive JSON serialization and UTF-8 staging; same interval as json_staging',
        'json_staging': 'exclusive JSON serialization and UTF-8 staging',
        'metadata_serialization': 'exclusive metadata JSON serialization and UTF-8 staging',
        'array_serialization': 'exclusive NPY serialization for small arrays',
        'array_serialization_and_write': 'inclusive NPY serialization and member write for large arrays',
        'archive_member_write': 'inclusive ZIP compression and member write; same interval as member_compression_write',
        'compression': 'inclusive ZIP compression and member write; compatibility alias',
        'member_compression_write': 'inclusive ZIP compression and member write',
        'zip_archive_write': 'inclusive ZIP body and finalization',
        'zip_finalization': 'exclusive ZipFile.close central-directory finalization',
        'filesystem_replace': 'exclusive atomic archive rename',
        'total': 'inclusive archive preparation through atomic archive rename; excludes timing-manifest write',
    }
    document_started = perf_counter()
    started = document_started
    objects, ids, types = _collect(document)
    profile['object_collection'] = perf_counter() - started
    _progress(progress, 'object_collection', len(objects), len(objects))
    started = perf_counter()
    for obj in objects:
        kind = types[type(obj)]
        normalize_vocabulary(obj)
        if kind in {'system', 'network', 'group', 'interface', 'session'} or kind in ADAPTERS:
            stable_id(obj)
        if kind == 'network':
            obj.network_scope = canonical_network_scope(getattr(obj, 'network_scope', None) or scope(obj))
        if kind == 'interface':
            caches = getattr(obj, 'representation_caches', None) or {}
            links = getattr(obj, 'representation_links', None) or {}
            metadata = getattr(getattr(obj, 'geometry_analysis', None), 'metadata', {}) or {}
            environment = canonical_environment(
                metadata.get('environment', 'solvent_competing' if metadata.get('solvent_context') else 'dry'))
            network_id = stable_id(obj.net)
            interface_id = stable_id(obj)
            provenance = {
                'network_id': network_id,
                'system_id': stable_id(obj.sys),
                'interface_id': interface_id,
                'network_scope': canonical_network_scope(getattr(obj.net, 'network_scope', None) or scope(obj.net)),
                'partition': getattr(obj.net, 'settings', {}).get('net_type'),
                'environment': environment,
            }
            for representation, cache in caches.items():
                link = links.get(representation)
                if isinstance(link, dict):
                    link.setdefault('cache_id', stable_id(cache) if type(cache) in types else None)
                    link.setdefault('network_id', network_id)
                    link.setdefault('interface_id', interface_id)
                    link.setdefault('representation', representation)
                    link.setdefault('environment', environment)
                    link.setdefault('network_scope', provenance['network_scope'])
                if cache is not None and not isinstance(cache, dict):
                    cache_provenance = dict(getattr(cache, 'archive_provenance', {}) or {})
                    cache_provenance.setdefault('network_id', network_id)
                    cache_provenance.setdefault('system_id', provenance['system_id'])
                    cache_provenance.setdefault('interface_id', interface_id)
                    cache_provenance.setdefault('network_scope', provenance['network_scope'])
                    cache_provenance.setdefault('partition', provenance['partition'])
                    cache_provenance.setdefault('environment', environment)
                    cache_provenance.setdefault('representation', representation)
                    if representation.startswith('molecular_contact_surface'):
                        selection_source = getattr(cache, 'contact_selection_source', None)
                        if selection_source:
                            cache_provenance.setdefault('contact_selection_source', selection_source)
                    if cache_provenance != getattr(cache, 'archive_provenance', None):
                        object.__setattr__(cache, 'archive_provenance', cache_provenance)
            _validate_cache_owners(obj)
    stable_values = [getattr(obj, 'archive_id', None) for obj in objects if hasattr(obj, 'archive_id')]
    if len(stable_values) != len(set(stable_values)):
        raise ArchiveError('Duplicate stable archive identity')
    for obj in objects:
        if types[type(obj)] == 'interface':
            analysis = getattr(obj, 'geometry_analysis', None)
            key = getattr(obj, '_geometry_analysis_key', None)
            if key is not None and analysis is not None:
                from vorpy.src.interface.geometry_analysis import geometry_cache_key
                curvature_requested = bool(analysis.metadata.get('curvature_requested', True))
                expected_key = (*geometry_cache_key(obj, analysis.metadata.get('alpha_value', 0.)), curvature_requested)
                obj.archive_cache_state = ('AVAILABLE' if key == expected_key
                                           else 'STALE')
            for name in ('geometry_analysis', '_geometry_dual', '_geometry_filtration', '_geometry_representation'):
                cache = getattr(obj, name, None)
                if cache is not None:
                    metadata = getattr(analysis, 'metadata', {})
                    provenance = {'network_id': stable_id(obj.net), 'system_id': stable_id(obj.sys),
                        'network_scope': canonical_network_scope(obj.net.network_scope),
                        'partition': obj.net.settings.get('net_type')}
                    if name in ('geometry_analysis', '_geometry_representation'):
                        filtration = getattr(obj, '_geometry_filtration', None)
                        alpha = metadata.get('alpha_value')
                        provenance.update(interface_id=stable_id(obj), alpha=alpha,
                            environment=canonical_environment(
                                metadata.get('environment',
                                             'solvent_competing' if metadata.get('solvent_context') else 'dry')),
                            representation='Voronoi', alpha_convention=metadata.get('alpha_convention'),
                            alpha_units=metadata.get('alpha_units'),
                            query_id=hashlib.sha256(f'{stable_id(filtration) if filtration else obj.net.archive_id}:{alpha!r}'.encode()).hexdigest())
                    if not hasattr(cache, 'archive_provenance'):
                        object.__setattr__(cache, 'archive_provenance', provenance)
        elif types[type(obj)] == 'network' and not hasattr(obj, 'generator_provenance'):
            boundary = set(getattr(obj, 'boundary_indices', ()))
            obj.generator_provenance = {
                int(index): {'generator_id': int(index), 'network_id': obj.archive_id,
                    'system_id': stable_id(obj.sys),
                    'system_atom_index': (None if int(index) in boundary else int(row.get('system_num', row['num']))),
                    'stable_generator_id': f'{obj.archive_id}:generator:{int(index)}',
                    'source_atom_id': row.get('stable_id'),
                    'is_boundary_generator': int(index) in boundary}
                for index, row in obj.balls.iterrows()}
    profile['scientific_object_preparation'] = perf_counter() - started
    destination = Path(path).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    source = root.sys.files.get('base_file') if root is not None else None
    metadata = {
        'format': FORMAT, 'format_version': FORMAT_VERSION, 'vorpy_version': __version__,
        'created_at': datetime.now(timezone.utc).isoformat(),
        'document_kind': document_kind, 'document_root': ids[id(document)],
        'root_network': ids[id(root)] if root is not None else None,
        'name': root.group_name if root is not None else 'Workbench session', 'original_structure_path': source,
        'original_structure_filename': Path(source).name if source else None,
        'units': {'length': 'angstrom', 'area': 'angstrom^2', 'volume': 'angstrom^3',
                  'mean_curvature': 'angstrom^-1', 'gaussian_curvature': 'angstrom^-2',
                  'integrated_mean_curvature': 'angstrom',
                  'integrated_gaussian_curvature': '1'},
        'counts': ({'atoms': len(root.sys.balls), **{key: len(getattr(root, key))
            if getattr(root, key, None) is not None else None for key in ('balls', 'surfs', 'edges', 'verts')}}
            if root is not None else {}),
        'solve_duration_seconds': root.metrics.get('tot') if root is not None else None,
        **{key: root.settings.get(setting) if root is not None else None for key, setting in (
            ('scheme', 'net_type'), ('surface_resolution', 'surf_res'), ('boundary_padding', 'box_size'),
            ('maximum_vertex_radius', 'max_vert'), ('probe_size', 'probe_size'))},
    }
    profile['archive_preparation'] = perf_counter() - document_started
    _progress(progress, 'archive_preparation', 1, 1)
    fd, temp_name = tempfile.mkstemp(prefix='.vorpy-', suffix='.tmp', dir=destination.parent)
    os.close(fd)
    writer = None
    try:
        started = perf_counter()
        with zipfile.ZipFile(temp_name, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=4, allowZip64=True) as archive:
            writer = Writer(archive, ids, profile, progress=progress)
            records = []
            traversal_started = perf_counter()
            for object_index, obj in enumerate(objects, 1):
                kind = types[type(obj)]
                type_profile = profile['object_type_timings'].setdefault(
                    kind, {'count': 0, 'seconds': 0.0, 'fields': {}})
                type_profile['count'] += 1
                object_started = perf_counter()
                fields = {}
                for name in _OBJECT_FIELDS[kind]:
                    try:
                        value = getattr(obj, name)
                    except AttributeError:
                        continue
                    field_started = perf_counter()
                    fields[name] = writer.encode(_archive_field_value(kind, name, value))
                    field_elapsed = perf_counter() - field_started
                    timing_name = (
                        'representation_cache_serialization'
                        if kind in ADAPTERS or name in _REPRESENTATION_FIELDS
                        else 'scientific_object_conversion'
                    )
                    profile[timing_name] = profile.get(timing_name, 0.0) + field_elapsed
                    type_profile['fields'][name] = type_profile['fields'].get(name, 0.0) + field_elapsed
                records.append({'id': ids[id(obj)], 'kind': kind, 'schema_version': 1, 'fields': fields})
                type_profile['seconds'] += perf_counter() - object_started
                if progress is not None and object_index % 4096 == 0:
                    progress('object_table_traversal', object_index, len(objects))
            profile['object_table_traversal'] = perf_counter() - traversal_started
            _progress(progress, 'object_table_traversal', len(objects), len(objects))
            metadata['solve_settings'] = writer.encode(root.settings if root is not None else {})
            writer.json('state.json', {'objects': records})
            writer.json('metadata.json', metadata)
            zip_body_finished = perf_counter()
        profile['zip_finalization'] = perf_counter() - zip_body_finished
        profile['zip_archive_write'] = perf_counter() - started
        _progress(progress, 'zip_finalization', 1, 1)
        started = perf_counter()
        os.replace(temp_name, destination)
        profile['filesystem_replace'] = perf_counter() - started
    finally:
        if os.path.exists(temp_name):
            os.unlink(temp_name)
    profile['total'] = perf_counter() - document_started
    profile['archive_bytes'] = destination.stat().st_size
    profile['archive_uncompressed_bytes'] = profile['uncompressed_member_bytes']
    profile['archive_compressed_bytes'] = profile['compressed_member_bytes']
    profile['object_count'] = len(objects)
    profile['array_count'] = writer.array_count if writer is not None else 0
    profile['table_count'] = len(writer.tables) if writer is not None else 0
    profile['member_count'] = writer.member_count if writer is not None else 0
    profile['peak_process_memory_bytes'], peak_source = _process_memory()
    profile['memory_source'] = peak_source
    profile['timing_manifest'] = str(_timing_manifest_path(destination))
    _write_timing_manifest(destination, profile)
    if root is not None and getattr(root, 'sys', None) is not None:
        root.sys.archive_timing = profile
    return destination


@dataclass(frozen=True)
class _Reference:
    id: int


def _resolve(value, objects):
    if isinstance(value, _Reference):
        return objects[value.id]
    if isinstance(value, pd.DataFrame):
        for column in value.columns:
            if value[column].dtype == object and any(isinstance(v, _Reference) for v in value[column]):
                value[column] = [objects[v.id] if isinstance(v, _Reference) else v for v in value[column]]
        return value
    if isinstance(value, dict):
        return {key: _resolve(item, objects) for key, item in value.items()}
    if isinstance(value, list):
        return [_resolve(item, objects) for item in value]
    if isinstance(value, tuple):
        return tuple(_resolve(item, objects) for item in value)
    if isinstance(value, (set, frozenset)):
        return type(value)(_resolve(item, objects) for item in value)
    return value


def _read_v1(reader, metadata):
    records = reader.json('state.json')['objects']
    if not isinstance(records, list) or not records:
        raise ArchiveError('Missing object records')
    for index, record in enumerate(records):
        kind = record['kind']
        if record['id'] != index or kind not in FIELDS or not isinstance(record['fields'], dict):
            raise ArchiveError('Invalid object identity/type')
        record_version = record.get('schema_version', 1 if metadata['format_version'] == 1 else None)
        if type(record_version) is not int or record_version != 1:
            raise ArchiveError(f'Unsupported {kind} schema version')
        if metadata['format_version'] == 1 and (kind in ADAPTERS or kind == 'session'):
            raise ArchiveError('Scientific session records require archive version 2')
        required_fields = (set(ADAPTERS[kind][1])
                           - set(OPTIONAL_ADAPTER_FIELDS.get(kind, ()))) if kind in ADAPTERS else set()
        if required_fields - set(record['fields']):
            raise ArchiveError(f'Missing required cached {kind} fields')
        allowed = (V1_FIELDS[kind] if metadata['format_version'] == 1
                   else V2_FIELDS[kind] if metadata['format_version'] == 2
                   else FIELDS[kind])
        if set(record['fields']) - set(allowed + LINKS[kind]):
            raise ArchiveError(f'Unknown fields for {kind}')
        reader.references[index] = _Reference(index)
    # Parse/validate all arrays and references before allocating core objects.
    decoded = [{key: reader.decode(value) for key, value in record['fields'].items()} for record in records]
    classes = _classes()
    objects = [classes[record['kind']].__new__(classes[record['kind']]) for record in records]
    for record, fields, obj in zip(records, decoded, objects):
        for key, value in fields.items():
            object.__setattr__(obj, key, _resolve(value, objects))
        kind = record['kind']
        if kind == 'molecular_contact_surface' and not hasattr(obj, 'curvature'):
            # Curvature was added after the original contact-surface archive
            # schema; absence means genuinely unavailable, never zero.
            obj.curvature = None
        if kind == 'molecular_contact_surface' and not hasattr(obj, 'contact_selection_source'):
            # Selection provenance was added after the original contact-surface
            # archive schema; absence remains a legacy limitation.
            obj.contact_selection_source = None
        if kind == 'system':
            obj.files = dict.fromkeys(('base_file', 'ball_file', 'verts_file', 'net_file', 'ndx_file', 'frame_files'))
            obj.files.update(dir=None, root_dir=str(Path.cwd()))
            obj.gui = None
            obj.atoms = None
            obj.user_atoms = None
            obj.simple = False
            obj.print_actions = False
            obj.verbose = False
            if not hasattr(obj, 'interface_geometry_cache'):
                obj.interface_geometry_cache = {'verts': {}, 'edges': {}, 'surfs': {}}
            obj.spatial_index_cache = None
            obj._spatial_revision = 0
            obj.run_start = obj.run_end = None
            obj.run_active = False
            obj.run_process = ''
            obj.run_progress = 0.0
            obj.run_network = None
            obj._progress_last_emit = 0.0
            obj._progress_last_process = None
            obj._progress_line_len = 0
            obj.cmnds = None
        elif kind == 'network':
            obj.progress_window = None
            obj.progress_network_name = None
            obj.progress_process_prefix = None
            obj.loaded_from_archive = True
        elif kind in {'group', 'interface'}:
            obj.dir = None
            if kind == 'group':
                obj.verts = None
            else:
                obj._buried_water_exports_complete = False
    for record, obj in zip(records, objects):
        kind = record['kind']
        if kind in {'system', 'network', 'group', 'interface', 'session'} or kind in ADAPTERS:
            stable_id(obj)
        if kind == 'network':
            obj.network_scope = canonical_network_scope(getattr(obj, 'network_scope', None) or scope(obj))
        normalize_vocabulary(obj)
        if kind in {'dual_complex', 'alpha_filtration', 'aw_alpha_filtration', 'apollonius_complex'}:
            if not isinstance(obj.network, classes['network']):
                raise ArchiveError(f'Invalid {kind} network reference')
    for record, obj in zip(records, objects):
        kind = record['kind']
        if kind == 'interface':
            _validate_cache_owners(obj)
            restore_runtime_keys(obj)
        elif kind == 'boundary_config':
            from vorpy.src.boundary import BoundaryMode
            object.__setattr__(obj, 'mode', BoundaryMode(obj.mode) if obj.mode is not None else None)
    identity_values = [getattr(obj, 'archive_id', None) for obj in objects if hasattr(obj, 'archive_id')]
    if any(not isinstance(value, str) or not value for value in identity_values):
        raise ArchiveError('Invalid stable archive identity')
    if len(identity_values) != len(set(identity_values)):
        raise ArchiveError('Duplicate stable archive identity')
    document_kind = metadata.get('document_kind', 'network' if metadata['format_version'] == 1 else None)
    if document_kind not in {'network', 'session'}:
        raise ArchiveError('Unknown document kind')
    root_id = metadata['root_network']
    if root_id is None and document_kind == 'session':
        document = _session_root(metadata, records, objects)
        _validate_context(objects, records)
        return document
    if not isinstance(root_id, int) or not 0 <= root_id < len(objects) or records[root_id]['kind'] != 'network':
        raise ArchiveError('Invalid root network ID')
    root = objects[root_id]
    for record, obj in zip(records, objects):
        if record['kind'] == 'network' and all(getattr(obj, key, None) is not None for key in ('balls', 'verts', 'edges', 'surfs')):
            validate_network(obj)
    for key in ('balls', 'surfs', 'edges', 'verts'):
        if metadata['counts'][key] != (len(getattr(root, key)) if getattr(root, key, None) is not None else None):
            raise ArchiveError(f'Entity count mismatch: {key}')
    if metadata['counts']['atoms'] != len(root.sys.balls):
        raise ArchiveError('Atom count mismatch')
    _validate_context(objects, records)
    root.archive_metadata = metadata
    root.sys.archive_metadata = metadata
    if document_kind == 'session':
        document = _session_root(metadata, records, objects)
        if all(root is not item for item in document.networks):
            raise ArchiveError('Session root network is not a session member')
        if document.active_network is not root:
            raise ArchiveError('Session active network disagrees with archive root')
        return document
    return root


def _session_root(metadata, records, objects):
    index = metadata.get('document_root')
    if type(index) is not int or not 0 <= index < len(objects) or records[index]['kind'] != 'session':
        raise ArchiveError('Invalid session document root')
    session = objects[index]
    if any(not isinstance(system, _classes()['system']) for system in session.systems):
        raise ArchiveError('Invalid session system member')
    if any(not isinstance(net, _classes()['network']) for net in session.networks):
        raise ArchiveError('Invalid session network member')
    if any(all(net.sys is not system for system in session.systems) for net in session.networks):
        raise ArchiveError('Session network belongs to an absent system')
    if session.active_network_id is not None and not any(
            net.archive_id == session.active_network_id for net in session.networks):
        raise ArchiveError('Invalid active network identity')
    return session


def _validate_cache_owners(iface):
    for name in ('geometry_analysis', '_geometry_dual', '_geometry_filtration', '_geometry_representation'):
        cache = getattr(iface, name, None)
        provenance = getattr(cache, 'archive_provenance', {})
        if provenance and provenance.get('network_id') != getattr(iface.net, 'archive_id', None):
            raise ArchiveError(f'{name} provenance references another network')
    for name in ('_geometry_dual', '_geometry_filtration'):
        cache = getattr(iface, name, None)
        if cache is not None and (cache.network is not iface.net
                or cache.scheme != iface.net.settings.get('net_type', 'aw')):
            raise ArchiveError(f'{name} belongs to another network or partition')
    source = getattr(getattr(iface, '_geometry_filtration', None), '_aw_source', None)
    if source is not None and (source.network is not iface.net or source.incidence.network is not iface.net):
        raise ArchiveError('AW source belongs to another network')
    caches = getattr(iface, 'representation_caches', None)
    if caches is not None and not isinstance(caches, dict):
        raise ArchiveError('representation_caches must be a mapping')
    links = getattr(iface, 'representation_links', None)
    if links is not None and not isinstance(links, dict):
        raise ArchiveError('representation_links must be a mapping')
    for key, link in (links or {}).items():
        if not isinstance(key, str) or not isinstance(link, dict):
            raise ArchiveError('Invalid representation link')
        source_network = link.get('network_id', link.get('source_network_id'))
        if source_network is not None and source_network != getattr(iface.net, 'archive_id', None):
            raise ArchiveError('representation link references another network')
        source_interface = link.get('interface_id', link.get('source_interface_id'))
        if source_interface is not None and source_interface != getattr(iface, 'archive_id', None):
            raise ArchiveError('representation link references another interface')
    for key, cache in (caches or {}).items():
        if not isinstance(key, str):
            raise ArchiveError('Invalid representation cache key')
        provenance = cache.get('archive_provenance', {}) if isinstance(cache, dict) else getattr(cache, 'archive_provenance', {})
        if provenance and provenance.get('network_id') != getattr(iface.net, 'archive_id', None):
            raise ArchiveError('representation cache references another network')


# Explicit version dispatch is the insertion point for future schema migrations.
LOADERS = {1: _read_v1, 2: _read_v1, 3: _read_v1}


def load_network(path):
    """Load a self-contained solved network. Never build/analyze geometry."""
    document = _load_document(path)
    if isinstance(document, Session):
        if document.active_network is None:
            raise ArchiveError('Session has no network; use load_session')
        return document.active_network
    return document


def load_session(path):
    """Load a session, or wrap a legacy/network archive without analysis."""
    document = _load_document(path)
    return document if isinstance(document, Session) else Session.from_network(document)


def _load_document(path):
    try:
        with zipfile.ZipFile(path, 'r') as archive:
            names = archive.namelist()
            if len(names) != len(set(names)):
                raise ArchiveError('Duplicate archive members')
            if any(name.startswith('/') or '..' in Path(name).parts or '\\' in name for name in names):
                raise ArchiveError('Unsafe archive member name')
            if any(not (name in {'metadata.json', 'state.json'} or name.startswith(('arrays/', 'tables/'))) for name in names):
                raise ArchiveError('Unexpected archive member')
            if sum(info.file_size for info in archive.infolist()) > 128 * 1024**3:
                raise ArchiveError('Archive exceeds the 128 GiB uncompressed safety limit')
            reader = Reader(archive)
            metadata = reader.json('metadata.json')
            if metadata.get('format') != FORMAT:
                raise ArchiveError('Not a VorPy Network Archive')
            version = metadata.get('format_version')
            if type(version) is not int or version not in LOADERS:
                raise ArchiveError(
                    f'Unsupported VorPy archive format version: {version!r}; '
                    f'supported: {", ".join(map(str, sorted(LOADERS)))}')
            document = LOADERS[version](reader, metadata)
            document.archive_path = str(Path(path).resolve())
            networks = document.networks if isinstance(document, Session) else [document]
            for root in networks:
                root.archive_metadata = metadata
                root.archive_path = str(Path(path).resolve())
                root.sys.archive_metadata = metadata
                root.sys.files['dir'] = str(Path(path).resolve().with_suffix('')) + '_exports'
            return document
    except ArchiveError:
        raise
    except (zipfile.BadZipFile, KeyError, TypeError, ValueError, AttributeError, IndexError, EOFError, OSError) as error:
        raise ArchiveError(f'Invalid VorPy network archive: {error}') from error


def _ids(values, valid, label):
    if values is None:
        return
    if not isinstance(values, (list, tuple, set, np.ndarray)):
        raise ArchiveError(f'{label}: expected a list of IDs')
    for value in values:
        if not isinstance(value, (int, np.integer)) or int(value) not in valid:
            raise ArchiveError(f'{label}: invalid/dangling ID {value!r}')


def validate_network(net):
    """Validate table shapes, topology ID spaces, and mesh-local connectivity."""
    for key in ('balls', 'verts', 'edges', 'surfs'):
        table = getattr(net, key, None)
        if not isinstance(table, pd.DataFrame):
            raise ArchiveError(f'Completed network requires the {key} table')
        if not table.index.is_unique:
            raise ArchiveError(f'Duplicate {key} table IDs')
    required = {'balls': ('loc', 'rad', 'num'), 'verts': ('loc', 'rad', 'balls'),
                'edges': ('balls', 'verts', 'points'), 'surfs': ('balls', 'verts', 'edges', 'points', 'tris')}
    for key, columns in required.items():
        if not set(columns).issubset(getattr(net, key).columns):
            raise ArchiveError(f'Missing required {key} columns')
    ball_ids = set(int(v) for v in net.balls['num'])
    if len(ball_ids) != len(net.balls):
        raise ArchiveError('Duplicate cell topology IDs')
    valid = {'balls': ball_ids, **{key: set(range(len(getattr(net, key)))) for key in ('verts', 'edges', 'surfs')}}
    for table_name in ('balls', 'verts', 'edges', 'surfs'):
        table = getattr(net, table_name)
        for relation in ('balls', 'verts', 'edges', 'surfs'):
            if relation in table:
                for values in table[relation]:
                    _ids(values, valid[relation], f'{table_name}.{relation}')
        for field, width in (('loc', 3), ('points', 3), ('tris', 3)):
            if field not in table:
                continue
            for value in table[field]:
                array = np.asarray(value)
                if array.size == 0 and field != 'loc':
                    continue
                expected = (3,) if field == 'loc' else (len(array), width)
                if array.shape != expected or array.dtype.kind not in 'biuf' or not np.all(np.isfinite(array)):
                    raise ArchiveError(f'Invalid {table_name}.{field} shape/dtype/values')
        if table_name == 'surfs':
            for points, triangles in zip(table['points'], table['tris']):
                array = np.asarray(triangles)
                if array.size and (array.dtype.kind not in 'iu' or np.any(array < 0) or np.any(array >= len(points))):
                    raise ArchiveError('Dangling surface triangle vertex')
    _ids(net.group, ball_ids, 'network.group')
    if net.iface_grps is not None:
        for values in net.iface_grps:
            _ids(values, ball_ids, 'network.iface_grps')


def _validate_context(objects, records):
    classes = _classes()
    for obj, record in zip(objects, records):
        kind = record['kind']
        if kind == 'system':
            if not isinstance(obj.balls, pd.DataFrame) or not {'loc', 'rad', 'num', 'res', 'chn'}.issubset(obj.balls.columns):
                raise ArchiveError('Invalid molecular atom table')
            for key, cls in (('res', classes['residue']), ('chn', classes['chain'])):
                if any(value is not None and not isinstance(value, cls) for value in obj.balls[key]):
                    raise ArchiveError(f'Invalid atom {key} reference')
        elif kind == 'group':
            if not isinstance(obj.sys, classes['system']) or (obj.net is not None and not isinstance(obj.net, classes['network'])):
                raise ArchiveError('Invalid group context')
            _ids(obj.ball_ndxs, set(range(len(obj.sys.balls))), 'group.ball_ndxs')
            for key, table_name in (('layer_atoms', None), ('layer_net_atoms', 'balls'), ('layer_surfs', 'surfs'), ('layer_edges', 'edges'), ('layer_verts', 'verts')):
                for values in getattr(obj, key, None) or []:
                    valid = set(range(len(obj.sys.balls))) if table_name is None else set(obj.net.balls['num']) if table_name == 'balls' else set(range(len(getattr(obj.net, table_name))))
                    _ids(values, valid, f'group.{key}')
        elif kind == 'interface':
            if not isinstance(obj.group1, classes['group']) or not isinstance(obj.group2, classes['group']):
                raise ArchiveError('Invalid interface group reference')


def group_from_network(network, atom_indices, name='selection'):
    """Create a selection sharing solved geometry, restricted to complete cells.

    No geometry is rebuilt. Group-specific layers/statistics are derived on
    demand by the ordinary Group analysis/export methods.
    """
    from copy import copy
    from vorpy.src.group import Group
    indices = sorted(set(int(index) for index in atom_indices))
    if not indices:
        raise ArchiveError('Select at least one atom')
    mapping = network.balls.get('system_num', network.balls['num'])
    complete = network.balls.get('complete', pd.Series(False, index=network.balls.index))
    available = {int(system_id): int(topology_id) for system_id, topology_id, done in
                 zip(mapping, network.balls['num'], complete) if bool(done)}
    if any(index not in available for index in indices):
        raise ArchiveError('Selection includes cells that were not completely solved in this archive')
    system = network.sys
    if system.files['dir'] is None:
        system.files['dir'] = str(Path.cwd() / 'VorPy_Output')
    group = Group(system, name=name, atoms=indices, settings=dict(network.settings), make_net=False)
    group.net = copy(network)
    # A new selection is a distinct network view; shallow-copying its archive
    # identity would create two objects with the same persisted identity.
    for name in ('archive_id', 'generator_provenance'):
        if hasattr(group.net, name):
            delattr(group.net, name)
    stable_id(group.net)
    group.net.group = [available[index] for index in indices]
    group.net.group_name = name
    group.net.settings = group.settings
    group.ball_ndxs = indices
    if group not in system.groups:
        system.groups.append(group)
    return group
