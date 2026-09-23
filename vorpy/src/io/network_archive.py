"""Versioned, self-contained solved-network archives (ZIP + JSON + NumPy).

Public API accepts a Network, Group, or System on save and returns the selected
Network on load. Its .sys retains all archived groups/interfaces and networks.
No solve, triangulation, or scientific analysis is invoked by the loader.
"""
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import os
import tempfile
import zipfile

import numpy as np
import pandas as pd

from vorpy.src.io.archive_codec import ArchiveError, Reader, Writer
from vorpy.src.io.archive_fields import (
    GROUP_FIELDS, INTERFACE_FIELDS, NETWORK_FIELDS, SYSTEM_FIELDS, RESIDUE_FIELDS, CHAIN_FIELDS,
)
from vorpy.src.version import __version__

FORMAT = 'VorPy Network Archive'
FORMAT_VERSION = 1
FIELDS = {'system': SYSTEM_FIELDS, 'network': NETWORK_FIELDS, 'group': GROUP_FIELDS,
          'interface': INTERFACE_FIELDS, 'residue': RESIDUE_FIELDS,
          'chain': CHAIN_FIELDS, 'solvent': CHAIN_FIELDS}
LINKS = {
    'system': ('groups', 'ifaces', 'residues', 'chains', 'sol'),
    'network': ('sys',),
    'group': ('sys', 'net', 'interfaces', 'parent_interface'),
    'interface': ('sys', 'net', 'group1', 'group2', 'partial_group1', 'partial_group2', 'water_groups', 'buried_water_groups'),
    'residue': ('sys', 'chain', 'mol'),
    'chain': ('residues',), 'solvent': ('residues',),
}


def _classes():
    from vorpy.src.system import System
    from vorpy.src.network import Network
    from vorpy.src.group import Group
    from vorpy.src.interface import Interface
    from vorpy.src.objects.residue import Residue
    from vorpy.src.objects.chain import Chain, Sol
    return {'system': System, 'network': Network, 'group': Group, 'interface': Interface,
            'residue': Residue, 'chain': Chain, 'solvent': Sol}


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
            if id(value) in ids:
                return
            ids[id(value)] = len(objects)
            objects.append(value)
            kind = types[type(value)]
            for name in LINKS[kind]:
                visit(getattr(value, name, None))
        elif isinstance(value, dict):
            for item in value.values():
                visit(item)
        elif isinstance(value, (list, tuple)):
            for item in value:
                visit(item)
    visit(root)
    return objects, ids, types


def save_network(network, path):
    """Atomically save the root network and its molecular/group/interface context."""
    root = _root(network)
    if root is None or root.sys is None:
        raise ArchiveError('A solved Network with a parent System is required')
    validate_network(root)
    objects, ids, types = _collect(root)
    destination = Path(path).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    source = root.sys.files.get('base_file')
    metadata = {
        'format': FORMAT, 'format_version': FORMAT_VERSION, 'vorpy_version': __version__,
        'created_at': datetime.now(timezone.utc).isoformat(), 'root_network': ids[id(root)],
        'name': root.group_name, 'original_structure_path': source,
        'original_structure_filename': Path(source).name if source else None,
        'units': {'length': 'angstrom', 'area': 'angstrom^2', 'volume': 'angstrom^3',
                  'mean_curvature': 'angstrom^-1', 'gaussian_curvature': 'angstrom^-2'},
        'counts': {'atoms': len(root.sys.balls), **{key: len(getattr(root, key)) for key in ('balls', 'surfs', 'edges', 'verts')}},
        'solve_duration_seconds': root.metrics.get('tot'),
        'scheme': root.settings.get('net_type'),
        'surface_resolution': root.settings.get('surf_res'),
        'boundary_padding': root.settings.get('box_size'),
        'maximum_vertex_radius': root.settings.get('max_vert'),
        'probe_size': root.settings.get('probe_size'),
    }
    fd, temp_name = tempfile.mkstemp(prefix='.vorpy-', suffix='.tmp', dir=destination.parent)
    os.close(fd)
    try:
        with zipfile.ZipFile(temp_name, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=4, allowZip64=True) as archive:
            writer = Writer(archive, ids)
            records = []
            for obj in objects:
                kind = types[type(obj)]
                fields = {name: writer.encode(getattr(obj, name)) for name in FIELDS[kind] + LINKS[kind]
                          if hasattr(obj, name)}
                records.append({'id': ids[id(obj)], 'kind': kind, 'fields': fields})
            metadata['solve_settings'] = writer.encode(root.settings)
            writer.json('state.json', {'objects': records})
            writer.json('metadata.json', metadata)
        os.replace(temp_name, destination)
    finally:
        if os.path.exists(temp_name):
            os.unlink(temp_name)
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
    return value


def _read_v1(reader, metadata):
    records = reader.json('state.json')['objects']
    if not isinstance(records, list) or not records:
        raise ArchiveError('Missing object records')
    for index, record in enumerate(records):
        kind = record['kind']
        if record['id'] != index or kind not in FIELDS or not isinstance(record['fields'], dict):
            raise ArchiveError('Invalid object identity/type')
        if set(record['fields']) - set(FIELDS[kind] + LINKS[kind]):
            raise ArchiveError(f'Unknown fields for {kind}')
        reader.references[index] = _Reference(index)
    # Parse/validate all arrays and references before allocating core objects.
    decoded = [{key: reader.decode(value) for key, value in record['fields'].items()} for record in records]
    classes = _classes()
    objects = [classes[record['kind']].__new__(classes[record['kind']]) for record in records]
    for record, fields, obj in zip(records, decoded, objects):
        for key, value in fields.items():
            setattr(obj, key, _resolve(value, objects))
        kind = record['kind']
        if kind == 'system':
            obj.files = dict.fromkeys(('base_file', 'ball_file', 'verts_file', 'net_file', 'ndx_file', 'frame_files'))
            obj.files.update(dir=None, root_dir=str(Path.cwd()))
            obj.gui = None
            obj.atoms = None
            obj.user_atoms = None
            obj.simple = False
            obj.print_actions = False
            obj.verbose = False
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
    root_id = metadata['root_network']
    if not isinstance(root_id, int) or not 0 <= root_id < len(objects) or records[root_id]['kind'] != 'network':
        raise ArchiveError('Invalid root network ID')
    root = objects[root_id]
    for record, obj in zip(records, objects):
        if record['kind'] == 'network' and all(getattr(obj, key, None) is not None for key in ('balls', 'verts', 'edges', 'surfs')):
            validate_network(obj)
    for key in ('balls', 'surfs', 'edges', 'verts'):
        if metadata['counts'][key] != len(getattr(root, key)):
            raise ArchiveError(f'Entity count mismatch: {key}')
    if metadata['counts']['atoms'] != len(root.sys.balls):
        raise ArchiveError('Atom count mismatch')
    _validate_context(objects, records)
    root.archive_metadata = metadata
    root.sys.archive_metadata = metadata
    return root


# Explicit version dispatch is the insertion point for future schema migrations.
LOADERS = {1: _read_v1}


def load_network(path):
    """Load a self-contained solved network. Never build/analyze geometry."""
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
                raise ArchiveError(f'Unsupported VorPy archive format version: {version!r}; supported: 1')
            root = LOADERS[version](reader, metadata)
            root.archive_path = str(Path(path).resolve())
            root.sys.files['dir'] = str(Path(path).resolve().with_suffix('')) + '_exports'
            return root
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
    group.net.group = [available[index] for index in indices]
    group.net.group_name = name
    group.net.settings = group.settings
    group.ball_ndxs = indices
    if group not in system.groups:
        system.groups.append(group)
    return group
