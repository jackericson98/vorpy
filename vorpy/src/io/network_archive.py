"""Versioned, self-contained solved-network archives (ZIP + JSON + NumPy).

Public API accepts a Network, Group, or System on save and returns the selected
Network on load. Its .sys retains all archived groups/interfaces and networks.
No solve, triangulation, or scientific analysis is invoked by the loader.
"""
from dataclasses import dataclass, fields, replace
from datetime import datetime, timezone
from collections.abc import Mapping
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


def _source_sha256(path):
    """Hash an available source file without making source loading a prerequisite."""
    if not path:
        return None
    try:
        digest = hashlib.sha256()
        with Path(path).open('rb') as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b''):
                digest.update(chunk)
        return digest.hexdigest()
    except (OSError, TypeError, ValueError):
        return None


def validate_network_source(network, source_path, *, expected_scheme=None,
                            expected_solvent_mode=None):
    """Validate a loaded network against its recorded molecular source.

    This is deliberately a source-table check only: it parses the source
    atoms, but never constructs or analyzes a Network.  Recovery callers can
    therefore fail closed before resolving selectors against stale geometry.
    """
    if not source_path:
        raise ArchiveError('Checkpoint does not record an original structure path')
    source = Path(source_path).expanduser().resolve()
    if not source.is_file():
        raise ArchiveError(f'Recovery source file does not exist: {source}')

    metadata = dict(getattr(network, 'archive_metadata', {}) or {})
    recorded_hash = metadata.get('source_sha256')
    actual_hash = _source_sha256(source)
    if not recorded_hash or actual_hash != recorded_hash:
        raise ArchiveError(
            'Recovery source hash does not match the checkpoint: '
            f'recorded={recorded_hash!r}, actual={actual_hash!r}'
        )

    settings = dict(getattr(network, 'settings', {}) or {})
    scheme = str(settings.get('net_type', '')).strip().lower()
    if expected_scheme is not None and scheme != str(expected_scheme).strip().lower():
        raise ArchiveError(
            f'Recovery requires {expected_scheme!r} geometry, found {scheme!r}'
        )
    solvent_mode = settings.get('solvent_mode', getattr(network, 'solvent_mode', None))
    solvent_mode = getattr(solvent_mode, 'value', solvent_mode)
    solvent_mode = str(solvent_mode).strip().upper() if solvent_mode is not None else None
    if expected_solvent_mode is not None and solvent_mode != str(expected_solvent_mode).strip().upper():
        raise ArchiveError(
            f'Recovery requires solvent mode {expected_solvent_mode!r}, '
            f'found {solvent_mode!r}'
        )

    # Local generator order and stable IDs are part of the trust boundary;
    # do not infer or repair them during recovery.
    balls = getattr(network, 'balls', None)
    generator_indices = getattr(network, 'generator_system_indices', None)
    stable_ids = getattr(network, 'generator_stable_ids', None)
    if balls is None or generator_indices is None or stable_ids is None:
        raise ArchiveError('Checkpoint lacks explicit generator identity metadata')
    generator_indices = tuple(
        None if value is None else int(value) for value in generator_indices
    )
    stable_ids = tuple(str(value) for value in stable_ids)
    if len(generator_indices) != len(balls) or len(stable_ids) != len(balls):
        raise ArchiveError('Checkpoint generator metadata disagrees with its balls table')
    if any(value is None for value in generator_indices):
        raise ArchiveError('Recovery source contains non-source boundary generators')
    if len(set(generator_indices)) != len(generator_indices):
        raise ArchiveError('Checkpoint generator system indices are not unique')

    # Import locally to keep archive loading independent of molecular input
    # readers and to make it explicit that this path does not build geometry.
    from vorpy.src.system.system import System

    source_system = System(
        str(source), make_dir=False, print_actions=False,
        solvent_mode=expected_solvent_mode,
    )
    source_balls = source_system.balls
    if source_balls is None or len(source_balls) != len(generator_indices):
        raise ArchiveError(
            'Checkpoint generator count does not match the recovery source: '
            f'{len(generator_indices)} != {0 if source_balls is None else len(source_balls)}'
        )
    if tuple(generator_indices) != tuple(range(len(source_balls))):
        raise ArchiveError('Checkpoint generators are not the ordered source PDB universe')

    expected_stable_ids = tuple(f'system:{index}' for index in generator_indices)
    if stable_ids != expected_stable_ids:
        raise ArchiveError('Checkpoint generator stable IDs do not match source order')

    # A checkpoint may retain the source-generator policy separately from the
    # ball table.  If it is present, it is an assertion about the exact source
    # universe, not advisory presentation metadata.
    provenance = dict(getattr(network, 'archive_provenance', {}) or {})
    generator_policy = provenance.get('source_generator_policy')
    if generator_policy is None:
        generator_policy = settings.get('generator_policy')
    if isinstance(generator_policy, Mapping):
        policy_count = generator_policy.get('source_generator_count')
        if policy_count is not None:
            try:
                policy_count = int(policy_count)
            except (TypeError, ValueError, OverflowError) as error:
                raise ArchiveError('Checkpoint source-generator policy count is unreadable') from error
            if policy_count != len(generator_indices):
                raise ArchiveError(
                    'Checkpoint source-generator policy count does not match its balls table'
                )
        active_indices = generator_policy.get('active_system_indices')
        if active_indices is not None:
            try:
                active_indices = tuple(int(value) for value in active_indices)
            except (TypeError, ValueError, OverflowError) as error:
                raise ArchiveError('Checkpoint active generator policy is unreadable') from error
            if active_indices != generator_indices:
                raise ArchiveError(
                    'Checkpoint active generator policy does not match source order'
                )

    identity_balls = getattr(getattr(network, 'sys', None), 'balls', None)
    if identity_balls is None or len(identity_balls) != len(source_balls):
        raise ArchiveError('Checkpoint lacks the archived source atom table')
    comparable = ('name', 'res_name', 'res_seq', 'chain_id', 'element')
    for local_index, system_index in enumerate(generator_indices):
        archived = identity_balls.iloc[system_index]
        geometry = balls.iloc[local_index]
        original = source_balls.iloc[system_index]
        for column in comparable:
            if str(archived.get(column, '')) != str(original.get(column, '')):
                raise ArchiveError(
                    f'Generator identity mismatch at source index {system_index}: {column}'
                )
        try:
            if not np.array_equal(
                np.asarray(geometry['loc'], dtype=float),
                np.asarray(original['loc'], dtype=float),
            ) or not np.isclose(float(geometry['rad']), float(original['rad'])):
                raise ArchiveError(
                    f'Generator coordinates/radius mismatch at source index {system_index}'
                )
        except (KeyError, TypeError, ValueError) as error:
            raise ArchiveError(
                f'Generator geometry identity is unreadable at source index {system_index}'
            ) from error

    complete = balls.get('complete')
    complete_count = (
        int(complete.astype(bool).sum())
        if complete is not None else None
    )
    return {
        'source_path': str(source),
        'source_sha256': actual_hash,
        'generator_count': len(generator_indices),
        'scheme': scheme,
        'solvent_mode': solvent_mode,
        'boundary_mode': getattr(network, 'boundary_mode', None)
        or settings.get('boundary_mode'),
        'complete_cell_count': complete_count,
        'complete_cell_evidence': complete is not None,
        'generator_policy': generator_policy,
        'geometry_status': provenance.get('geometry_status'),
        'analysis_status': provenance.get('analysis_status'),
        'settings': settings,
    }


def _write_timing_manifest(path, profile, *, fail_if_exists=False):
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
        if fail_if_exists:
            try:
                # Hard-link creation is the Windows/POSIX exclusive-create
                # primitive: a racing writer cannot replace the destination.
                os.link(temp_name, destination)
            except FileExistsError as error:
                raise FileExistsError(
                    f'Refusing to overwrite existing archive timing sidecar: {destination}'
                ) from error
        else:
            os.replace(temp_name, destination)
    finally:
        if os.path.exists(temp_name):
            os.unlink(temp_name)
    return destination

FORMAT = 'VorPy Network Archive'
FORMAT_VERSION = 3
GROUP_LAYER_ATOMS_CONTRACT = 'system_global_nonnegative_v1'
ARCHIVE_PURPOSES = frozenset({
    'interface_results',
    'user_archive',
    'recovery_checkpoint',
    'cache_export',
})
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
# Boundary adapter records are referenced through the object table, not by
# scientific archive identity. Assigning runtime IDs to these frozen config
# values would change solve-context digests during archive preparation.
IDENTITY_ADAPTERS = frozenset(ADAPTERS) - {'boundary_config', 'boundary_generator'}
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
# v3 optional extension: a subordinate network may own boundary generators
# distinct from its parent System. Older v1/v2 and v3 archives omit this field.
FIELDS['network'] += ('boundary_generators',)
# v3 optional extension: preserve producer-attested solvent mode and the
# local-to-System generator membership used by the solved network. These are
# deliberately after the v1/v2 snapshots; absent fields remain unresolved.
FIELDS['system'] += ('solvent_mode',)
FIELDS['network'] += (
    'generator_system_indices', 'generator_stable_ids',
    'solvent_local_indices', 'solvent_generator_ids', 'solvent_mode',
    'solvent_cells_complete', 'solvent_generator_digest',
    'retained_system_indices', 'retained_generator_ids',
    'retained_generator_digest',
)
# v3 optional extension: preserve the producer namespace separately from the
# archive object ID. Missing legacy values remain unavailable on load.
FIELDS['network'] += ('source_network_id',)
# v3 optional extension: preserve checkpoint/analysis state without treating
# a geometry checkpoint as an analysis-complete Result.
FIELDS['network'] += ('archive_provenance',)
# v3 optional extension: preserve the canonical solve dependency for native
# interface vertex state; missing legacy values must fail closed on reuse.
FIELDS['network'] += ('interface_vertex_state_dependency',)
FIELDS['system'] += ('retained_system_indices',)
FIELDS['interface'] += ('representation_caches', 'representation_links')
# v3 optional extension: preserve the typed non-group interface handoff.  The
# payload is a strict detached snapshot, not a reference to executable state.
FIELDS['interface'] += ('typed_results', 'results_dependency')

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


def _archive_provenance_value(value):
    """Convert typed lineage nested in provenance to JSON-safe archive state."""
    from vorpy.src.results import GeometryLineage

    if isinstance(value, GeometryLineage):
        return value.as_dict()
    if isinstance(value, Mapping):
        return {key: _archive_provenance_value(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_archive_provenance_value(item) for item in value]
    if isinstance(value, tuple):
        return tuple(_archive_provenance_value(item) for item in value)
    if isinstance(value, (set, frozenset)):
        return type(value)(_archive_provenance_value(item) for item in value)
    return value


_INTERFACE_RESULT_SNAPSHOT_FORMAT = 'vorpy.interface_result_snapshot'
_RESULTS_DEPENDENCY_SNAPSHOT_FORMAT = 'vorpy.results_dependency_snapshot'


def _result_dictionary(value):
    from vorpy.src.results.model import dictionary
    return dictionary(value)


def _result_snapshot(value):
    from vorpy.src.results import InterfaceResult
    if not isinstance(value, InterfaceResult):
        raise ArchiveError('interface.typed_results must contain InterfaceResult values')
    payload = _result_dictionary(value)
    return {
        'format': _INTERFACE_RESULT_SNAPSHOT_FORMAT,
        'schema_version': 1,
        'result_schema_version': payload.get('schema_version'),
        'payload': payload,
    }


def _dependency_snapshot(value):
    from vorpy.src.results import ResultsDependency
    if not isinstance(value, ResultsDependency):
        raise ArchiveError('interface.results_dependency must be ResultsDependency or None')
    return {
        'format': _RESULTS_DEPENDENCY_SNAPSHOT_FORMAT,
        'schema_version': 1,
        'payload': _result_dictionary(value),
    }


def _mapping(value, label):
    if not isinstance(value, Mapping):
        raise ArchiveError(f'{label} must be a mapping')
    return dict(value)


def _dataclass_payload(cls, value, label, *, converters=None):
    payload = _mapping(value, label)
    allowed = {item.name for item in fields(cls)}
    unknown = set(payload).difference(allowed)
    if unknown:
        raise ArchiveError(f'{label} has unknown fields: {sorted(unknown)}')
    converters = converters or {}
    try:
        return cls(**{
            name: converters[name](item) if name in converters else item
            for name, item in payload.items()
        })
    except (TypeError, ValueError, KeyError) as error:
        raise ArchiveError(f'Invalid {label}: {error}') from error


def _restore_result_provenance(value, label):
    if value is None:
        return None
    from vorpy.src.results import Provenance
    return _dataclass_payload(Provenance, value, label)


def _restore_result_quantity(value, label):
    from vorpy.src.results import Quantity
    return _dataclass_payload(
        Quantity, value, label,
        converters={'status': lambda item: item,
                    'provenance': lambda item: _restore_result_provenance(item, f'{label}.provenance')},
    )


def _restore_result_quantities(value, label):
    payload = _mapping(value, label)
    return {
        name: _restore_result_quantity(item, f'{label}.{name}')
        for name, item in payload.items()
    }


def _restore_interface_result(value):
    from vorpy.src.results import (
        AtomIdentity, ContactRecord, Environment, GaussianCurvature, InterfaceResult,
        MeanCurvature, RepresentationExtension, ResultIdentity, Status, Topology,
    )
    snapshot = _mapping(value, 'interface result snapshot')
    if snapshot.get('format') != _INTERFACE_RESULT_SNAPSHOT_FORMAT:
        raise ArchiveError('Unsupported interface result snapshot format')
    if snapshot.get('schema_version') != 1:
        raise ArchiveError('Unsupported interface result snapshot schema version')
    payload = _mapping(snapshot.get('payload'), 'interface result snapshot payload')
    if snapshot.get('result_schema_version') != payload.get('schema_version'):
        raise ArchiveError('Interface result snapshot schema version mismatch')
    expected = {item.name for item in fields(InterfaceResult)}
    unknown = set(payload).difference(expected)
    if unknown:
        raise ArchiveError(f'Interface result snapshot has unknown fields: {sorted(unknown)}')

    def _quantity(value, name):
        return _restore_result_quantity(value, f'interface result {name}')

    def _quantity_map(value, name):
        return _restore_result_quantities(value, f'interface result {name}')

    def _environment(value):
        if value is None:
            return None
        return _dataclass_payload(
            Environment, value, 'interface result environment',
            converters={'provenance': lambda item: _restore_result_provenance(
                item, 'interface result environment.provenance')},
        )

    def _topology(value):
        if value is None:
            return None
        return _dataclass_payload(
            Topology, value, 'interface result topology',
            converters={'quantities': lambda item: _quantity_map(item, 'topology.quantities')},
        )

    def _mean(value):
        if value is None:
            return None
        return _dataclass_payload(
            MeanCurvature, value, 'interface result mean_curvature',
            converters={name: lambda item, name=name: _quantity(item, f'mean_curvature.{name}')
                        for name in (
                            'integrated_mean_curvature', 'complete_total',
                            'smooth_face_contribution', 'edge_contribution',
                            'boundary_contribution',
                        )},
        )

    def _gaussian(value):
        if value is None:
            return None
        return _dataclass_payload(
            GaussianCurvature, value, 'interface result gaussian_curvature',
            converters={name: lambda item, name=name: _quantity(item, f'gaussian_curvature.{name}')
                        for name in (
                            'integrated_gaussian_curvature', 'smooth_face_contribution',
                            'intrinsic_discrete_contribution', 'intrinsic_seam_contribution',
                            'interior_vertex_defects', 'boundary_geodesic_contribution',
                            'boundary_corner_contribution', 'gauss_bonnet_completion',
                            'gauss_bonnet_expected', 'gauss_bonnet_residual',
                        )},
        )

    def _contact(value):
        payload = _mapping(value, 'interface result contact')
        return _dataclass_payload(
            ContactRecord, payload, 'interface result contact',
            converters={
                'identity': lambda item: _dataclass_payload(
                    ResultIdentity, item, 'interface result contact.identity'),
                'atom_A': lambda item: _dataclass_payload(
                    AtomIdentity, item, 'interface result contact.atom_A'),
                'atom_B': lambda item: _dataclass_payload(
                    AtomIdentity, item, 'interface result contact.atom_B'),
                'quantities': lambda item: _quantity_map(item, 'contact.quantities'),
                'provenance': lambda item: _restore_result_provenance(
                    item, 'interface result contact.provenance'),
            },
        )

    def _extension(value):
        return _dataclass_payload(
            RepresentationExtension, value, 'interface result extension',
            converters={'provenance': lambda item: _restore_result_provenance(
                item, 'interface result extension.provenance')},
        )

    converters = {
        'identity': lambda item: _dataclass_payload(ResultIdentity, item, 'interface result identity'),
        'status': lambda item: item,
        'provenance': lambda item: _restore_result_provenance(item, 'interface result provenance'),
        'environment': _environment,
        'metrics': lambda item: _quantity_map(item, 'metrics'),
        'topology': _topology,
        'mean_curvature': _mean,
        'gaussian_curvature': _gaussian,
        'contacts': lambda item: tuple(_contact(entry) for entry in item)
        if isinstance(item, (list, tuple)) else (_ for _ in ()).throw(
            ArchiveError('interface result contacts must be a list')),
        'selection': lambda item: _quantity_map(item, 'selection'),
        'extensions': lambda item: tuple(_extension(entry) for entry in item)
        if isinstance(item, (list, tuple)) else (_ for _ in ()).throw(
            ArchiveError('interface result extensions must be a list')),
    }
    return _dataclass_payload(InterfaceResult, payload, 'interface result', converters=converters)


def _restore_dependency(value):
    from vorpy.src.results import ResultsDependency, SharedNetworkIdentity, Status
    if value is None:
        return None
    snapshot = _mapping(value, 'results dependency snapshot')
    if snapshot.get('format') != _RESULTS_DEPENDENCY_SNAPSHOT_FORMAT:
        raise ArchiveError('Unsupported results dependency snapshot format')
    if snapshot.get('schema_version') != 1:
        raise ArchiveError('Unsupported results dependency snapshot schema version')
    payload = _mapping(snapshot.get('payload'), 'results dependency snapshot payload')
    return _dataclass_payload(
        ResultsDependency, payload, 'results dependency',
        converters={
            'identity': lambda item: _dataclass_payload(
                SharedNetworkIdentity, item, 'results dependency identity'),
            'status': lambda item: Status(item),
        },
    )


def _archive_field_value(kind, name, value):
    if kind == 'interface' and name == 'typed_results':
        if value is None:
            return ()
        if not isinstance(value, (list, tuple)):
            raise ArchiveError('interface.typed_results must be a sequence')
        return tuple(_result_snapshot(item) for item in value)
    if kind == 'interface' and name == 'results_dependency':
        return _dependency_snapshot(value) if value is not None else None
    if name == 'solvent_mode' and value is not None:
        # System stores the enum while Network stores its canonical string.
        return getattr(value, 'value', value)
    if kind == 'interface' and name == 'water_topology':
        return _water_topology_snapshot(value)
    if kind == 'interface' and name == 'water_geometries' and value is not None:
        return [dict(geometry, topology_analysis=_water_topology_snapshot(geometry['topology_analysis']))
                if isinstance(geometry, dict) and 'topology_analysis' in geometry else geometry
                for geometry in value]
    if name in {'archive_provenance', 'analysis_provenance'}:
        return _archive_provenance_value(value)
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


def _adopt_typed_result_archive_ids(objects, types):
    """Bind a producer-issued archive ID before generic IDs are allocated."""
    from vorpy.src.results import InterfaceResult

    for iface in objects:
        if types[type(iface)] != 'interface':
            continue
        network = getattr(iface, 'net', None)
        candidates = {
            result.identity.network_id
            for result in (getattr(iface, 'typed_results', ()) or ())
            if isinstance(result, InterfaceResult)
        }
        if not candidates:
            continue
        if len(candidates) != 1:
            raise ArchiveError('interface typed results disagree on archive identity')
        candidate = candidates.pop()
        existing = getattr(network, 'archive_id', None)
        if existing is None:
            object.__setattr__(network, 'archive_id', candidate)


def save_network(network, path, *, progress=None, fail_if_exists=False,
                 write_timing_manifest=None, archive_purpose=None):
    """Atomically save the root network and its molecular/group/interface context."""
    if archive_purpose is not None and archive_purpose not in ARCHIVE_PURPOSES:
        raise ValueError(f'Unknown archive purpose: {archive_purpose!r}')
    root = _root(network)
    if root is None or root.sys is None:
        raise ArchiveError('A solved Network with a parent System is required')
    if write_timing_manifest is None:
        write_timing_manifest = getattr(root.sys, '_write_archive_timing_manifest', True)
    validate_network(root)
    return _save_document(
        root, path, root, 'network', progress=progress,
        fail_if_exists=fail_if_exists,
        write_timing_manifest=write_timing_manifest,
        archive_purpose=archive_purpose,
    )


def save_session(session, path, *, progress=None, fail_if_exists=False,
                 write_timing_manifest=None, archive_purpose=None):
    """Save a complete session in the same .vpy object graph format."""
    if archive_purpose is not None and archive_purpose not in ARCHIVE_PURPOSES:
        raise ValueError(f'Unknown archive purpose: {archive_purpose!r}')
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
    if write_timing_manifest is None:
        owner = getattr(getattr(session, 'active_network', None), 'sys', None)
        write_timing_manifest = getattr(owner, '_write_archive_timing_manifest', True)
    return _save_document(
        session, path, session.active_network, 'session', progress=progress,
        fail_if_exists=fail_if_exists,
        write_timing_manifest=write_timing_manifest,
        archive_purpose=archive_purpose,
    )


def _save_document(document, path, root, document_kind, *, progress=None,
                   fail_if_exists=False, write_timing_manifest=True,
                   archive_purpose=None):
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
    destination = Path(path).expanduser().resolve()
    if fail_if_exists and (
        destination.exists()
        or _timing_manifest_path(destination).exists()
    ):
        raise FileExistsError(
            f'Refusing to overwrite existing archive or timing sidecar: {destination}'
        )
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
    _validate_group_layer_atoms(objects, types, allow_legacy_sentinel=False)
    _adopt_typed_result_archive_ids(objects, types)
    started = perf_counter()
    for obj in objects:
        kind = types[type(obj)]
        normalize_vocabulary(obj)
        if kind in {'system', 'network', 'group', 'interface', 'session'} or kind in IDENTITY_ADAPTERS:
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
            _validate_typed_result_owners(obj)
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
    destination.parent.mkdir(parents=True, exist_ok=True)
    source = root.sys.files.get('base_file') if root is not None else None
    if not source and root is not None:
        # Reloaded systems may not retain runtime file handles in ``sys.files``.
        # Preserve an already-attested source path; do not infer a new one.
        source = (getattr(root, 'archive_metadata', {}) or {}).get(
            'original_structure_path'
        )
    solve_universe_id = None
    if root is not None and hasattr(root, 'get_solve_context_key'):
        try:
            solve_universe_id = root.get_solve_context_key().solve_universe_id
        except (AttributeError, TypeError, ValueError):
            # Legacy/incomplete archives remain serializable; recovery will
            # fail closed if a caller requires a solve-universe identity.
            solve_universe_id = None
    metadata = {
        'format': FORMAT, 'format_version': FORMAT_VERSION, 'vorpy_version': __version__,
        'producer_contracts': {'group.layer_atoms': GROUP_LAYER_ATOMS_CONTRACT},
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
    if solve_universe_id is not None:
        metadata['solve_universe_id'] = str(solve_universe_id)
    metadata['source_sha256'] = _source_sha256(source)
    if archive_purpose is not None:
        metadata['archive_purpose'] = archive_purpose
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
        if fail_if_exists:
            try:
                # os.link() fails atomically when another producer wins the
                # same destination between preflight and finalization.
                os.link(temp_name, destination)
            except FileExistsError as error:
                raise FileExistsError(
                    f'Refusing to overwrite existing archive: {destination}'
                ) from error
        else:
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
    profile['timing_manifest'] = (
        str(_timing_manifest_path(destination)) if write_timing_manifest else None
    )
    if write_timing_manifest:
        _write_timing_manifest(
            destination, profile, fail_if_exists=fail_if_exists
        )
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
            archived_solvent_mode = hasattr(obj, 'solvent_mode')
            if archived_solvent_mode:
                from vorpy.src.boundary import normalize_solvent_mode
                obj.solvent_mode = normalize_solvent_mode(obj.solvent_mode)
            else:
                # Legacy archives do not attest a solvent generator universe.
                obj.solvent_mode = None
            if not hasattr(obj, 'retained_system_indices'):
                obj.retained_system_indices = ()
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
            if not hasattr(obj, 'source_network_id'):
                obj.source_network_id = None
            if not hasattr(obj, 'generator_system_indices'):
                obj.generator_system_indices = None
            if not hasattr(obj, 'generator_stable_ids'):
                obj.generator_stable_ids = None
            if not hasattr(obj, 'solvent_local_indices'):
                obj.solvent_local_indices = ()
            if not hasattr(obj, 'solvent_generator_ids'):
                obj.solvent_generator_ids = ()
            if not hasattr(obj, 'solvent_mode'):
                obj.solvent_mode = None
            if not hasattr(obj, 'solvent_cells_complete'):
                obj.solvent_cells_complete = None
            if not hasattr(obj, 'solvent_generator_digest'):
                obj.solvent_generator_digest = None
            if not hasattr(obj, 'retained_system_indices'):
                obj.retained_system_indices = ()
            if not hasattr(obj, 'retained_generator_ids'):
                obj.retained_generator_ids = ()
            if not hasattr(obj, 'retained_generator_digest'):
                obj.retained_generator_digest = None
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
                obj.typed_results = tuple(
                    _restore_interface_result(item)
                    for item in getattr(obj, 'typed_results', ()) or ()
                )
                obj.results_dependency = _restore_dependency(
                    getattr(obj, 'results_dependency', None)
                )
    for record, obj in zip(records, objects):
        kind = record['kind']
        if kind in {'system', 'network', 'group', 'interface', 'session'} or kind in IDENTITY_ADAPTERS:
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
            _validate_typed_result_owners(obj, fail_closed=True)
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
        _validate_context(objects, records, metadata)
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
    _validate_context(objects, records, metadata)
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


def _validate_typed_result_owners(iface, *, fail_closed=False):
    """Validate typed interface ownership and optionally fail closed on load."""
    from vorpy.src.results import (
        InterfaceResult,
        ResultsDependency,
        Status,
        validate_results_dependency,
    )

    results = getattr(iface, 'typed_results', ()) or ()
    if not isinstance(results, (list, tuple)):
        raise ArchiveError('interface.typed_results must be a sequence')
    archive_id = getattr(iface.net, 'archive_id', None)
    source_id = getattr(iface.net, 'source_network_id', None)
    unresolved = []
    for result in results:
        if not isinstance(result, InterfaceResult):
            raise ArchiveError('interface.typed_results must contain InterfaceResult values')
        if result.identity.network_id != archive_id:
            raise ArchiveError('interface typed result network identity disagrees with archive')
        if result.identity.source_network_id == archive_id:
            raise ArchiveError('interface typed result producer identity aliases archive identity')
        if source_id is not None and result.identity.source_network_id != source_id:
            raise ArchiveError('interface typed result producer identity disagrees with network')
        if result.identity.interface_id != getattr(iface, 'interface_id', None):
            raise ArchiveError('interface typed result interface identity disagrees with archive')
        if result.identity.source_interface_id not in (None, getattr(iface, 'interface_id', None)):
            raise ArchiveError('interface typed result source interface disagrees with archive')
        group1_id = getattr(getattr(iface, 'group1', None), 'group_id', None)
        group2_id = getattr(getattr(iface, 'group2', None), 'group_id', None)
        if group1_id is not None and result.identity.group_A != group1_id:
            raise ArchiveError('interface typed result Group A disagrees with archive')
        if group2_id is not None and result.identity.group_B != group2_id:
            raise ArchiveError('interface typed result Group B disagrees with archive')
        if not all((
            result.identity.group_A,
            result.identity.group_B,
            result.identity.orientation,
            result.identity.selection_id,
        )):
            unresolved.append(result)
            continue
        lineage = result.provenance.geometry_lineage
        if lineage is None or lineage.origin.value == 'UNRESOLVED':
            unresolved.append(result)
            continue
        if lineage.source_network_id != result.identity.source_network_id:
            unresolved.append(result)
        if lineage.solve_universe_id != result.identity.solve_universe_id:
            unresolved.append(result)
        if lineage.source_interface_id not in (None, getattr(iface, 'interface_id', None)):
            unresolved.append(result)
    dependency = getattr(iface, 'results_dependency', None)
    if dependency is not None:
        if not isinstance(dependency, ResultsDependency):
            raise ArchiveError('interface.results_dependency must be ResultsDependency or None')
        if results and dependency.identity.source_network_id != results[0].identity.source_network_id:
            raise ArchiveError('interface dependency producer identity disagrees with result')
        if results and dependency.identity.solve_universe_id != results[0].identity.solve_universe_id:
            raise ArchiveError('interface dependency solve identity disagrees with result')
        if results and dependency.identity.selection_id != results[0].identity.selection_id:
            unresolved.extend(results)
        if results:
            expected_orientation = {
                'group1 -> group2': 'A -> B',
                'group2 -> group1': 'B -> A',
            }.get(results[0].identity.orientation)
            if expected_orientation is None or dependency.identity.orientation != expected_orientation:
                unresolved.extend(results)
            if validate_results_dependency(dependency, dependency) in {
                Status.UNRESOLVED,
                Status.NOT_CALCULATED,
                Status.NOT_SUPPORTED,
            }:
                unresolved.extend(results)
    elif results:
        unresolved.extend(results)
    if fail_closed and unresolved:
        unresolved_ids = {id(result) for result in unresolved}
        iface.typed_results = tuple(
            result if id(result) not in unresolved_ids
            else replace(result, status=Status.UNRESOLVED)
            for result in results
        )


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
        if (isinstance(value, (bool, np.bool_))
                or not isinstance(value, (int, np.integer)) or int(value) not in valid):
            raise ArchiveError(f'{label}: invalid/dangling ID {value!r}')


def _validate_generator_metadata(net):
    """Validate explicit producer membership without repairing it on load."""
    ball_count = len(net.balls)
    boundary_values = getattr(net, 'boundary_indices', None)
    boundary_indices = set()
    if boundary_values is not None:
        if not isinstance(boundary_values, (list, tuple, set, np.ndarray)):
            raise ArchiveError('network.boundary_indices: expected a list of IDs')
        for value in boundary_values:
            if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
                raise ArchiveError(f'network.boundary_indices: invalid ID {value!r}')
            value = int(value)
            if value < 0 or value >= ball_count:
                raise ArchiveError(f'network.boundary_indices: invalid/dangling ID {value!r}')
            boundary_indices.add(value)
        if len(boundary_indices) != len(boundary_values):
            raise ArchiveError('network.boundary_indices: duplicate IDs')

    boundary_generators = getattr(net, 'boundary_generators', None)
    if boundary_generators is not None and not isinstance(
            boundary_generators, (list, tuple, np.ndarray)):
        raise ArchiveError('network.boundary_generators: expected a list of generators')
    if boundary_generators is not None and len(boundary_generators):
        if len(boundary_generators) != len(boundary_indices):
            raise ArchiveError('network.boundary_generators: count disagrees with boundary_indices')
        for generator in boundary_generators:
            position = getattr(generator, 'position', None)
            radius = getattr(generator, 'radius', None)
            position_array = np.asarray(position) if position is not None else np.asarray(())
            if (position_array.shape != (3,) or position_array.dtype.kind not in 'biuf'
                    or not np.all(np.isfinite(position_array))):
                raise ArchiveError('network.boundary_generators: invalid position')
            if (isinstance(radius, bool) or not isinstance(radius, (int, float, np.number))
                    or not np.isfinite(radius) or float(radius) <= 0):
                raise ArchiveError('network.boundary_generators: invalid radius')

    generator_system_indices = getattr(net, 'generator_system_indices', None)
    if generator_system_indices is not None:
        if not isinstance(generator_system_indices, (list, tuple, np.ndarray)):
            raise ArchiveError('network.generator_system_indices: expected a list of IDs')
        if len(generator_system_indices) != ball_count:
            raise ArchiveError('network.generator_system_indices: length disagrees with balls')
        system_numbers = net.balls.get('system_num')
        for local_index, value in enumerate(generator_system_indices):
            is_boundary = local_index in boundary_indices
            if value is None:
                if not is_boundary:
                    raise ArchiveError(
                        f'network.generator_system_indices: missing physical generator ID at {local_index}'
                    )
                expected = -1
            else:
                if (is_boundary or isinstance(value, bool)
                        or not isinstance(value, (int, np.integer)) or int(value) < 0):
                    raise ArchiveError(
                        f'network.generator_system_indices: invalid physical generator ID {value!r}'
                    )
                expected = int(value)
            if system_numbers is not None:
                actual = system_numbers.iloc[local_index]
                if pd.isna(actual):
                    # Older/manual producers used a missing boundary marker;
                    # preserve it rather than manufacturing the F10 -1 value.
                    if not is_boundary:
                        raise ArchiveError(
                            f'network.generator_system_indices: disagrees with balls at {local_index}'
                        )
                elif int(actual) != expected:
                    raise ArchiveError(
                        f'network.generator_system_indices: disagrees with balls at {local_index}'
                    )

    generator_stable_ids = getattr(net, 'generator_stable_ids', None)
    if generator_stable_ids is not None:
        if (not isinstance(generator_stable_ids, (list, tuple, np.ndarray))
                or len(generator_stable_ids) != ball_count
                or any(not isinstance(value, str) or not value for value in generator_stable_ids)):
            raise ArchiveError('network.generator_stable_ids: invalid generator IDs')
        if len(set(generator_stable_ids)) != len(generator_stable_ids):
            raise ArchiveError('network.generator_stable_ids: duplicate generator IDs')

    retained_indices = getattr(net, 'retained_system_indices', None)
    if retained_indices is not None:
        for value in retained_indices:
            if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or int(value) < 0:
                raise ArchiveError(f'network.retained_system_indices: invalid ID {value!r}')


def validate_network(net):
    """Validate table shapes, topology ID spaces, and mesh-local connectivity."""
    source_network_id = getattr(net, 'source_network_id', None)
    if source_network_id is not None and (
        not isinstance(source_network_id, str) or not source_network_id.strip()
    ):
        raise ArchiveError('network.source_network_id must be a non-empty string or None')
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
    _validate_generator_metadata(net)


def _allow_legacy_layer_sentinel(metadata):
    if 'producer_contracts' not in metadata:
        return True
    contracts = metadata['producer_contracts']
    if not isinstance(contracts, dict):
        raise ArchiveError('Invalid producer contract metadata')
    if 'group.layer_atoms' not in contracts:
        raise ArchiveError('Missing group.layer_atoms producer contract')
    contract = contracts['group.layer_atoms']
    if contract != GROUP_LAYER_ATOMS_CONTRACT:
        raise ArchiveError('Unsupported group.layer_atoms producer contract')
    return False


def _validate_group_layer_atoms(objects, kinds, *, allow_legacy_sentinel):
    system_class = _classes()['system']
    for obj in objects:
        if kinds[type(obj)] != 'group':
            continue
        if not isinstance(getattr(obj, 'sys', None), system_class):
            raise ArchiveError('Invalid group context')
        valid = set(range(len(obj.sys.balls)))
        if allow_legacy_sentinel:
            valid.add(-1)
        for values in getattr(obj, 'layer_atoms', None) or []:
            _ids(values, valid, 'group.layer_atoms')


def _validate_context(objects, records, metadata):
    classes = _classes()
    allow_legacy_layer_sentinel = _allow_legacy_layer_sentinel(metadata)
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
                    if key == 'layer_atoms':
                        valid = set(range(len(obj.sys.balls)))
                        if allow_legacy_layer_sentinel:
                            valid.add(-1)
                    elif table_name == 'balls':
                        valid = set(obj.net.balls['num'])
                    else:
                        valid = set(range(len(getattr(obj.net, table_name))))
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
