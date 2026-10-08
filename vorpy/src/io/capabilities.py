"""Read-only availability/provenance inventory over scientific session objects."""
from dataclasses import dataclass
from enum import Enum

from .scientific_adapters import canonical_environment, canonical_network_scope


class CapabilityState(str, Enum):
    AVAILABLE = 'AVAILABLE'
    PARTIAL = 'PARTIAL'
    UNRESOLVED = 'UNRESOLVED'
    NOT_CALCULATED = 'NOT_CALCULATED'
    NOT_SUPPORTED = 'NOT_SUPPORTED'
    STALE = 'STALE'


@dataclass(frozen=True)
class Capability:
    network: str | None
    interface: str | None
    environment: str
    partition: str
    representation: str
    analysis_query: float | None
    result: str
    state: CapabilityState
    reason: str = ''
    network_scope: str = ''
    analysis_id: str | None = None
    provenance_issues: tuple[str, ...] = ()


def _capability_state(value):
    if isinstance(value, CapabilityState):
        return value
    try:
        return CapabilityState(str(value).upper())
    except (TypeError, ValueError):
        return None


def _representation_info(cache, link):
    info = dict(link) if isinstance(link, dict) else {}
    metadata = cache.get('metadata') if isinstance(cache, dict) else getattr(cache, 'metadata', None)
    if isinstance(metadata, dict):
        for key, value in metadata.items():
            info.setdefault(key, value)
    return info


def _dependency_status(dependencies):
    if not isinstance(dependencies, dict):
        return None, ''
    priority = (CapabilityState.NOT_SUPPORTED, CapabilityState.NOT_CALCULATED,
                CapabilityState.UNRESOLVED, CapabilityState.PARTIAL, CapabilityState.STALE)
    for wanted in priority:
        for name, detail in dependencies.items():
            raw = detail.get('state') if isinstance(detail, dict) else detail
            state = _capability_state(raw)
            if state is wanted:
                reason = detail.get('reason') if isinstance(detail, dict) else None
                if not reason:
                    reason = {
                        'selected_ab_contacts': 'selected contacts unavailable',
                        'spherical_clipping': 'clipping geometry unavailable',
                        'atom_centers_radii': 'atom centers/radii unavailable',
                    }.get(name, f'{name} unavailable')
                return state, reason
    return None, ''


def _cached_representation_capabilities(iface, net, environment, partition, scope_name):
    caches = getattr(iface, 'representation_caches', None) or {}
    links = getattr(iface, 'representation_links', None) or {}
    for key in sorted(set(caches) | set(links), key=str):
        if not isinstance(key, str):
            continue
        cache = caches.get(key)
        info = _representation_info(cache, links.get(key, {}))
        state = _capability_state(info.get('state', info.get('status')))
        if state is None:
            state = CapabilityState.AVAILABLE if cache is not None else CapabilityState.NOT_CALCULATED
        if cache is None and state is CapabilityState.AVAILABLE:
            state = CapabilityState.NOT_CALCULATED
        dependency_state, dependency_reason = _dependency_status(info.get('dependencies'))
        if dependency_state is not None and (cache is None or state is CapabilityState.AVAILABLE):
            state = dependency_state
        reason = info.get('reason', '') or dependency_reason
        if info.get('stale'):
            state, reason = CapabilityState.STALE, reason or 'Cached representation is stale'
        issues = info.get('provenance_issues', ())
        yield Capability(
            getattr(net, 'archive_id', None), getattr(iface, 'interface_id', None),
            canonical_environment(info.get('environment', environment)),
            info.get('partition', partition), key,
            info.get('analysis_query', info.get('alpha_value')),
            'geometry', state, reason,
            canonical_network_scope(info.get('network_scope', scope_name)),
            info.get('analysis_id'), tuple(issues) if isinstance(issues, (list, tuple)) else (),
        )


def inventory(session, **filters):
    """Inspect cached records only. Filters never create a solve/query/result.

    Missing requested environments/partitions are NOT_CALCULATED; Dual A/B
    remain NOT_SUPPORTED. Existing incidence duals use representation='dual'.
    """
    from .scientific_adapters import canonical_environment, canonical_network_scope, scope
    filters = dict(filters)
    if 'environment' in filters:
        filters['environment'] = canonical_environment(filters['environment'])
    if 'network_scope' in filters:
        filters['network_scope'] = canonical_network_scope(filters['network_scope'])
    rows = []
    seen = set()
    for system in session.systems:
        for iface in getattr(system, 'ifaces', None) or []:
            if id(iface) in seen:
                continue
            seen.add(id(iface))
            net = getattr(iface, 'net', None)
            analysis = getattr(iface, 'geometry_analysis', None)
            metadata = getattr(analysis, 'metadata', {})
            raw_environment = metadata.get('environment') if analysis else None
            if raw_environment is None and analysis:
                raw_environment = 'solvent_competing' if metadata.get('solvent_context') else 'dry'
            env = canonical_environment(raw_environment)
            partition = {'aw': 'AW', 'pow': 'Power', 'prm': 'Primitive'}.get(
                getattr(net, 'settings', {}).get('net_type'), 'unknown')
            alpha = metadata.get('alpha_value')
            selection = getattr(analysis, 'alpha_selection', {})
            side = getattr(analysis, 'voronoi_side_1', {})
            rep = getattr(iface, '_geometry_representation', None)
            dual = getattr(iface, '_geometry_dual', None)
            filtration = getattr(iface, '_geometry_filtration', None)
            issues = []
            for name, cache in (('dual', dual), ('filtration', filtration)):
                if cache is not None and cache.network is not net:
                    issues.append(f'{name} references another network')
            if rep is not None and rep.alpha != alpha:
                issues.append('representation and analysis alpha differ')
            if rep is not None and selection.get('selected_physical_pairs') is not None:
                if selection['selected_physical_pairs'] != len(rep.selected_surfaces):
                    issues.append('cached selected physical count differs from representation')
            key = getattr(iface, '_geometry_analysis_key', None)
            stale = getattr(iface, 'archive_cache_state', None) == 'STALE'
            if key is not None and net is not None and alpha is not None:
                from vorpy.src.interface.geometry_analysis import geometry_cache_key
                # All stable identities exist after archive load. Avoid assigning one
                # during an inventory of a never-saved live object.
                if getattr(net, 'archive_id', None) is not None:
                    curvature_requested = bool(metadata.get('curvature_requested', True))
                    expected_key = (*geometry_cache_key(iface, alpha), curvature_requested)
                    stale = stale or key != expected_key

            def add(representation, result, exists, status=None, certified=None, apply_stale=True):
                state = CapabilityState.NOT_CALCULATED
                reason = 'Required cached result is absent'
                if exists:
                    state, reason = CapabilityState.AVAILABLE, ''
                    if status == 'unresolved':
                        state = CapabilityState.UNRESOLVED
                    elif status == 'partial' or certified is False:
                        state = CapabilityState.PARTIAL
                    if stale and apply_stale:
                        state, reason = CapabilityState.STALE, 'Cached analysis does not match current source/query'
                rows.append(Capability(getattr(net, 'archive_id', None),
                    getattr(iface, 'interface_id', None), env, partition, representation, alpha,
                    result, state, reason,
                    canonical_network_scope(getattr(net, 'network_scope', scope(net))) if net else '',
                    getattr(analysis, 'archive_id', None), tuple(issues)))

            add('Voronoi', 'network_geometry', net is not None and getattr(net, 'surfs', None) is not None,
                apply_stale=False)
            add('Voronoi', 'geometry', rep is not None, selection.get('status'))
            add('Voronoi', 'analysis', analysis is not None, selection.get('status'))
            add('Voronoi', 'mapping', rep is not None, selection.get('status'))
            add('Voronoi', 'mean_curvature', analysis is not None and bool(side),
                'unresolved' if side.get('partial_integrated_mean_curvature', side.get('partial_H')) is None
                and side.get('total_integrated_mean_curvature', side.get('total_H')) is None else None,
                side.get('mean_curvature_certified', side.get('mean_certified')))
            add('Voronoi', 'gaussian_curvature', analysis is not None and bool(side),
                'unresolved' if side.get('total_integrated_gaussian_curvature', side.get('total_K')) is None else None,
                side.get('gaussian_curvature_certified', side.get('gaussian_total_certified')))
            add('dual', 'geometry', dual is not None,
                'partial' if dual is not None and dual.unresolved_features else None)
            add('dual', 'filtration', filtration is not None,
                'partial' if filtration is not None and (filtration.blocked or
                    any(not row.supported or row.filtration_birth is None
                        for records in filtration.records.values() for row in records.values())) else None)
            for representation in ('Dual A', 'Dual B'):
                rows.append(Capability(getattr(net, 'archive_id', None), getattr(iface, 'interface_id', None),
                    env, partition, representation, alpha, 'geometry', CapabilityState.NOT_SUPPORTED,
                'Representation is not implemented',
                canonical_network_scope(getattr(net, 'network_scope', scope(net))) if net else ''))
            rows.extend(_cached_representation_capabilities(
                iface, net, env, partition,
                canonical_network_scope(getattr(net, 'network_scope', scope(net))) if net else ''))
    interface_networks = {id(getattr(iface, 'net', None)) for system in session.systems
                          for iface in getattr(system, 'ifaces', None) or []}
    for net in session.networks:
        if id(net) in interface_networks:
            continue
        partition = {'aw': 'AW', 'pow': 'Power', 'prm': 'Primitive'}.get(net.settings.get('net_type'), 'unknown')
        for representation, result, exists in (('Voronoi', 'geometry', getattr(net, 'surfs', None) is not None),
                ('Voronoi', 'analysis', False), ('dual', 'geometry', False), ('dual', 'filtration', False)):
            rows.append(Capability(getattr(net, 'archive_id', None), None, 'unknown', partition,
                representation, None, result, CapabilityState.AVAILABLE if exists else CapabilityState.NOT_CALCULATED,
                '' if exists else 'Required cached result is absent',
                canonical_network_scope(getattr(net, 'network_scope', scope(net)))))
    selected = [row for row in rows if all(getattr(row, name) == value for name, value in filters.items())]
    if not selected and filters:
        return [Capability(filters.get('network'), filters.get('interface'),
            canonical_environment(filters.get('environment', 'unknown')),
            filters.get('partition', 'unknown'), filters.get('representation', 'Voronoi'),
            filters.get('analysis_query'), filters.get('result', 'analysis'),
            CapabilityState.NOT_SUPPORTED if filters.get('representation') in {'Dual A', 'Dual B'}
            else CapabilityState.NOT_CALCULATED, 'Requested cached combination is absent')]
    return selected
