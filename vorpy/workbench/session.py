"""Workbench bridge to the shared scientific .vpy session archive."""
from pathlib import Path

from vorpy.src.io import Session, load_session, save_session
from vorpy.src.io.scientific_adapters import stable_id
from vorpy.workbench.project import Project


def session_from_workbench(project, results, active_source=None, structure_states=None,
                           trajectory_results=None, network_sizes_by_frame=None):
    systems, networks = [], []
    cached_results = list(results.values()) + list((trajectory_results or {}).values())
    for result in cached_results:
        owner = result.export_group
        system = getattr(owner, 'sys', None)
        if system is None:
            continue
        if all(system is not item for item in systems):
            systems.append(system)
        for group in (system.groups or []) + (system.ifaces or []):
            net = getattr(group, 'net', None)
            if net is not None and all(net is not item for item in networks):
                networks.append(net)
    current = results.get(active_source)
    net = getattr(getattr(current, 'export_group', None), 'net', None)
    session = Session(systems, networks, active_network_id=stable_id(net) if net is not None else None)
    # Live scientific references stay in AnalysisResult.export_group; arrays are
    # packed by the same codec, not duplicated into result_state JSON.
    project.result_state = None
    session.workbench = {'project': project, 'results': dict(results),
        'active_source': active_source, 'structure_states': dict(structure_states or {}),
        'trajectory_results': dict(trajectory_results or {}),
        'network_sizes_by_frame': dict(network_sizes_by_frame or {})}
    return session


def workbench_from_session(session, source):
    """Restore display state, wrapping old network archives without analysis."""
    state = session.workbench
    if 'project' in state and 'results' in state:
        for result in state['results'].values():
            network = getattr(result.export_group, 'net', None)
            if network is not None:
                result.info_sections.update(cached_info_sections(network))
        return state['project'], state['results'], state.get('active_source'), state.get('structure_states', {})
    from vorpy.workbench.domain import AnalysisResult
    from vorpy.workbench.services.vorpy_backend import _atoms_from_system, _layers_from_network
    project = Project(name=Path(source).stem)
    results = {}
    active_source = None
    for index, network in enumerate(session.networks):
        key = Path(source).resolve() if index == 0 else Path(source).resolve().with_name(
            f'{Path(source).stem}__{network.archive_id}.vpy')
        owner = next((group for group in (network.sys.groups or []) if group.net is network), None)
        result = AnalysisResult(key, network.group_name, atoms=_atoms_from_system(network.sys),
            layers=_layers_from_network(network), export_group=owner,
            complete_cells=int(network.balls.get('complete', []).sum()) if 'complete' in network.balls else 0,
            surface_count=len(network.surfs) if network.surfs is not None else 0,
            info_sections=cached_info_sections(network))
        results[key] = result
        if network is session.active_network:
            active_source = key
    return project, results, active_source or next(iter(results), None), {}


def cached_info_sections(network):
    """Expose archived analysis values directly; never write/read an info file."""
    sections = {'Voronoi Network': [(name, str(len(getattr(network, name))))
        for name in ('balls', 'verts', 'edges', 'surfs') if getattr(network, name, None) is not None]}
    for iface in network.sys.ifaces or []:
        if iface.net is not network:
            continue
        analysis = getattr(iface, 'geometry_analysis', None)
        if analysis is None:
            sections[f'Interface: {iface.name}'] = [('Analysis', 'Not calculated')]
            continue
        sections.setdefault('Group Geometry', []).extend(
            (f'{iface.name}: {key}', str(value)) for key, value in analysis.voronoi_side_1.items()
            if key.startswith('topology_') or key in ('area', 'euler_characteristic', 'boundary_loops', 'connected_components'))
        sections.setdefault('Surface Curvature', []).extend(
            (f'{iface.name}: {key}', str(value)) for key, value in analysis.voronoi_side_1.items())
        sections['Surface Curvature'].extend((f'{iface.name}: coverage {key}', str(value))
            for key, value in analysis.coverage.items())
    return sections
