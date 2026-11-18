"""Molecule-level water shells, solved against the complete input geometry."""
import math
import gc
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from vorpy.src.boundary import WATER_RESIDUES, CellTarget, SolventMode
from vorpy.src.group.group import Group
from vorpy.src.group.layers import _build_index_maps
from vorpy.src.output.layout import unique_output_directory


def complete_layer_count(group):
    """Count consecutive nonempty layers whose water cells are all complete."""
    mapping = group.net.balls.get('system_num', group.net.balls['num'])
    complete = {int(atom) for atom, done in zip(mapping, group.net.balls['complete'].eq(True).fillna(False)) if done}
    for index, atoms in enumerate(group.solvent_layers['atoms']):
        if not atoms or not set(atoms).issubset(complete):
            return index
    return len(group.solvent_layers['atoms'])


def water_neighbors(network, system, frontier, visited):
    """Return complete water residues sharing a positive-area face with frontier."""
    _, mapping, _ = _build_index_maps(SimpleNamespace(net=network))
    frontier = set(frontier)
    neighbors = set()
    for face in network.surfs.itertuples():
        if not math.isfinite(face.sa) or face.sa <= 0:
            continue
        first, second = (mapping.get(int(atom), -1) for atom in face.balls)
        for inside, outside in ((first, second), (second, first)):
            if inside not in frontier or outside in visited or not 0 <= outside < len(system.balls):
                continue
            residue = system.balls.iloc[outside]['res']
            if residue.name.strip().upper() in WATER_RESIDUES:
                neighbors.update(int(atom) for atom in residue.atoms)
    return neighbors - visited


def boundary_indices(network, seed_atoms, layers):
    """Select only faces between layer n-1 and n, plus their edges/vertices."""
    _, mapping, _ = _build_index_maps(SimpleNamespace(net=network))
    membership = {int(atom): 0 for atom in seed_atoms}
    membership.update({int(atom): layer for layer, atoms in enumerate(layers, 1) for atom in atoms})
    surfaces = [[] for _ in layers]
    for face in network.surfs.itertuples():
        if not math.isfinite(face.sa) or face.sa <= 0:
            continue
        first, second = (membership.get(mapping.get(int(atom))) for atom in face.balls)
        if first is not None and second is not None and abs(first - second) == 1:
            surfaces[max(first, second) - 1].append(int(face.Index))
    edges, vertices = [], []
    for selected in surfaces:
        edge_ids = {int(edge) for face in selected for edge in network.surfs.loc[face, 'edges']}
        vertex_ids = {int(vertex) for face in selected for vertex in network.surfs.loc[face, 'verts']}
        vertex_ids.update(int(vertex) for edge in edge_ids for vertex in network.edges.loc[edge, 'verts'])
        edges.append(sorted(edge_ids))
        vertices.append(sorted(vertex_ids))
    return surfaces, edges, vertices


def build_solvent_layers(system, count, *, calculate_curvature=True, seed_groups=None):
    """Discover each frontier once, then solve one combined solvent export network."""
    original_groups = list(system.groups)
    groups = list(seed_groups) if seed_groups is not None else [
        group for group in original_groups if not getattr(group, 'solvent_layers', None)]
    seed_atoms = sorted({int(atom) for group in groups for atom in group.ball_ndxs
                         if system.balls.iloc[int(atom)]['res'].name.strip().upper() not in WATER_RESIDUES})
    if not seed_atoms:
        raise ValueError('Water layers require a non-water solute selection.')
    saved_layers = [group for group in original_groups if getattr(group, 'solvent_layers', None)]
    for saved in saved_layers:
        if saved.solvent_layers['seed_atoms'] == seed_atoms and len(saved.solvent_layers['atoms']) == count:
            print(f'Reusing {count} saved water layers; no geometry solve is needed.', flush=True)
            system.groups = [group for group in original_groups if group not in saved_layers] + [saved]
            return saved
    settings = dict(groups[0].settings)
    for key in ('generator_policy', 'source_generator_count', 'active_generator_count',
                'target_cell_count', 'target_system_indices', 'excluded_generator_count',
                'excluded_water_generator_count'):
        settings.pop(key, None)
    settings['solvent_mode_resolution'] = 'explicit_layer_context'
    settings['solvent_mode'] = SolventMode.SOLVENT_CONTEXT.value
    settings['cell_target'] = CellTarget.MOLECULAR.value
    directory = unique_output_directory(Path(system.files['dir']), 'solvent',
        [group.dir for group in groups if group.dir])
    directory.mkdir(parents=True)
    # Frontier vertex checkpoints must not overwrite the solute's vertex file.
    settings['sys_dir'] = str(directory)
    visited = set(seed_atoms)
    frontier = set(seed_atoms)
    cached = []
    physical_atoms = set(range(len(system.balls)))
    for group in original_groups:
        network = group.net or getattr(group, 'source_network', None)
        if network is None or network.settings['net_type'] != settings['net_type']:
            continue
        _, mapping, _ = _build_index_maps(SimpleNamespace(net=network))
        if physical_atoms.issubset(mapping.values()):
            targets = {mapping[int(atom)] for atom in network.group if int(atom) in mapping}
            if all(network is not previous for previous, _ in cached):
                cached.append((network, targets))

    def cached_for(atoms, *, require_analysis=False):
        fields = {'vol', 'vdw_vol', 'com', 'moi', 'int_mean_curv', 'int_gauss_curv'}
        return next((network for network, targets in cached if set(atoms).issubset(targets)
                     and (not require_analysis or fields.issubset(network.balls.columns))), None)

    covered = set().union(*(targets for _, targets in cached))
    if not set(seed_atoms).issubset(covered):
        missing = sorted(set(seed_atoms) - covered)
        print(f'Saved geometry lacks wet-context seed cells; solving {len(missing)} seed cells '
              'in solvent context. The saved network is preserved.', flush=True)
        probe = Group(system, name='water_seed_context', atoms=missing,
                      settings=dict(settings), print_metrics=False)
        probe.build(calculate_curvature=False)
        cached.append((probe.net, set(missing)))
    networks = [network for network, _ in cached]
    layers = []
    water_atoms = {i for i, atom in system.balls.iterrows()
                   if atom['res'].name.strip().upper() in WATER_RESIDUES}
    for layer in range(1, count + 1):
        following = set()
        for network in networks:
            following.update(water_neighbors(network, system, frontier, visited))
        layers.append(sorted(following))
        visited.update(following)
        if layer == count or not following or water_atoms.issubset(visited):
            layers.extend([] for _ in range(count - layer))
            break
        reused = cached_for(following)
        if reused is not None:
            networks = [reused]
            frontier = following
            continue
        print(f'Solving water frontier {layer}: {len(following)} atoms', flush=True)
        probe = Group(system, name=f'water_frontier_{layer}', atoms=sorted(following),
                      settings=dict(settings), print_metrics=False)
        probe.build(calculate_curvature=False)
        frontier = following
        networks = [probe.net]
        cached.append((probe.net, set(following)))

    selected = sorted({atom for layer in layers for atom in layer})
    reused = cached_for(selected, require_analysis=calculate_curvature) if selected else networks[0]
    name = 'solvent'
    identifiers = {group.group_id or group.name for group in original_groups}
    suffix = 2
    while name in identifiers:
        name = f'solvent_{suffix}'
        suffix += 1
    solvent = Group(system, name=name, atoms=selected or list(seed_atoms),
                    settings=settings, make_net=False, print_metrics=False)
    if reused is not None:
        solvent.net = reused
    else:
        solvent.build(calculate_curvature=calculate_curvature)
    if not selected:
        solvent.ball_ndxs = []
    surfaces, edges, vertices = boundary_indices(solvent.net, seed_atoms, layers)
    solvent.solvent_layers = dict(atoms=layers, surfs=surfaces, edges=edges, verts=vertices,
                                  seed_atoms=seed_atoms)
    solvent.dir = str(directory)
    system.groups = [group for group in original_groups if group not in saved_layers] + [solvent]
    return solvent


def stream_solvent_layers(system, count, export_layer, *, calculate_curvature=True,
                          seed_groups=None):
    """Solve, export, and release one water shell at a time.

    All input generators remain available as geometric competitors.  Only the
    current shell is a target, so its topology and analysis can be released
    after ``export_layer`` returns.  Volumes are in A^3, areas in A^2, and
    integrated curvatures use the network's existing A / dimensionless units.
    """
    original_groups = list(system.groups)
    groups = list(seed_groups) if seed_groups is not None else [
        group for group in original_groups if not getattr(group, 'solvent_layers', None)]
    seed_atoms = sorted({int(atom) for group in groups for atom in group.ball_ndxs
                         if system.balls.iloc[int(atom)]['res'].name.strip().upper()
                         not in WATER_RESIDUES})
    if not seed_atoms:
        raise ValueError('Water layers require a non-water solute selection.')

    settings = dict(groups[0].settings)
    for key in ('generator_policy', 'source_generator_count', 'active_generator_count',
                'target_cell_count', 'target_system_indices', 'excluded_generator_count',
                'excluded_water_generator_count'):
        settings.pop(key, None)
    settings.update(solvent_mode_resolution='explicit_layer_context',
                    solvent_mode=SolventMode.SOLVENT_CONTEXT.value,
                    cell_target=CellTarget.MOLECULAR.value)
    directory = unique_output_directory(Path(system.files['dir']), 'solvent',
        [group.dir for group in groups if group.dir])
    directory.mkdir(parents=True)
    settings['sys_dir'] = str(directory)

    physical_atoms = set(range(len(system.balls)))
    source = None
    source_targets = set()
    for group in original_groups:
        network = group.net or getattr(group, 'source_network', None)
        if network is None or network.settings['net_type'] != settings['net_type']:
            continue
        _, mapping, _ = _build_index_maps(SimpleNamespace(net=network))
        if not physical_atoms.issubset(mapping.values()):
            continue
        targets = {mapping[int(atom)] for atom in network.group if int(atom) in mapping}
        if set(seed_atoms).issubset(targets):
            source, source_targets = network, targets
            break
    if source is None:
        print(f'Solving wet-context seed cells: {len(seed_atoms)} atoms', flush=True)
        seed = Group(system, name='water_seed_context', atoms=list(seed_atoms),
                     settings=dict(settings), print_metrics=False)
        seed.build(calculate_curvature=False)
        source, source_targets = seed.net, set(seed_atoms)
        del seed

    # Main outputs/checkpoint already exist; retain only the source needed to find layer 1.
    for group in original_groups:
        group.net = None
        group.source_network = None
        group.source_view = None
        group.source_view_analysis = None

    def partial_volume(network, atoms):
        mapping = network.balls.get('system_num', network.balls['num'])
        rows = {int(system_id): row for system_id, (_, row)
                in zip(mapping, network.balls.iterrows())}
        return sum(float(rows[atom]['vol']) for atom in atoms
                   if atom in rows and bool(rows[atom].get('complete', False))
                   and math.isfinite(float(rows[atom].get('vol', math.nan))))

    contained_volume = partial_volume(source, seed_atoms)
    visited, frontier = set(seed_atoms), set(seed_atoms)
    water_atoms = {i for i, atom in system.balls.iterrows()
                   if atom['res'].name.strip().upper() in WATER_RESIDUES}
    summaries = []
    analysis_fields = {'vol', 'vdw_vol', 'com', 'moi', 'int_mean_curv', 'int_gauss_curv'}
    for number in range(1, count + 1):
        atoms = sorted(water_neighbors(source, system, frontier, visited)) if frontier else []
        if atoms:
            can_reuse = set(atoms).issubset(source_targets)
            can_reuse &= not calculate_curvature or analysis_fields.issubset(source.balls.columns)
            # Group normalizes its atom input in place; keep frontier membership immutable.
            layer = Group(system, name=f'layer_{number:02d}', atoms=list(atoms),
                          settings=dict(settings), make_net=False, print_metrics=False)
            if can_reuse:
                layer.net = source
            else:
                source = None
                gc.collect()
                # shortcut: reuse prior outer vertices once Network accepts them as bounded seeds.
                print(f'Solving water layer {number}: {len(atoms)} atoms', flush=True)
                layer.build(calculate_curvature=calculate_curvature)
            surfaces, edges, vertices = boundary_indices(layer.net, frontier, [atoms])
        else:
            layer = Group(system, name=f'layer_{number:02d}', atoms=list(seed_atoms),
                          settings=dict(settings), make_net=False, print_metrics=False)
            layer.ball_ndxs = []
            layer.net = source
            surfaces, edges, vertices = [[]], [[]], [[]]

        layer.solvent_layers = dict(atoms=[atoms], surfs=surfaces, edges=edges,
                                    verts=vertices, seed_atoms=sorted(frontier),
                                    numbers=[number], allow_partial=True)
        layer.dir = str(directory / f'layer_{number:02d}')
        layer.get_info()
        contained_volume += float(layer.vol)
        layer.contained_volume = contained_volume
        layer.solvent_layer_number = number
        mapping = layer.net.balls.get('system_num', layer.net.balls['num'])
        complete = {int(atom) for atom, done in zip(
            mapping, layer.net.balls['complete'].eq(True).fillna(False)) if done}
        layer.complete_cell_count = sum(atom in complete for atom in atoms)
        layer.total_cell_count = len(atoms)
        layer.layer_complete = bool(atoms) and layer.complete_cell_count == len(atoms)
        layer.metrics_status = 'complete' if layer.layer_complete else 'partial'
        export_layer(layer)
        summaries.append(dict(layer=number, atoms=len(atoms),
            complete_cells=layer.complete_cell_count, complete=layer.layer_complete,
            metrics_status=layer.metrics_status, contained_volume=layer.contained_volume,
            layer_volume=layer.vol, layer_surface_area=layer.sa,
            integrated_mean_curvature=layer.int_mean_curv,
            integrated_gaussian_curvature=layer.int_gauss_curv))

        visited.update(atoms)
        frontier = set(atoms)
        if atoms:
            source = layer.net
            _, mapping, _ = _build_index_maps(SimpleNamespace(net=source))
            source_targets = {mapping[int(atom)] for atom in source.group if int(atom) in mapping}
        if water_atoms.issubset(visited):
            frontier = set()
        del layer

    pd.DataFrame(summaries).to_csv(directory / 'layers.csv', index=False)
    (directory / 'README.txt').write_text(
        'Each layer_NN directory contains one independently exported water-shell solve.\n'
        'layer_volume and contained_volume are in A^3; layer_surface_area is in A^2.\n'
        'contained_volume is the solute volume plus every shell through the listed layer.\n'
        'Integrated mean curvature uses A; integrated Gaussian curvature is dimensionless.\n'
        'metrics_status=partial means totals use available complete cells and boundary geometry.\n'
        'All input atoms remain geometric competitors, but only one layer topology is retained at a time.\n',
        encoding='utf-8')
    return summaries
