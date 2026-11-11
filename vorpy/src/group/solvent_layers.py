"""Molecule-level water shells, solved against the complete input geometry."""
import math
from pathlib import Path
from types import SimpleNamespace

from vorpy.src.boundary import WATER_RESIDUES, CellTarget
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
        first, second = (mapping[int(atom)] for atom in face.balls)
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
        first, second = (membership.get(mapping[int(atom)]) for atom in face.balls)
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


def build_solvent_layers(system, count, *, calculate_curvature=True):
    """Discover each frontier once, then solve one combined solvent export network."""
    groups = list(system.groups)
    seed_atoms = sorted({int(atom) for group in groups for atom in group.ball_ndxs
                         if system.balls.iloc[int(atom)]['res'].name.strip().upper() not in WATER_RESIDUES})
    if not seed_atoms:
        raise ValueError('Water layers require a non-water solute selection.')
    settings = dict(groups[0].settings)
    settings['cell_target'] = CellTarget.MOLECULAR.value
    directory = unique_output_directory(Path(system.files['dir']), 'solvent',
        [group.dir for group in groups if group.dir])
    directory.mkdir(parents=True)
    # Frontier vertex checkpoints must not overwrite the solute's vertex file.
    settings['sys_dir'] = str(directory)
    visited = set(seed_atoms)
    frontier = set(seed_atoms)
    networks = [group.net or group.source_network for group in groups]
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
        print(f'Solving water frontier {layer}: {len(following)} atoms', flush=True)
        probe = Group(system, name=f'water_frontier_{layer}', atoms=sorted(following),
                      settings=dict(settings), print_metrics=False)
        probe.build(calculate_curvature=False)
        frontier = following
        networks = [probe.net]

    selected = sorted({atom for layer in layers for atom in layer})
    solvent = Group(system, name='solvent', atoms=selected or list(seed_atoms),
                    settings=settings, make_net=bool(selected), print_metrics=False)
    if selected:
        solvent.build(calculate_curvature=calculate_curvature)
    else:
        # An empty cohort still has four explicit empty exports and a readable checkpoint.
        from copy import copy
        solvent.ball_ndxs = []
        solvent.net = copy(networks[0])
        solvent.net.group = []
    surfaces, edges, vertices = boundary_indices(solvent.net, seed_atoms, layers)
    solvent.solvent_layers = dict(atoms=layers, surfs=surfaces, edges=edges, verts=vertices,
                                  seed_atoms=seed_atoms)
    solvent.dir = str(directory)
    system.groups.append(solvent)
    return solvent
