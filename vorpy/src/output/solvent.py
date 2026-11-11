"""The normal export presets applied to disjoint molecular water layers."""
from pathlib import Path

import pandas as pd

from vorpy.src.output.edges import write_edges
from vorpy.src.output.mesh import MeshData, write_mesh
from vorpy.src.output.pdb import write_pdb
from vorpy.src.output.surfs import write_surfs
from vorpy.src.output.verts import write_off_verts
from vorpy.src.group.solvent_layers import complete_layer_count


def export_solvent_layers(group, plan):
    data = group.solvent_layers
    directory = Path(group.dir)
    directory.mkdir(parents=True, exist_ok=True)
    network = group.net
    complete_count = complete_layer_count(group)
    requested = len(data['atoms'])
    if complete_count < requested:
        excluded = str(requested) if complete_count + 1 == requested else f'{complete_count + 1}-{requested}'
        print(f'\nWarning: {requested} layers were requested; only {complete_count} complete '
              f'{"layer exists" if complete_count == 1 else "layers exist"}. '
              f'{"Layer" if complete_count + 1 == requested else "Layers"} {excluded} '
              'will not be included in logs.', flush=True)
    mapping = network.balls.get('system_num', network.balls['num'])
    cells = {int(system_id): row for system_id, (_, row) in zip(mapping, network.balls.iterrows())}
    complete_atoms = {int(atom) for atom, done in zip(mapping, network.balls['complete'].eq(True).fillna(False)) if done}
    summaries = []
    for layer, atoms in enumerate(data['atoms'], 1):
        molecules = {}
        for atom in atoms:
            residue = group.sys.balls.iloc[atom]['res']
            molecules.setdefault(min(residue.atoms), []).append(atom)
        summaries.append(dict(layer=layer, atoms=len(atoms), water_molecules=len(molecules),
            complete=bool(atoms) and set(atoms).issubset(complete_atoms),
            included_in_logs=layer <= complete_count,
            complete_cells=sum(atom in complete_atoms for atom in atoms),
            complete_molecules=sum(all(atom in complete_atoms for atom in members)
                                   for members in molecules.values()),
            boundary_surfaces=len(data['surfs'][layer - 1]),
            boundary_edges=len(data['edges'][layer - 1])))
    pd.DataFrame(summaries).to_csv(directory / 'layers.csv', index=False)
    (directory / 'README.txt').write_text(
        'Layer 1 touches the selected solute; later layers follow water-water Voronoi faces.\n'
        'Each water molecule appears once, at its shortest observed face distance. Ions are excluded.\n'
        'layer_NN_surfs/edges/verts show the interface from layer NN-1 to NN (solute = 0).\n'
        'The single native *_logs.csv includes only consecutive complete layers, starting at layer 1.\n'
        'A layer is complete only if it is nonempty and every water atom cell is marked complete.\n'
        'The first incomplete/empty layer and all later layers are excluded from log rows and totals.\n'
        'Layer geometry and checkpoint membership remain available for inspection.\n'
        'Missing faces can affect deeper layer assignments. Empty layers have explicit empty files.\n',
        encoding='utf-8')
    exported = set()
    for _, options in plan:
        option = next(iter(options))
        if option in {'surr_atoms', 'surr_resids', 'ext_atoms'}:
            continue
        kind = next((kind for kind in ('surfs', 'edges', 'verts') if kind in option), None)
        if option == 'atoms' or (kind and not option.startswith('atom_')):
            kind = 'pdb' if option == 'atoms' else kind
            if kind in exported:
                continue
            exported.add(kind)
            for index, atoms in enumerate(data['atoms']):
                name = f'layer_{index + 1:02d}'
                if kind == 'pdb':
                    if atoms:
                        write_pdb(atoms, name, group.sys, directory=str(directory))
                    else:
                        (directory / (name + '.pdb')).write_text('HEADER  empty water layer\nEND\n')
                    continue
                selected = data[kind][index]
                name += '_' + kind
                file_type = group.settings.get('file_type', 'off')
                if not selected:
                    write_mesh(MeshData([], []), name, file_type, str(directory))
                    continue
                writer = {'surfs': write_surfs, 'edges': write_edges, 'verts': write_off_verts}[kind]
                inside = set(data['seed_atoms'])
                inside.update(atom for layer in data['atoms'][:index] for atom in layer)
                kwargs = dict(directory=str(directory), file_type=file_type,
                    target_cells=[int(cells[atom]['num']) for atom in sorted(inside)])
                if kind == 'edges':
                    kwargs['radius'] = group.settings['edge_width']
                    if group.settings.get('_edge_col_explicit'):
                        kwargs.update(color=group.settings['edge_col'], color_scheme='fixed')
                elif kind == 'verts':
                    kwargs['vert_rad'] = group.settings['vertex_size']
                    if group.settings.get('_vert_col_explicit'):
                        kwargs.update(color=group.settings['vert_col'], color_scheme='fixed')
                writer(network, selected, name, **kwargs)
        else:
            group.exports(**options)
