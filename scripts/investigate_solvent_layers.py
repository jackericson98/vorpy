"""Audit water-only Voronoi layers in a solved archive, without resolving geometry."""
import argparse
from collections import deque
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from vorpy.src.boundary import WATER_RESIDUES
from vorpy.src.group import Group
from vorpy.src.group.layers import _build_index_maps
from vorpy.src.io import load_network
from vorpy.src.output.pdb import write_pdb


def shortest_layers(adjacency, seeds, water_ids, depth):
    """A water belongs only to its shortest face-connected distance from the seed."""
    distances = {seed: 0 for seed in seeds}
    queue = deque(seeds)
    while queue:
        current = queue.popleft()
        if distances[current] == depth:
            continue
        for neighbor in adjacency.get(current, ()):
            if neighbor in water_ids and neighbor not in distances:
                distances[neighbor] = distances[current] + 1
                queue.append(neighbor)
    return {water: layer for water, layer in distances.items() if water in water_ids}


def self_check():
    # Shortest paths stay disjoint through cycles and cannot bridge through an ion.
    graph = {0: {1, 2, 9}, 1: {0, 2}, 2: {0, 1, 3}, 3: {2, 4}, 4: {3}, 9: {0, 5}, 5: {9}}
    assert shortest_layers(graph, {0}, {1, 2, 3, 4, 5}, 4) == {1: 1, 2: 1, 3: 2, 4: 3}
    assert shortest_layers(graph, {0}, {1, 2, 3, 4, 5}, 1) == {1: 1, 2: 1}


def load_archive_tables(path, structure):
    """Read only numeric geometry tables from an older archive, leaving its graph untouched."""
    import zipfile
    from vorpy.src.io.archive_codec import Reader
    from vorpy.src.system import System
    system = System(str(structure.resolve()))
    with zipfile.ZipFile(path) as archive:
        reader = Reader(archive)
        metadata = reader.json('metadata.json')
        if metadata.get('format') != 'VorPy Network Archive' or metadata.get('format_version') != 3:
            raise ValueError('Table audit supports only version 3 network archives')
        objects = reader.json('state.json')['objects']
        fields = next(obj['fields'] for obj in objects
                      if obj['kind'] == 'network' and obj['id'] == metadata['root_network'])
        tables = {name: reader.decode(fields[name]) for name in ('balls', 'surfs', 'edges', 'verts')}
        settings = reader.decode(fields['settings'])
    mapping = tables['balls'].get('system_num', tables['balls']['num']).astype(int)
    if set(mapping) != set(system.balls.index) or len(mapping) != len(system.balls):
        raise ValueError('Table audit requires all PDB atoms to be represented exactly once')
    if not np.allclose(np.stack(tables['balls']['loc']),
                       np.stack(system.balls.loc[mapping, 'loc']), atol=.002, rtol=0):
        raise ValueError('Archived coordinates do not match the supplied PDB')
    return SimpleNamespace(sys=system, settings=settings, **tables)


def audit(network, seed_name, depth, destination):
    system = network.sys
    atoms = system.balls
    _, topology_to_system, _ = _build_index_maps(SimpleNamespace(net=network))
    # Residue object membership distinguishes repeated/wrapped PDB residue numbers.
    residue_ids = {i: min(atom.res.atoms) for i, atom in atoms.iterrows()}
    water_ids = {residue_ids[i] for i, atom in atoms.iterrows()
                 if atom.res.name.strip().upper() in WATER_RESIDUES}
    seed_atoms = [i for i, atom in atoms.iterrows() if atom.res.name.strip().upper() == seed_name]
    if not seed_atoms or seed_name in WATER_RESIDUES:
        raise ValueError('Select an existing non-water seed residue name')
    if len(set(topology_to_system.values())) != len(topology_to_system):
        raise ValueError('Ambiguous topology-to-system mapping')
    mapped_atoms = set(topology_to_system.values())
    missing = [i for i in atoms.index if residue_ids[i] in water_ids and i not in mapped_atoms]
    if missing:
        raise ValueError(f'Archive omits {len(missing)} water atoms; rebuild with solvent included')
    if not set(seed_atoms).issubset(mapped_atoms):
        raise ValueError('Archive omits seed atoms')

    adjacency = {}
    faces = []
    skipped = 0
    for row_id, face in network.surfs.iterrows():
        a, b = map(int, face.balls)
        if a not in topology_to_system or b not in topology_to_system:
            raise ValueError('Surface endpoint is absent from the network mapping')
        a, b = topology_to_system[a], topology_to_system[b]
        if a not in residue_ids or b not in residue_ids:
            continue  # Virtual boundary generators are not water molecules.
        if not np.isfinite(face.sa) or face.sa <= 0:
            skipped += 1
            continue
        ra, rb = residue_ids[a], residue_ids[b]
        adjacency.setdefault(ra, set()).add(rb)
        adjacency.setdefault(rb, set()).add(ra)
        faces.append(dict(surface_index=row_id, atom_i=a, atom_j=b,
                          molecule_i=ra, molecule_j=rb, surface_area=face.sa))
    seeds = {residue_ids[i] for i in seed_atoms}
    distances = shortest_layers(adjacency, seeds, water_ids, depth)
    # Independent set-frontier traversal verifies the queue result on the real graph.
    seen, frontier = set(seeds), set(seeds)
    for layer in range(1, depth + 1):
        frontier = {b for a in frontier for b in adjacency.get(a, ()) if b in water_ids} - seen
        assert frontier == {water for water, distance in distances.items() if distance == layer}
        seen.update(frontier)

    destination.mkdir(parents=True, exist_ok=False)
    solvent = destination / 'solvent'
    solvent.mkdir()
    # Saved trajectory archives may reference an expired temporary PDB; use embedded atoms.
    if not system.files.get('base_file') or not Path(system.files['base_file']).is_file():
        system.files['base_file'] = None
    write_pdb(seed_atoms, 'seed', system, directory=str(destination))
    face_table = pd.DataFrame(faces, columns=['surface_index', 'atom_i', 'atom_j',
                                             'molecule_i', 'molecule_j', 'surface_area'])
    topology_rows = {topology_to_system[int(row.num)]: row for _, row in network.balls.iterrows()}
    summaries = []
    viewer = [f'load "{(destination / "seed.pdb").as_posix()}", seed', 'hide everything',
              'show sticks, seed', 'color gray70, seed']
    for layer in range(1, depth + 1):
        members = {water for water, distance in distances.items() if distance == layer}
        selected = [i for i in atoms.index if residue_ids[i] in members]
        folder = solvent / f'layer_{layer:02d}'
        folder.mkdir()
        write_pdb(selected, 'water', system, directory=str(folder))
        if selected:
            assert sum(line.startswith(('ATOM  ', 'HETATM')) for line in
                       (folder / 'water.pdb').read_text().splitlines()) == len(selected)
            viewer += [f'load "{(folder / "water.pdb").as_posix()}", water_layer_{layer:02d}',
                       f'show lines, water_layer_{layer:02d}',
                       f'color {("cyan", "green", "yellow", "salmon")[(layer - 1) % 4]}, water_layer_{layer:02d}']
        records = []
        for i in selected:
            atom, cell = atoms.loc[i], topology_rows[i]
            records.append(dict(system_index=i, molecule_id=residue_ids[i], layer=layer,
                name=atom['name'], residue=atom.res.name, residue_sequence=atom.res.seq,
                element=atom.element, x=atom['loc'][0], y=atom['loc'][1], z=atom['loc'][2],
                volume=cell.vol, surface_area=cell.sa, complete=bool(cell.complete)))
        pd.DataFrame(records, columns=['system_index', 'molecule_id', 'layer', 'name', 'residue',
            'residue_sequence', 'element', 'x', 'y', 'z', 'volume', 'surface_area', 'complete']).to_csv(
            folder / 'atoms.csv', index=False)
        face_table[face_table.atom_i.isin(selected) | face_table.atom_j.isin(selected)].to_csv(
            folder / 'incident_surfaces.csv', index=False)
        summaries.append(dict(layer=layer, water_molecules=len(members), atoms=len(selected),
            complete_cells=sum(bool(topology_rows[i].complete) for i in selected),
            complete_molecules=sum(all(bool(topology_rows[i].complete) for i in selected
                                      if residue_ids[i] == water) for water in members)))
    pd.DataFrame(summaries).to_csv(solvent / 'summary.csv', index=False)
    (destination / 'view_layers.pml').write_text('\n'.join(viewer + ['group solvent, water_layer_*',
        'bg_color white', 'zoom all']), encoding='utf-8')

    legacy = Group(system, name='legacy_seed', atoms=list(seed_atoms), settings=dict(network.settings),
                   make_net=False)
    legacy.net = network
    legacy.get_layers(max_layers=depth)
    legacy_rows = []
    seen = set(seed_atoms)
    for layer in range(1, depth + 1):
        selected = set(legacy.layer_atoms[layer]) if layer < len(legacy.layer_atoms) else set()
        expected = {i for i in atoms.index if distances.get(residue_ids[i]) == layer}
        legacy_rows.append(dict(layer=layer, legacy_atoms=len(selected), expected_atoms=len(expected),
            repeated_atoms=len(selected & seen), non_water_atoms=sum(residue_ids[i] not in water_ids for i in selected),
            missing_expected_atoms=len(expected - selected), unexpected_atoms=len(selected - expected)))
        seen.update(selected)
    pd.DataFrame(legacy_rows).to_csv(destination / 'legacy_comparison.csv', index=False)
    complete_seed = sum(bool(topology_rows[i].complete) for i in seed_atoms)
    complete_inner = all(row['complete_cells'] == row['atoms'] for row in summaries[:-1])
    report = dict(seed_residue=seed_name, seed_atoms=len(seed_atoms), complete_seed_cells=complete_seed,
        status='observed_graph_only' if complete_seed != len(seed_atoms) or not complete_inner else 'inner_cells_complete',
        total_atoms=len(atoms),
        network_atoms=len(network.balls), total_water_molecules=len(water_ids),
        assigned_water_molecules=len(distances), unassigned_water_molecules=len(water_ids - distances.keys()),
        excluded_nonpositive_or_nonfinite_faces=skipped, layers=summaries, legacy=legacy_rows)
    (destination / 'report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--network', type=Path)
    source.add_argument('--archive-tables', type=Path, help='Read-only numeric-table audit of an older v3 archive')
    parser.add_argument('--structure', type=Path, help='Matching PDB required with --archive-tables')
    parser.add_argument('--seed-residue', required=True)
    parser.add_argument('--layers', type=int, default=4)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.layers < 1:
        parser.error('--layers must be positive')
    if args.archive_tables and not args.structure:
        parser.error('--archive-tables requires --structure')
    self_check()
    loader = None
    try:
        if args.archive_tables:
            network = load_archive_tables(args.archive_tables, args.structure)
        else:
            from vorpy.src.command.vpy_cmnd import Command
            if Command._is_canonical_checkpoint(args.network):
                loader = Command()
                network = loader._load_canonical_checkpoint(args.network)
            else:
                network = load_network(args.network)
        audit(network, args.seed_residue.strip().upper(), args.layers, args.output.resolve())
    finally:
        if loader is not None:
            loader._cleanup_checkpoint_source()


if __name__ == '__main__':
    main()
