"""Read-only, selective checkpoint diagnostics for the six PEG first-layer holds."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import csv
import io
import json
from pathlib import Path
import sys
from types import SimpleNamespace
import zipfile

import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull, cKDTree

REPOSITORY_ROOT = next(parent for parent in Path(__file__).resolve().parents
                       if (parent / 'pyproject.toml').is_file())
sys.path.insert(0, str(REPOSITORY_ROOT))
from vorpy.src.io.archive_codec import Reader
from vorpy.src.inputs.frames import is_pdb_virtual_site
from vorpy.src.network.analyze import diagnose_geometry_completeness


CASES = [('aw', 'PEG2000', 9), ('aw', 'PEG4000', 2), ('aw', 'PEG4000', 3),
         ('aw', 'PEG4000', 7), ('aw', 'PEG4000', 11), ('prm', 'PEG4000', 9)]


@contextmanager
def geometry_archive(path, owner='solvent'):
    with zipfile.ZipFile(path) as outer:
        if 'checkpoint/manifest.json' in outer.namelist():
            manifest = json.loads(outer.read('checkpoint/manifest.json'))
            state = next(n for n in manifest['document']['networks']
                         if n['network_id'].endswith(':' + owner))
            payload = next(p for p in manifest['payloads'] if p['payload_id'] == state['payload_ref'])
            with zipfile.ZipFile(io.BytesIO(outer.read(payload['member']))) as inner:
                yield inner
        else:
            yield outer


def table(reader, path, columns):
    node = reader.json(path)
    return pd.DataFrame({c: reader.column(v, node['count']) for c, v in node['columns'].items()
                         if c in columns}, index=reader.decode(node['index']))


def read_geometry(archive):
    reader = Reader(archive)
    meta = reader.json('metadata.json')
    records = reader.json('state.json')['objects']
    root = records[meta['root_network']]['fields']
    fields = {
        'balls': {'num', 'system_num', 'loc', 'rad', 'complete', 'surfs', 'edges', 'verts'},
        'verts': {'balls', 'edges', 'surfs', 'loc', 'rad', 'dub'},
        'edges': {'balls', 'verts', 'surfs'}, 'surfs': {'balls', 'edges', 'verts', 'sa'},
    }
    net = SimpleNamespace(boundary_indices=(), group=())
    for kind, columns in fields.items():
        setattr(net, kind, table(reader, root[kind]['table'], columns))
    net.settings = {'net_type': meta['scheme']}
    return net, meta


def source_frame(path, ordinal):
    atoms, cryst, current = [], None, 1
    with path.open() as stream:
        for line in stream:
            if line.startswith('CRYST1'):
                cryst = line.rstrip()
            if line.startswith(('ATOM  ', 'HETATM')) and not is_pdb_virtual_site(line):
                atoms.append(line)
            if line[:6].strip() in {'END', 'ENDMDL'} and atoms:
                if current == ordinal:
                    return atoms, cryst
                atoms = []
                current += 1
    if atoms and current == ordinal:
        return atoms, cryst
    raise ValueError(f'Missing frame {ordinal}: {path}')


def atom_identity(line):
    return (line[12:16].strip(), line[17:20].strip(), line[21:22], int(line[22:26]),
            *(float(line[a:b]) for a, b in [(30, 38), (38, 46), (46, 54)]))


def layer_membership(source, layer):
    lookup = {atom_identity(line): index for index, line in enumerate(source)}
    if len(lookup) != len(source):
        raise ValueError('Nonunique source identities/coordinates need disambiguation.')
    selected = []
    for line in layer:
        index = lookup[atom_identity(line)]
        if line[6:11] != source[index][6:11]:
            raise ValueError('Layer PDB serial disagrees with its independently matched source atom.')
        selected.append(index)
    if len(selected) != len(set(selected)):
        raise ValueError('Duplicate atom in layer PDB.')
    return selected


def lattice(cryst):
    if cryst is None:
        raise ValueError('CRYST1 is unavailable; no box-edge hypothesis can be tested.')
    a, b, c = [float(cryst[x:y]) for x, y in [(6, 15), (15, 24), (24, 33)]]
    alpha, beta, gamma = np.deg2rad([float(cryst[x:y]) for x, y in [(33, 40), (40, 47), (47, 54)]])
    cx = c * np.cos(beta)
    cy = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    return np.array([[a, 0, 0], [b * np.cos(gamma), b * np.sin(gamma), 0],
                     [cx, cy, np.sqrt(c * c - cx * cx - cy * cy)]])


def investigate(base, output):
    output.mkdir(parents=True, exist_ok=False)
    summaries = []
    for scheme, molecule, frame in CASES:
        label = f'{scheme}_{molecule}_frame_{frame:04d}'
        print(f'Inspecting {label} (no solve)', flush=True)
        root = base / f'Reruns/layers_4/{scheme}/{molecule}/frame_{frame:04d}/{molecule}'
        source, cryst = source_frame(base / f'Frames/{molecule}.pdb', frame)
        layer_lines = [line for line in (root / 'solvent/layer_01.pdb').read_text().splitlines()
                       if line.startswith(('ATOM  ', 'HETATM'))]
        layer = set(layer_membership(source, layer_lines))
        xyz = np.array([atom_identity(line)[4:] for line in source])
        h = lattice(cryst)
        inverse = np.linalg.inv(h)
        fractional = xyz @ inverse
        wrapped = fractional % 1
        face_distances = np.minimum(wrapped, 1 - wrapped) / np.linalg.norm(inverse, axis=0)
        polymer = np.array([i for i, line in enumerate(source) if line[17:20].strip() == 'PEG'])
        hull = ConvexHull(xyz)
        with geometry_archive(root / f'{molecule}_checkpoint.vpy') as archive:
            net, metadata = read_geometry(archive)
        mapping = net.balls.set_index('num')['system_num'].to_dict()
        bad = [int(row.num) for row in net.balls.itertuples()
               if int(row.system_num) in layer and not row.complete]
        net.group = bad
        diagnoses = diagnose_geometry_completeness(net)
        cells = []
        bad_molecules = set()
        for cell in diagnoses['cells']:
            index = mapping[cell['local_id']]
            line = source[index]
            identity = atom_identity(line)
            residue = (identity[1], identity[2], identity[3])
            bad_molecules.add(residue)
            members = [i for i, candidate in enumerate(source) if
                       (candidate[17:20].strip(), candidate[21:22], int(candidate[22:26])) == residue]
            assert np.array_equal(np.asarray(net.balls.loc[cell['local_id'], 'loc']), xyz[index])
            displacement = xyz[polymer] - xyz[index]
            fractional_delta = displacement @ inverse
            fractional_delta -= np.round(fractional_delta)
            # Enumerate neighboring lattice translations after nearest fractional reduction.
            translations = np.array(list(np.ndindex(3, 3, 3))) - 1
            periodic = min(float(np.linalg.norm((fractional_delta + shift) @ h, axis=1).min())
                           for shift in translations)
            cell.update(source_index=int(index), pdb_serial=int(line[6:11]), original_record=line.rstrip(),
                        original_molecule_records=[source[i].rstrip() for i in members],
                        molecule_source_indices=members, fractional_coordinates=fractional[index].tolist(),
                        periodic_face_distance_A=float(face_distances[index].min()),
                        finite_cloud_hull_distance_A=float(-(hull.equations[:, :3] @ xyz[index] + hull.equations[:, 3]).max()),
                        nearest_polymer_finite_A=float(np.linalg.norm(displacement, axis=1).min()),
                        nearest_polymer_periodic_A=periodic)
            cells.append(cell)
        molecules = {(source[i][17:20].strip(), source[i][21:22], int(source[i][22:26])) for i in layer}
        summary = dict(case=label, layer1_atoms=len(layer), incomplete_atoms=len(bad),
                       incomplete_atom_fraction=len(bad) / len(layer), layer1_molecules=len(molecules),
                       affected_molecules=len(bad_molecules), affected_molecule_fraction=len(bad_molecules) / len(molecules),
                       polymer_min_periodic_face_distance_A=float(face_distances[polymer].min()),
                       cryst1=cryst, lattice_rows_A=h.tolist(), checkpoint_metadata=metadata, cells=cells)
        (output / f'{label}.json').write_text(json.dumps(summary, indent=2, default=str), encoding='utf-8')
        summaries.append({k: v for k, v in summary.items() if k not in {'cryst1', 'lattice_rows_A', 'checkpoint_metadata', 'cells'}})
        print(f'  {len(bad)}/{len(layer)} cells; {len(bad_molecules)}/{len(molecules)} waters affected', flush=True)
    with (output / 'incomplete_summary.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    investigate(args.root, args.output)
