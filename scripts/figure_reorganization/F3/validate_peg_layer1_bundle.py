"""Independently validate a recovered layer bundle and its scientific reload."""
import argparse
import csv
import io
import json
from pathlib import Path
import sys
from unittest.mock import patch
import zipfile

import numpy as np

REPOSITORY_ROOT = next(parent for parent in Path(__file__).resolve().parents
                       if (parent / 'pyproject.toml').is_file())
sys.path.insert(0, str(REPOSITORY_ROOT))
if 'figure_reorganization' in Path(__file__).parts:
    from scripts.figure_reorganization.F3.investigate_peg_first_layer import atom_identity, layer_membership, source_frame, table
    from scripts.figure_reorganization.F3.reexport_peg_layer1 import cell_digest, digest
else:
    from scripts.investigate_peg_first_layer import atom_identity, layer_membership, source_frame, table
    from scripts.reexport_peg_layer1 import cell_digest, digest
from vorpy.src.command.vpy_cmnd import Command
from vorpy.src.group.solvent_layers import complete_layer_count
from vorpy.src.io.archive_codec import Reader
from vorpy.src.network import Network
from vorpy.src.output.checkpoint import _load_checkpoint
from types import SimpleNamespace


def raw_fingerprints(archive, roots_only=False):
    reader = Reader(archive)
    metadata = reader.json('metadata.json')
    records = reader.json('state.json')['objects']
    roots = [records[metadata['root_network']]] if roots_only else [r for r in records if r['kind'] == 'network']
    result = {}
    for record in roots:
        fields = record['fields']
        balls = table(reader, fields['balls']['table'], {'loc', 'rad', 'num', 'system_num', 'complete', 'vol', 'sa', 'int_gauss_curv'})
        result[cell_digest(SimpleNamespace(balls=balls))] = {
            'name': fields['group_name'], 'generators': len(balls),
            'targets': len(reader.decode(fields['group'])), 'settings': reader.decode(fields['settings'])}
    return result


def atom_log(path):
    with path.open(newline='', encoding='utf-8') as stream:
        reader = csv.reader(stream)
        for row in reader:
            if row and row[0].lower() == 'group information':
                group = dict(zip(next(reader), next(reader)))
            if row and row[0] == 'Atoms':
                headers = next(reader)
                atoms = []
                for values in reader:
                    if values and values[0] == 'Surfaces':
                        break
                    if values:
                        atoms.append(dict(zip(headers, values)))
                return atoms, group
    raise ValueError(f'No atom section: {path}')


def validate(bundle, original, trajectory, frame):
    bundle, original = bundle.resolve(), original.resolve()
    original_digest = digest(original)
    lines, cryst = source_frame(trajectory, frame)
    layer_path = bundle / 'solvent/layer_01.pdb'
    layer = [line for line in layer_path.read_text().splitlines() if line.startswith(('ATOM  ', 'HETATM'))]
    membership = layer_membership(lines, layer)
    atoms, solvent = atom_log(bundle / 'solvent/aw_logs.csv')
    assert set(membership) == {int(row['Index']) for row in atoms}
    assert len(atoms) == len(membership)
    for row in atoms:
        index = int(row['Index'])
        assert row['Complete Cell?'].lower() == 'true'
        assert row['Name'].strip() == lines[index][12:16].strip()
        assert row['Residue'].strip() == lines[index][17:20].strip()
        assert int(row['Residue Sequence']) == int(lines[index][22:26])
        assert np.array_equal([float(row[k]) for k in ['X', 'Y', 'Z']], atom_identity(lines[index])[4:])
    assert not (bundle / 'solvent/layer_02.pdb').exists()
    with (bundle / 'source_identity_mapping.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['system_index', 'original_pdb_serial', 'atom', 'residue', 'chain', 'residue_sequence', 'x', 'y', 'z'])
        for index in membership:
            writer.writerow([index, int(lines[index][6:11]), *atom_identity(lines[index])])
    checkpoint = next(bundle.glob('*_checkpoint.vpy'))
    print('Validating canonical payload hashes and preserved cell scalar tables...', flush=True)
    with zipfile.ZipFile(original) as archive:
        previous = raw_fingerprints(archive)
    container = _load_checkpoint(checkpoint)
    replayed = {}
    for payload in container.payloads:
        with zipfile.ZipFile(io.BytesIO(payload.content)) as archive:
            replayed.update(raw_fingerprints(archive, roots_only=True))
    assert previous.keys() == replayed.keys(), 'Scientific cell scalar tables changed during replay.'
    evidence = dict(original_archive=str(original), original_archive_sha256=original_digest,
        canonical_checkpoint=str(checkpoint), canonical_checkpoint_sha256=digest(checkpoint),
        source_trajectory=str(trajectory.resolve()), source_frame=frame, source_cryst1=cryst,
        layer1_atoms=len(membership), original_source_serials_verified_by_identity_and_coordinates=True,
        original_cell_fingerprints=previous, replayed_cell_fingerprints=replayed,
        geometry_cell_fingerprints_unchanged=True, geometry_solves=0,
        original_source_byte_hash_available=False,
        source_bytes_origin='Verified trajectory atom records/CRYST1; recovered representation, not a claim of original legacy bytes.',
        solvent_group=solvent)
    del container
    print('Reloading canonical scientific geometry with Network.build forbidden...', flush=True)
    command = Command()
    try:
        with patch.object(Network, 'build', side_effect=AssertionError('Validation must not solve geometry.')):
            net = command._load_canonical_checkpoint(checkpoint)
            groups = [g for g in net.sys.groups if getattr(g, 'solvent_layers', None)]
            assert len(groups) == 1
            assert len(groups[0].solvent_layers['atoms']) == 1
            assert complete_layer_count(groups[0]) == 1
            assert set(groups[0].solvent_layers['atoms'][0]) == set(membership)
        evidence['scientific_checkpoint_reload_passed'] = True
    finally:
        command._cleanup_checkpoint_source()
    evidence['original_archive_unchanged_during_validation'] = digest(original) == original_digest
    assert evidence['original_archive_unchanged_during_validation']
    (bundle / 'recovery_evidence.json').write_text(json.dumps(evidence, indent=2), encoding='utf-8')
    print(f'Validated {len(membership)} one-layer atoms and exact source serials: {bundle}', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', required=True, type=Path)
    parser.add_argument('--original', required=True, type=Path)
    parser.add_argument('--trajectory', required=True, type=Path)
    parser.add_argument('--frame', required=True, type=int)
    args = parser.parse_args()
    validate(args.bundle, args.original, args.trajectory, args.frame)
