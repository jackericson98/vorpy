"""Verified, solve-forbidden one-layer replay into a new directory."""
from __future__ import annotations

import argparse
from copy import copy
import hashlib
import json
import os
from pathlib import Path
import sys
from unittest.mock import patch
import zipfile

import numpy as np

REPOSITORY_ROOT = next(parent for parent in Path(__file__).resolve().parents
                       if (parent / 'pyproject.toml').is_file())
sys.path.insert(0, str(REPOSITORY_ROOT))
if 'figure_reorganization' in Path(__file__).parts:
    from scripts.figure_reorganization.F3.investigate_peg_first_layer import atom_identity, source_frame
else:
    from scripts.investigate_peg_first_layer import atom_identity, source_frame
from vorpy.src.boundary import SolventMode, SOLVENT_RESIDUES
from vorpy.src.command.vpy_cmnd import Command
from vorpy.src.command.calculation_plan import build_export_plan
from vorpy.src.command.command_export import argv_export
from vorpy.src.group.solvent_layers import build_solvent_layers, complete_layer_count
from vorpy.src.io import load_network
from vorpy.src.network import Network


def digest(path):
    hasher = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            hasher.update(block)
    return hasher.hexdigest()


def cell_digest(network):
    hasher = hashlib.sha256()
    for column in ['loc', 'rad', 'num', 'system_num', 'complete', 'vol', 'sa', 'int_gauss_curv']:
        if column in network.balls:
            hasher.update(np.asarray(network.balls[column].tolist(), dtype=float).tobytes())
    return hasher.hexdigest()


def boundary_summary(group):
    fields = ['int_gauss_curv', 'int_gauss_curv_face', 'int_gauss_curv_edge', 'int_gauss_curv_vertex',
              'euler_characteristic', 'boundary_component_count', 'boundary_face_count',
              'boundary_edge_count', 'boundary_vertex_count', 'boundary_is_orientable',
              'boundary_is_complete', 'gauss_bonnet_expected', 'gauss_bonnet_error',
              'gauss_bonnet_relative_error', 'gaussian_curvature_status']
    return dict(name=group.name, atoms=len(group.ball_ndxs),
                face_euler_sum=sum(getattr(group, 'boundary_face_euler_characteristics', {}).values()),
                **{key: getattr(group, key, None) for key in fields})


def recover(archive, trajectory, frame, output, compare_original=False):
    original_directory = Path.cwd()
    archive, trajectory = archive.resolve(), trajectory.resolve()
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    command = Command()
    evidence = dict(archive=str(archive.resolve()), archive_sha256=digest(archive),
                    trajectory=str(trajectory.resolve()), trajectory_sha256=digest(trajectory), frame=frame)
    with zipfile.ZipFile(archive) as source_zip:
        canonical = 'checkpoint/manifest.json' in source_zip.namelist()
    print(f'Loading saved geometry: {archive}', flush=True)
    try:
        network = command._load_canonical_checkpoint(archive) if canonical else load_network(archive)
        system = network.sys
        lines, cryst = source_frame(trajectory, frame)
        if len(lines) != len(system.balls):
            raise ValueError('Source frame and archived System have different atom counts.')
        for index, (line, row) in enumerate(zip(lines, system.balls.itertuples())):
            identity = atom_identity(line)
            if (not np.array_equal(row.loc, identity[4:]) or row.name != identity[0]
                    or row.res.name.strip().upper() != identity[1].upper() or row.res_seq != identity[3]):
                raise ValueError(f'Unverified frame/atom identity at source index {index}.')
            raw_chain = system.balls.iloc[index].get('chain_name')
            if raw_chain is None or str(raw_chain).strip() != identity[2].strip():
                raise ValueError(f'Unverified chain identity at source index {index}.')
        print(f'Verified all {len(lines)} source atom identities and coordinates.', flush=True)
        evidence['verified_atoms'] = len(lines)
        networks = []
        for group in system.groups:
            candidate = group.net or getattr(group, 'source_network', None)
            if candidate is not None and all(candidate is not other for other in networks):
                networks.append(candidate)
        snapshots = [cell_digest(net) for net in networks]
        original_settings = [dict(net.settings) for net in networks]
        for candidate in networks:
            if candidate.settings['net_type'] != 'aw':
                raise ValueError('This isolated recovery is restricted to the audited AW cases.')
            numbers = candidate.balls.get('system_num', candidate.balls['num']).astype(int).to_numpy()
            if set(numbers) != set(range(len(lines))) or len(numbers) != len(lines):
                raise ValueError('The archived generator universe is not exactly the verified physical source atoms.')
            for position, system_id in enumerate(numbers):
                if (not np.array_equal(candidate.balls.iloc[position]['loc'], system.balls.iloc[system_id]['loc'])
                        or candidate.balls.iloc[position]['rad'] != system.balls.iloc[system_id]['rad']):
                    raise ValueError('Network generators disagree with archived System geometry/radii.')
            if getattr(candidate, 'solvent_mode', None) is None:
                candidate.generator_system_indices = tuple(map(int, numbers))
                candidate.solvent_local_indices = tuple(position for position, system_id in enumerate(numbers)
                    if system.balls.iloc[system_id]['res'].name.strip().upper() in SOLVENT_RESIDUES)
                candidate.solvent_mode = SolventMode.SOLVENT_CONTEXT
                candidate.settings['solvent_mode_resolution'] = 'verified_archive_generator_universe'
                candidate.cell_target = 'MOLECULAR'
                candidate.refresh_solvent_provenance()
        evidence['original_network_settings'] = original_settings
        evidence['recovered_network_settings'] = [dict(net.settings) for net in networks]
        evidence['generator_universe'] = 'Every physical source atom, same coordinates and archived radii; no images or artificial generators.'
        source_dir = output / 'verified_source'
        source_dir.mkdir()
        source = source_dir / trajectory.name
        if canonical:
            source.write_bytes(Path(system.files['base_file']).read_bytes())
            evidence['source_bytes_origin'] = 'Existing canonical checkpoint source, independently matched to trajectory frame.'
        else:
            source.write_text((cryst + '\n' if cryst else '') + ''.join(lines) + 'END\n')
            evidence['source_bytes_origin'] = 'Recovered representation from original trajectory atom records/CRYST1; original legacy source-byte hash unavailable.'
        evidence['recovered_source_sha256'] = digest(source)
        system.files['base_file'] = str(source)
        system.files['dir'] = str(output)
        system.solvent_mode = SolventMode.SOLVENT_CONTEXT
        for group in system.groups:
            group.dir = None
        with patch.object(Network, 'build', side_effect=AssertionError('No scientific geometry rebuild is authorized by this replay.')):
            if compare_original:
                evidence['original_logged_cohort_component_aware_reduction'] = []
                for original_group in system.groups:
                    if original_group.net is None:
                        continue
                    group = copy(original_group)
                    if getattr(group, 'solvent_layers', None):
                        included = complete_layer_count(group)
                        group.ball_ndxs = [atom for layer in group.solvent_layers['atoms'][:included] for atom in layer]
                    group._export_info_ready = False
                    group.layer_surfs = None
                    group._export_topology_indices = None
                    print(f'Reducing original logged cohort: {group.name} ({len(group.ball_ndxs)} atoms)...', flush=True)
                    group.get_info()
                    evidence['original_logged_cohort_component_aware_reduction'].append(boundary_summary(group))
            print('Selecting one solvent layer from saved geometry...', flush=True)
            solvent = build_solvent_layers(system, 1)
            if complete_layer_count(solvent) != 1:
                raise ValueError('The saved first layer is incomplete; no analysis bundle was exported.')
            evidence['layer1_atoms'] = len(solvent.solvent_layers['atoms'][0])
            for group in system.groups:
                if group.net is not None:
                    group._export_info_ready = False
                    group.get_info()
            evidence['one_layer_boundary_reduction'] = [boundary_summary(g) for g in system.groups if g.net is not None]
            evidence['geometry_cell_fingerprints_unchanged'] = snapshots == [cell_digest(net) for net in networks]
            if not evidence['geometry_cell_fingerprints_unchanged']:
                raise ValueError('Saved cell geometry changed during replay.')
            system._export_plan = build_export_plan([['large']], scope='group')
            print('Exporting verified one-layer bundle...', flush=True)
            argv_export(system, [['large']])
        evidence['geometry_solves'] = 0
        evidence['input_archive_unchanged'] = digest(archive) == evidence['archive_sha256']
        if not evidence['input_archive_unchanged']:
            raise ValueError('Input archive bytes changed unexpectedly.')
        (output / 'recovery_evidence.json').write_text(json.dumps(evidence, indent=2, default=str), encoding='utf-8')
        print(f'Recovery complete: {output}', flush=True)
    finally:
        command._cleanup_checkpoint_source()
        os.chdir(original_directory)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', required=True, type=Path)
    parser.add_argument('--trajectory', required=True, type=Path)
    parser.add_argument('--frame', required=True, type=int)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--compare-original', action='store_true')
    args = parser.parse_args()
    recover(args.archive, args.trajectory, args.frame, args.output, args.compare_original)
