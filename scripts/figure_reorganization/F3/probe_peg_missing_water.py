"""Isolate one missing water's cells with the original full generator context."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

REPOSITORY_ROOT = next(parent for parent in Path(__file__).resolve().parents
                       if (parent / 'pyproject.toml').is_file())
sys.path.insert(0, str(REPOSITORY_ROOT))
if 'figure_reorganization' in Path(__file__).parts:
    from scripts.figure_reorganization.F3.investigate_peg_first_layer import atom_identity, source_frame, geometry_archive, read_geometry
else:
    from scripts.investigate_peg_first_layer import atom_identity, source_frame, geometry_archive, read_geometry
from vorpy.src.group import Group
from vorpy.src.system import System
from vorpy.src.network.analyze import diagnose_geometry_completeness
from vorpy.src.io import save_network


def probe(base, scheme, frame, output, control_water=False, molecule='PEG4000', indices=None):
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    lines, cryst = source_frame(base / f'Frames/{molecule}.pdb', frame)
    source = output / f'{molecule}.pdb'
    source.write_text(cryst + '\n' + ''.join(lines) + 'END\n')
    system = System(file=str(source), make_dir=False, solvent_mode='SOLVENT_CONTEXT')
    system.files['dir'] = str(output.resolve())
    archive = base / f'Reruns/layers_4/{scheme}/{molecule}/frame_{frame:04d}/{molecule}/{molecule}_checkpoint.vpy'
    with geometry_archive(archive) as saved:
        old, metadata = read_geometry(saved)
    assert np.array_equal(np.asarray(old.balls['loc'].tolist()), np.asarray(system.balls['loc'].tolist()))
    assert np.array_equal(old.balls['rad'], system.balls['rad'])
    settings = dict(net_type=scheme, surf_res=metadata['surface_resolution'],
                    max_vert=metadata['maximum_vertex_radius'], box_size=metadata['boundary_padding'],
                    sys_dir=str(output.resolve()))
    selected = list(indices) if indices is not None else [8796, 8797, 8798]
    reference = selected[0]
    if control_water:
        identities = {atom_identity(line): i for i, line in enumerate(lines)}
        root = archive.parent
        layer_ids = {identities[atom_identity(line)] for line in (root / 'solvent/layer_01.pdb').read_text().splitlines()
                     if line.startswith(('ATOM  ', 'HETATM'))}
        for index in sorted(layer_ids):
            if lines[index][12:16].strip() != 'OW' or index in selected:
                continue
            members = [i for i in layer_ids if lines[i][17:27] == lines[index][17:27]]
            if len(members) != 3 or not old.balls.loc[members, 'complete'].all():
                continue
            if np.linalg.norm(system.balls.iloc[index]['loc'] - system.balls.iloc[reference]['loc']) > 10:
                selected = sorted(members) + selected
                break
        assert len(selected) == 6, 'No independently complete distant control water was found.'
    group = Group(system, name='missing_SOL2723', atoms=selected, settings=settings,
                  solvent_mode='SOLVENT_CONTEXT', cell_target='MOLECULAR')
    system.groups = [group]
    group.build()
    assert len(group.net.balls) == len(system.balls)
    assert np.array_equal(np.asarray(group.net.balls['loc'].tolist()), np.asarray(system.balls['loc'].tolist()))
    rows = group.net.balls.loc[group.net.balls['system_num'].isin(selected)]
    result = dict(scheme=scheme, molecule=molecule, frame=frame, unchanged_physical_generators=len(system.balls),
                  original_settings=settings, control_water=control_water, cells=[dict(source_index=int(row.system_num),
                  complete=bool(row.complete), vertices=len(row.verts), edges=len(row.edges),
                  faces=len(row.surfs), volume=float(row.vol), area=float(row.sa)) for row in rows.itertuples()],
                  search_timing=group.net.vert_timing)
    if not rows['complete'].all():
        result['geometry_diagnostics'] = diagnose_geometry_completeness(group.net)
        result['boundedness_caveat'] = 'The legacy homogeneous recession diagnostic is not an independent additive-weighted boundedness certificate.'
        save_network(group.net, output / 'diagnostic_incomplete_network.vpy', write_timing_manifest=False)
    (output / 'result.json').write_text(json.dumps(result, indent=2, default=str))
    print(json.dumps(result['cells'], indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--scheme', required=True, choices=['aw', 'prm'])
    parser.add_argument('--frame', required=True, type=int)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--control-water', action='store_true')
    parser.add_argument('--molecule', default='PEG4000')
    parser.add_argument('--indices', nargs='+', type=int)
    args = parser.parse_args()
    probe(args.root, args.scheme, args.frame, args.output, args.control_water, args.molecule, args.indices)
