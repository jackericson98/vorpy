"""Read-only component-aware topology audit of flagged PEG AW checkpoints."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import sys
from types import SimpleNamespace
import zipfile

import numpy as np

REPOSITORY_ROOT = next(parent for parent in Path(__file__).resolve().parents
                       if (parent / 'pyproject.toml').is_file())
sys.path.insert(0, str(REPOSITORY_ROOT))
if 'figure_reorganization' in Path(__file__).parts:
    from scripts.figure_reorganization.F3.investigate_peg_first_layer import table
    from scripts.figure_reorganization.F3.validate_peg_layer1_bundle import atom_log
else:
    from scripts.investigate_peg_first_layer import table
    from scripts.validate_peg_layer1_bundle import atom_log
from vorpy.src.calculations.boundary_complex import boundary_complex
from vorpy.src.io.archive_codec import Reader


def inspect_archive(path, native_root):
    with zipfile.ZipFile(path) as outer:
        manifest = json.loads(outer.read('checkpoint/manifest.json'))
        state = max(manifest['document']['networks'], key=lambda n: sum(n['completeness'].values()))
        payload = next(p for p in manifest['payloads'] if p['payload_id'] == state['payload_ref'])
        content = outer.read(payload['member'])
        if hashlib.sha256(content).hexdigest() != payload['sha256']:
            raise ValueError('Geometry payload hash mismatch.')
    results = []
    with zipfile.ZipFile(io.BytesIO(content)) as archive:
        reader = Reader(archive)
        for record in reader.json('state.json')['objects']:
            if record['kind'] != 'network':
                continue
            fields = record['fields']
            role = 'solvent' if fields['group_name'].startswith('solvent') else 'solute'
            logs = list((native_root / 'groups').glob('*/*_logs.csv')) if role == 'solute' else list((native_root / 'solvent').glob('*_logs.csv'))
            if len(logs) != 1:
                raise ValueError(f'Expected one {role} log, found {len(logs)}.')
            atoms, group = atom_log(logs[0])
            result = dict(role=role, native_atoms=len(atoms), original_chi=group['Euler Characteristic'],
                          measured_gaussian=group['Integrated Gaussian Curvature'],
                          original_expected=group['Gauss-Bonnet Expected'],
                          original_relative_error=group['Gauss-Bonnet Relative Error'],
                          original_orientable=group['Boundary Orientable'])
            if not atoms:
                results.append(dict(result, status='EMPTY_LOGGED_COHORT_NOT_A_NONORIENTABLE_SURFACE'))
                continue
            balls = table(reader, fields['balls']['table'], {'num', 'system_num'})
            mapping = dict(zip(balls['num'], balls.get('system_num', balls['num'])))
            body = {int(atom['Index']) for atom in atoms}
            if not body.issubset(mapping.values()):
                raise ValueError('Native atom selection is outside archived physical generators.')
            net = SimpleNamespace(
                edges=table(reader, fields['edges']['table'], {'verts'}),
                surfs=table(reader, fields['surfs']['table'], {'balls', 'edges', 'tris'}))
            surfaces = [int(s.Index) for s in net.surfs.itertuples()
                        if sum(mapping[int(ball)] in body for ball in s.balls) == 1]
            edges = {int(e) for s in surfaces for e in net.surfs.loc[s, 'edges']}
            vertices = {int(v) for e in edges for v in net.edges.loc[e, 'verts']}
            topology = boundary_complex(net, surfaces)
            measured = float(result['measured_gaussian'])
            unresolved = topology['unresolved_faces']
            chi = len(vertices) - len(edges) + sum(topology['face_euler'].values()) if not unresolved else None
            expected = 2 * np.pi * chi if chi is not None else None
            relative = abs(measured - expected) / abs(expected) if expected else None
            result.update(status='UNRESOLVED_FACE_TOPOLOGY' if unresolved else 'MEASURED_CURVATURE_UNCHANGED',
                surface_records=len(surfaces), disk_patches=topology['patch_count'],
                vertices=len(vertices), edges=len(edges), corrected_chi=chi,
                corrected_expected=expected, corrected_absolute_error=measured - expected if expected is not None else None,
                corrected_relative_error=relative, corrected_orientable=topology['orientable'],
                corrected_components=topology['components'] if not unresolved else None, unresolved_faces=unresolved,
                above_1pct_after_correction=(relative > .01) if relative is not None and not unresolved else None)
            results.append(result)
    return results


def audit(base, output, limit=None):
    output.mkdir(parents=True, exist_ok=False)
    with (base / 'Reruns/completeness_audit_2026-11-12/runs.csv').open(newline='') as stream:
        rows = [r for r in csv.DictReader(stream) if r['scheme'] == 'aw' and r['curvature_review'] == 'True']
    if limit:
        rows = rows[:limit]
    all_results = []
    for index, row in enumerate(rows, 1):
        native = base / f"Reruns/layers_4/aw/{row['molecule']}/{row['frame']}/{row['molecule']}"
        path = native / f"{row['molecule']}_checkpoint.vpy"
        label = f"aw_{row['molecule']}_{row['frame']}"
        print(f'{index}/{len(rows)}: {label} (read-only)', flush=True)
        if not path.exists():
            all_results.append(dict(case=label, status='ORIGINAL_CHECKPOINT_UNAVAILABLE'))
            continue
        results = inspect_archive(path, native)
        for result in results:
            result['case'] = label
            result['native_included_layers'] = row['included_layers'] if result['role'] == 'solvent' else 0
            all_results.append(result)
        (output / f'{label}.json').write_text(json.dumps(results, indent=2), encoding='utf-8')
    (output / 'all_results.json').write_text(json.dumps(all_results, indent=2), encoding='utf-8')
    columns = sorted({k for r in all_results for k in r if k != 'unresolved_faces'})
    with (output / 'boundary_corrections.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(all_results)
    print(f'Finished {len(rows)} saved-geometry inspections; no solves or re-exports.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--limit', type=int)
    args = parser.parse_args()
    audit(args.root, args.output, args.limit)
