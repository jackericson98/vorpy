"""Join the preserved PEG audit to explicitly verified one-layer recoveries."""
import csv
from collections import Counter
import hashlib
import json
from pathlib import Path


ROOT = next(parent for parent in Path(__file__).resolve().parents
            if (parent / 'pyproject.toml').is_file())
DATA = Path(r'C:\Users\Ryzen-3600\Versiti Dropbox\John Ericson\Versiti lab\Jack\Ericson Voronoi DNA\P2\NewModels\PEG')
OUT = ROOT / 'output/peg_one_layer_investigation_20261112'
REPLAYS = {
    ('aw', 'PEG3000', 'frame_0005'): 'recovered_PEG3000_aw_frame0005_v2',
    ('aw', 'PEG6000', 'frame_0002'): 'representative_PEG6000_aw_frame0002',
    ('aw', 'PEG2000', 'frame_0005'): 'representative_PEG2000_aw_frame0005',
}
ACTIVE = {('aw', 'PEG300', 'frame_0003'), ('prm', 'PEG1000', 'frame_0001')}


def write_rows(path, rows):
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    with (DATA / 'Reruns/completeness_audit_2026-11-12/runs.csv').open(encoding='utf-8-sig') as stream:
        original = list(csv.DictReader(stream))
    rows = []
    for entry in original:
        key = entry['scheme'], entry['molecule'], entry['frame']
        bundle = OUT / REPLAYS[key] if key in REPLAYS else None
        screened = entry['export_membership_checks_pass'] == 'True'
        if key in ACTIVE:
            status = 'EXCLUDED_ACTIVE_RERUN_NOT_INSPECTED'
        elif bundle:
            evidence = json.loads((bundle / 'recovery_evidence.json').read_text())
            assert evidence['geometry_solves'] == 0 and evidence['geometry_cell_fingerprints_unchanged']
            screened = True
            status = 'VERIFIED_ONE_LAYER_REEXPORT'
        elif entry['layer1_complete'] != 'True':
            status = 'HOLD_INCOMPLETE_FIRST_LAYER'
        else:
            assert screened, f'Unclassified original hold: {key}'
            status = 'SCREENED_ATOM_DATA_SELECT_LAYER1'
        usable = screened and key not in ACTIVE
        native_scope = 'ONE_LAYER' if bundle else ('POOLED_TWO_LAYERS' if entry['included_layers'] == '2' else 'ONE_LAYER' if entry['included_layers'] == '1' else 'EMPTY_OR_MISSING')
        if not usable:
            curvature = 'UNAVAILABLE_FOR_COMPLETE_LAYER1_ANALYSIS'
        elif bundle:
            curvature = 'REEXPORTED_COMPONENT_AWARE; residual_is_not_a_convergence_certificate'
        elif entry['scheme'] == 'aw':
            curvature = 'REEXPORT_LAYER1_WITH_COMPONENT_AWARE_CODE_BEFORE_GROUP_CURVATURE_USE'
        elif native_scope == 'POOLED_TWO_LAYERS':
            curvature = 'REEXPORT_LAYER1_BEFORE_GROUP_METRIC_USE'
        else:
            curvature = 'ORIGINAL_SCREEN_ONLY; no_independent_mesh_or_convergence_certification'
        run = DATA / 'Reruns/layers_4' / entry['scheme'] / entry['molecule'] / entry['frame'] / entry['molecule']
        rows.append(dict(scheme=key[0], molecule=key[1], frame=key[2], status=status,
            usable_for_screened_layer1_atom_analysis=usable,
            layer1_atoms=entry['layer1_atoms'], layer1_complete_cells=entry['layer1_complete_cells'],
            native_group_scope=native_scope, group_curvature_status=curvature,
            analysis_directory=str(bundle or run), original_directory=str(run),
            atom_selection='Match layer_01 PDB identity/coordinates to native log Index; do not assume serial minus one',
            limitations='Finite nonperiodic generator model; sum of atom SA is not group union SA; screening is not scientific certification',
            original_issues=entry['issues']))
    rows.sort(key=lambda row: (row['molecule'], row['frame'], row['scheme']))
    counts = Counter(row['status'] for row in rows)
    assert len(rows) == 231 and sum(r['usable_for_screened_layer1_atom_analysis'] for r in rows) == 223
    assert counts['HOLD_INCOMPLETE_FIRST_LAYER'] == 6 and counts['EXCLUDED_ACTIVE_RERUN_NOT_INSPECTED'] == 2
    write_rows(OUT / 'one_layer_run_status.csv', rows)
    common = []
    for molecule, frame in sorted({(r['molecule'], r['frame']) for r in rows}):
        selected = [r for r in rows if (r['molecule'], r['frame']) == (molecule, frame)]
        if len(selected) == 3 and all(r['usable_for_screened_layer1_atom_analysis'] for r in selected):
            common.append(dict(molecule=molecule, frame=frame, schemes='aw;pow;prm',
                scope='Screened layer1 atom data only; match atoms explicitly because shell membership is scheme-dependent'))
    assert len(common) == 69
    write_rows(OUT / 'common_three_scheme_frames.csv', common)

    paths = [(DATA / 'Frames/PEG3000.pdb', 'Original 11-frame trajectory; frame 5 independently matched'),
             (DATA / 'Checkpoints/aw/PEG3000/frame_0005/PEG3000.vpy', 'Original legacy solute geometry; source bytes unavailable')]
    for name in ('PEG3000', 'PEG3000_1', 'PEG3000_2'):
        for path in sorted((ROOT / 'output' / name).rglob('*')):
            if path.is_file():
                paths.append((path, 'Preserved alternate; only PEG3000_2.vpy was used for verified solvent geometry'))
    old = DATA / 'Reruns/layers_4/aw/PEG3000/frame_0005/PEG3000'
    for path in sorted(old.rglob('*')):
        if path.is_file():
            paths.append((path, 'Original conflicted exports; not repaired in place'))
    new = OUT / REPLAYS[('aw', 'PEG3000', 'frame_0005')]
    for path in sorted(new.rglob('*')):
        if path.is_file():
            paths.append((path, 'Verified one-layer recovery; see recovery_evidence.json'))
    inventory = []
    for path, note in paths:
        digest = ''
        if path.suffix in {'.vpy', '.pdb', '.json'}:
            hasher = hashlib.sha256()
            with path.open('rb') as stream:
                for block in iter(lambda: stream.read(1024 * 1024), b''):
                    hasher.update(block)
            digest = hasher.hexdigest()
        inventory.append(dict(path=str(path), bytes=path.stat().st_size, sha256=digest, note=note))
    write_rows(OUT / 'inventory_PEG3000_frame5.csv', inventory)
    print(json.dumps(dict(statuses=counts, screened_layer1_atom_runs=223, common_frames=len(common), inventory_files=len(inventory)), indent=2))


if __name__ == '__main__':
    main()
