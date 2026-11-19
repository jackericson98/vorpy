"""Prepare a recovery ledger from saved audits; never launches calculations."""
import argparse
from collections import Counter
import csv
import importlib.util
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--examples',type=Path,required=True)
    parser.add_argument('--selection',type=Path,required=True)
    parser.add_argument('--runs',type=Path,default=Path('output/Figure3B_revalidation'))
    args = parser.parse_args()
    repository = next(parent for parent in Path(__file__).resolve().parents
                      if (parent / 'pyproject.toml').is_file())
    spec = importlib.util.spec_from_file_location('figure3b',args.examples/'generate_figure3b_pymol.py')
    figure = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(figure)
    baseline = list(csv.DictReader((args.examples/'figure3b_validation.csv').read_text().splitlines()))
    selection = list(csv.DictReader(args.selection.read_text().splitlines()))
    recovered = json.loads((args.runs/'real_pair_validation.json').read_text())
    diagnostics = json.loads((args.runs/'original_incomplete_diagnostics.json').read_text())
    diagnostic_lookup = {(r['example'],r['scheme']):r['diagnostic'] for r in diagnostics}
    selected = {}
    for row in selection:
        molecule = row['System'].split('_')[1]
        name = molecule+'_'+str(int(row['Index']))
        source = next(p for p in (repository/'vorpy/data').glob('*.pdb') if p.stem.casefold()==molecule.casefold())
        atom = figure.read_pdb(source)[int(row['Index'])]
        assert (atom['name'],atom['resn'],atom['resi'],atom['chain']) == (row['Atom'],row['Residue'],row['Residue Sequence'],row['Chain']), name
        selected[name] = dict(example=name,molecule=molecule,central_atom_index=int(row['Index']),
                              source_selection_system=row['System'],source_pdb=str(source),
                              selection_direction=row['Selection Direction'],selection_rank=row['Selection Rank'],
                              identity_validation='source position and chemical identity matched',
                              existing_directory=(args.examples/name).is_dir())
    figure.write_csv(args.examples/'figure3b_priority_atoms.csv',list(selected.values()))
    known = {r['example'] for r in baseline}
    rows = list(baseline)
    for name,row in selected.items():
        if name not in known:
            for scheme in figure.SCHEMES:
                rows.append(dict(example=name,molecule=row['molecule'],central_atom_index=row['central_atom_index'],scheme=scheme,
                                 cell_completeness='unknown',output_availability='missing directory',mapping_validity='unknown'))
    ledger = []
    for row in rows:
        name, scheme = row['example'], row['scheme']
        category = '0 passed existing audit'
        cause = 'Complete saved cell and valid OFF boundary'
        action, recalculate = 'Retain scientific outputs', 'No'
        diag = diagnostic_lookup.get((name,scheme),{})
        if row['output_availability'] != 'available':
            category, cause = '5 missing calculation outputs', 'Only vertex artifacts or no example directory; no complete calculation export'
            action, recalculate = 'Diagnose/run one bounded selected-cell calculation; do not infer completeness from vertex files', 'Yes'
        elif row.get('log_atomic_records') == '0' or row.get('log_boundary_closed') == 'False':
            category, cause = '1 genuinely open cell', 'Missing constructed topology; exact mechanism not established for this atom'
            action, recalculate = 'Inspect search bounds and all-generator vertex validity before another bounded solve', 'Undetermined'
            if diag.get('bounded') is False:
                cause = 'UNBOUNDED in MOLECULE_ONLY according to weighted recession-cone diagnostic'
                action, recalculate = 'Choose the intended physical generator universe; no finite cutoff can close an unbounded cell', 'Not sufficient'
            if (name,scheme) == ('Hammerhead_51505','aw'):
                cause = 'Valid vertex clearance 5.363195 A excluded by max_vert=5 A'
                action, recalculate = 'Use validated mv=6 run; no topology fabrication', 'Yes (completed)'
            elif (name,scheme) == ('Hammerhead_51505','pow'):
                cause = 'Vertex search cutoff; mv=6 recovers topology, plus planar meshing/reduction defects'
                action, recalculate = 'Use validated mv=6 run and corrected planar layers', 'Yes (completed)'
            elif (name,scheme) == ('EDTA_6','aw'):
                cause = 'Valid vertices lie outside search box at box_size=1.25; below the vertex-radius cutoff'
                action, recalculate = 'Use validated bs=2 run, keeping mv=5 and all original generators', 'Yes (completed)'
        elif row.get('log_boundary_complete') == 'False':
            category, cause = '4 incorrect completeness detection', 'AW-only curvature records required for planar boundaries; flat tessellation also has export defects'
            action, recalculate = 'Retain archived vertices/topology; rebuild planar faces and re-export with corrected code, then validate', 'No vertex solve required'
        elif row['cell_completeness'] != 'complete':
            if not int(row.get('open_edges') or 0) and not int(row.get('nonmanifold_edges') or 0):
                category, cause = '6 multiple-loop AW face topology', 'One pairwise AW surface has two boundary cycles; one-face Euler reduction is invalid'
                action, recalculate = 'Retain archive for multi-patch face investigation; do not override orientability/completeness flags', 'Not yet justified'
            else:
                category, cause = '2 invalid exported boundary', 'Atomic topology marked complete but OFF contains open/nonmanifold mesh segments'
                action, recalculate = 'Retain archive; inspect canonical boundary samples and export rounding before deciding to rerun solver', 'Undetermined'
        secondary = '3 missing identity sidecar; recover from saved archive, never a serial offset' if row.get('mapping_validity') == 'unresolved' else ''
        if secondary:
            action += '; reserialize archived topology/system identity mapping'
        status = 'NOT REVALIDATED'
        verified = recovered.get(name,{}).get(scheme)
        if verified:
            status = 'PASS: independent geometry, mapping, CSV/OFF, Figure3B generator'
            action = 'Use new validated output; retain original artifact unchanged'
        molecule = row['molecule']
        pdb = next(p for p in (repository/'vorpy/data').glob('*.pdb') if p.stem.casefold() == molecule.casefold())
        limit, box = (6 if name == 'Hammerhead_51505' else 5), (2 if name == 'EDTA_6' else 1.25)
        destination = (args.runs/'remaining'/name/scheme).resolve()
        command = (f"& '{repository/'.venv/Scripts/python.exe'}' '{repository/'vorpy'}' '{pdb}' -g a {row['central_atom_index']} "
                   f"--solvent-mode MOLECULE_ONLY -s nt {scheme} -s sr 0.05 -s mv {limit} -s bs {box} -e large -e dir '{destination}'")
        role = 'validated parameter replay' if verified else 'bounded fresh-build option' if category.startswith(('0','4')) else 'diagnostic command only; recovery parameters/universe not yet approved by evidence'
        ledger.append(dict(molecule=molecule,central_atom_index=row['central_atom_index'],scheme=scheme,example=name,
                           priority='1 original selection' if name in selected else '2 historical/control',
                           required_for_initial_figure=scheme in ('aw','pow'),
                           existing_status=row['cell_completeness'],failure_category=category,secondary_failure=secondary,
                           root_cause=cause,recalculation_required=recalculate,recommended_action=action,validation_status=status,
                           new_validated_directory=str((args.runs/'validated'/name).resolve()) if verified else '',
                           cli_command=command,command_role=role))
    assert len(ledger) == len(baseline)+3*len(set(selected)-known)
    assert len({(r['example'],r['scheme']) for r in ledger}) == len(ledger)
    ledger.sort(key=lambda r:(r['priority'],r['example'],figure.SCHEMES.index(r['scheme'])))
    figure.write_csv(args.examples/'figure3b_recovery_plan.csv',ledger)
    updated = []
    for original in baseline:
        row = dict(original)
        row['original_cell_completeness'] = original['cell_completeness']
        row['source_directory'] = str((args.examples/original['example']).resolve())
        verified = recovered.get(original['example'],{}).get(original['scheme'])
        if verified:
            row.update(verified['export_audit'])
            row['source_directory'] = str((args.runs/'validated'/original['example']).resolve())
        updated.append(row)
    figure.write_csv(args.runs/'figure3b_validation_recovery.csv',updated)
    counts = Counter(r['failure_category'] for r in ledger if r['example'] in known)
    summary = dict(original_outputs=len(baseline),primary_categories=counts,
                   secondary_missing_mapping=sum(bool(r['secondary_failure']) for r in ledger),
                   original_selected_atoms=len(selected),missing_selected_directories=sorted(set(selected)-known),
                   ledger_rows=len(ledger),validated_pairs=sorted(recovered),
                   decision='GO for bounded finite cases with validated parameters/corrected planar layers; NO-GO for an unattended full batch')
    (args.runs/'recovery_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
