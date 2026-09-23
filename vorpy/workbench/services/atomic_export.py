"""Write current-frame atom selections in PDB or XYZ format."""
from pathlib import Path

from vorpy.src.output.pdb import make_pdb_line


def write_atomic_structure(result, indices, destination, name, file_type):
    if file_type not in {'pdb', 'xyz'}:
        raise ValueError('Atom structure format must be PDB or XYZ')
    selected = set(int(index) for index in indices)
    atoms = [atom for atom in result.atoms if atom.index in selected]
    if not atoms:
        return
    path = Path(destination) / f'{name}.{file_type}'
    if file_type == 'xyz':
        lines = [f'{len(atoms)}\n', f'{result.name} frame {result.frame_index}; coordinates in angstroms\n']
        lines.extend(f'{atom.element} {atom.position[0]:.8f} {atom.position[1]:.8f} {atom.position[2]:.8f}\n' for atom in atoms)
    else:
        records = []
        if result.source and result.source.suffix.lower() == '.pdb' and result.source.is_file():
            from vorpy.workbench.services.trajectory import index_pdb_frames, read_frame_lines
            ranges = result.frame_ranges or index_pdb_frames(result.source)
            records = [line for line in read_frame_lines(result.source, ranges[result.frame_index - 1])
                       if line[:6].strip() in {'ATOM', 'HETATM'}]
        lines = [f'HEADER  VorPy {result.name} frame {result.frame_index} {name}\n']
        for atom in atoms:
            if atom.index < len(records):
                lines.append(records[atom.index].rstrip('\n') + '\n')
            else:
                lines.append(make_pdb_line(ser_num=atom.serial, name=atom.name,
                             res_name=atom.residue_name, chain=atom.chain or ' ',
                             res_seq=atom.residue_sequence, elem=atom.element,
                             x=atom.position[0], y=atom.position[1], z=atom.position[2]))
        lines.append('END\n')
    path.write_text(''.join(lines))


def export_atomic_selections(result, options, destination, file_type):
    from vorpy.src.group.export import _surface_neighbor_system_indices
    group = result.export_group
    selections = {}
    if options.get('atoms'):
        selections['group_atoms'] = group.ball_ndxs
    if options.get('surr_atoms'):
        selections['surr_atoms'] = _surface_neighbor_system_indices(group)
    if options.get('surr_resids'):
        group.get_layers(max_layers=1)
        selections['surr_resids'] = group.layer_atoms[1] if len(group.layer_atoms or []) > 1 else []
    if options.get('full_molecule'):
        selections['full_molecule'] = [atom.index for atom in result.atoms]
    for name, indices in selections.items():
        write_atomic_structure(result, indices, destination, name, file_type)
