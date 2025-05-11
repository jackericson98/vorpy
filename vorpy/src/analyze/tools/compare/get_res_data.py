import csv

from Data.Analyze.tools.compare.read_logs import read_logs
from Data.Analyze.tools.compare.pdb_names import proteins, nucleics, ions, other, sols


def residue_data(sys, logs, get_all=False, get_vol=False, get_sa=False, get_curv=False, read_file=None, output_file=None):
    """
    Function that takes in a system and logs and creates a list of sorted residues that can be analyzed for the values
    """
    # Check to see that type of logs it is
    if type(logs) is dict:
        pass
    elif type(logs) is str and logs[-3:] == 'csv':
        logs = read_logs(logs, return_dict=False)
        # print(logs)
    else:
        print("Log File Error: Logs must be in the form of a dictionary from \'read_logs()\' or a \'.csv\' log file "
              "address")
        return
    # print('get_res_data 19: {}'.format(logs))
    # Define the different residue types in their respective dictionary
    protein_dict = {_: {} for _ in proteins}
    nucleic_dict = {_: {} for _ in nucleics}
    ion_dict = {_: {} for _ in ions}
    other_dict = {_: {} for _ in other}
    if read_file is not None:
        # The lines will go: sys_name, res_name, res_seq, vol, sa, max_surf, max_curv, avg_curv
        with open(read_file, 'r') as my_res_file:
            res_reader = csv.reader(my_res_file)
            for line in res_reader:
                if len(line) == 0:
                    continue
                if sys.name != line[0]:
                    continue
                res_info = {'vol': float(line[3]), 'sa': float(line[4]), 'max_surf': line[5],
                            'max_curv': float(line[6]), 'avg_curv': float(line[7])}
                if line[1] in protein_dict:
                    protein_dict[line[1]][line[2]] = res_info
                elif line[1] in nucleic_dict:
                    nucleic_dict[line[1]][line[2]] = res_info
                elif line[1] in ion_dict:
                    ion_dict[line[1]][line[2]] = res_info
                else:
                    if line[1] not in other_dict:
                        other_dict[line[1]] = {}
                    other_dict[line[1]][line[2]] = res_info
        return {'aminos': protein_dict, 'nucs': nucleic_dict, 'ions': ion_dict, 'other': other_dict}
    # Go through the residues in the system and analyze what's inside.
    for res in sys.residues:
        # Create the list of atoms and a surface dictionary lists separated into exterior and interior
        res_atoms, res_surfs = [], {'in': [], 'out': []}

        for i, atom_info_line in logs['atoms'].iterrows():
            # print('get_res_data: 33 {}'.format(atom_info_line))
            # Get the residue atoms information
            if atom_info_line['num'] in res.atoms:
                res_atoms.append(atom_info_line)

        for i, surf_info_line in logs['surfs'].iterrows():
            # Check of one of the surfaces atoms is in the residue
            if surf_info_line['atoms'][0] in res.atoms:
                # print('get_res_data 41: {}'.format(surf_info_line))
                # Check for the other surface
                if surf_info_line['atoms'][1] in res.atoms:
                    res_surfs['in'].append(surf_info_line)
                else:
                    res_surfs['out'].append(surf_info_line)
        # Set the variables to None
        vol, sa, max_curv, avg_curv, max_surf = None, None, None, None, None

        # Calculate the volume
        if get_all or get_vol:
            vol = sum([_['volume'] for _ in res_atoms])
            # print('get_res_data 53: {}'.format(vol))
        # Get the SA
        if get_all or get_sa:
            sa = sum([_['sa'] for _ in res_surfs['out']])
        # Get the curvatures
        if get_all or get_curv:
            if len(res_atoms) == 0:
                max_curv, max_surf, avg_curv = 0, 0, 0
            else:
                max_curv = max([_['max curv'] for _ in res_atoms])
                avg_curv = sum([_['curvature'] for _ in res_surfs['out'] + res_surfs['in']])/len(res_surfs)
                # Get the maximum curvature between residue and other atom
                try:
                    max_surf = max(res_surfs['out'], key=lambda x: x['curvature'])
                except ValueError:
                    max_surf = 0

        # Add the res info for analysis
        res_info = {'vol': vol, 'sa': sa, 'max_surf': max_surf, 'max_curv': max_curv, 'avg_curv': avg_curv}
        # print('get_res_data 72: {}'.format(res_info))
        # Add the information into the appropriate residue dictionary
        if res.name.lower() in protein_dict:
            protein_dict[res.name.lower()][res.seq] = res_info
        elif res.name.lower() in nucleic_dict:
            nucleic_dict[res.name.lower()][res.seq] = res_info
        elif res.name.lower() in ion_dict:
            ion_dict[res.name.lower()][res.seq] = res_info
        elif res.name.lower() in other_dict:
            other_dict[res.name.lower()][res.seq] = res_info
        else:
            other_dict[res.name.lower()] = {res.seq: res_info}
        # Write into the file
        if output_file is not None:
            with open(output_file, 'a') as res_writer:
                my_res_writer = csv.writer(res_writer)
                my_res_writer.writerow([sys.name, res.name.lower(), res.seq, vol, sa, max_surf, max_curv, avg_curv])
    # Return the sorted dictionary with the values
    return {'aminos': protein_dict, 'nucs': nucleic_dict, 'ions': ion_dict, 'other': other_dict}

