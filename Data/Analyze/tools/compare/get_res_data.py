from Data.Analyze.tools.compare.read_logs import read_logs
from Data.Analyze.tools.compare.pdb_names import proteins, nucleics, ions, other, sols


def residue_data(sys, logs, get_all=False, get_vol=False, get_sa=False, get_curv=False):
    """
    Function that takes in a system and logs and creates a list of sorted residues that can be analyzed for the values
    """
    # Check to see that type of logs it is
    if type(logs) is dict:
        pass
    elif type(logs) is str and logs[-3:] == 'csv':
        logs = read_logs(logs, new_logs=True)
    else:
        print("Log File Error: Logs must be in the form of a dictionary from \'read_logs()\' or a \'.csv\' log file "
              "address")
        return

    # Define the different residue types in their respective dictionary
    protein_dict = {_: {} for _ in proteins}
    nucleic_dict = {_: {} for _ in nucleics}
    ion_dict = {_: {} for _ in ions}
    other_dict = {_: {} for _ in other}

    # Go through the residues in the system and analyze what's inside.
    for res in sys.residues:
        # Create the list of atoms and a surface dictionary lists separated into exterior and interior
        res_atoms, res_surfs = [], {'in': [], 'out': []}

        for atom_info_line in logs['atoms']:
            # Get the residue atoms information
            if atom_info_line['num'] in res.atoms:
                res_atoms.append(atom_info_line)

        for surf_info_line in logs['surfs']:
            # Check of one of the surfaces atoms is in the residue
            if surf_info_line['atoms'][0] in res.atoms:
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

        # Add the information into the appropriate residue dictionary
        if res.name in protein_dict:
            protein_dict[res.name][res.seq] = res_info
        elif res.name in nucleic_dict:
            nucleic_dict[res.name][res.seq] = res_info
        elif res.name in ion_dict:
            ion_dict[res.name][res.seq] = res_info
        elif res.name in other_dict:
            other_dict[res.name][res.seq] = res_info
        else:
            other_dict[res.name] = {res.seq: res_info}

    # Return the sorted dictionary with the values
    return {'aminos': protein_dict, 'nucs': nucleic_dict, 'ions': ions, 'other': other}

