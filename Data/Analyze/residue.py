

def residue_data(sys, logs, sa=False, vol=False, curv=False, dists=False):
    # We need to sort the atoms into their respective residues
    amino_acids = {'ALA': {}, 'ARB': {}, 'ASN': {}, 'ASP': {}, 'CYS': {}, 'GLN': {}, 'GLU': {}, 'HIS': {}, 'ILE': {},
                   'LEU': {}, 'LYS': {}, 'MET': {}, 'PHE': {}, 'PRO': {}, 'SER': {}, 'THR': {}, 'TRP': {}, 'TYR': {},
                   'VAL': {}, 'GLY': {}, 'ARG': {}}

    nucleic_acids = {'DT': {}, 'DA': {}, 'DG': {}, 'DC': {}, 'DU': {}, 'U': {}, 'G': {}, 'A': {}, 'T': {}, 'C': {}}

    ions = {}

    other = {}
    # We want to collect data from this: Residue volume, residue surface area, residue average curvature,
    # residue maximum curvature, inter-residue sa

    for res in sys.residues:
        # Get the logs atoms for the residue
        res_atoms, res_surfs = [], {'in': [], 'out': []}
        for i, atom_info_line in logs['atoms'].iterrows():
            # Get the residue atoms information
            if atom_info_line['num'] in res.atoms:
                res_atoms.append(atom_info_line)
        # Get the residue surfaces
        for i, surf_info_line in logs['surfs'].iterrows():
            # Check of one of the surfaces atoms is in the residue
            if surf_info_line['atoms'][0] in res.atoms:
                # Check for the other surface
                if surf_info_line['atoms'][1] in res.atoms:
                    res_surfs['in'].append(surf_info_line)
                else:
                    res_surfs['out'].append(surf_info_line)
        # Calculate the volume and surface area
        vol = sum([_['volume'] for _ in res_atoms])
        sa = sum([_['sa'] for _ in res_surfs['out']])
        # Get the curvatures
        # max_curv = max([_['max curv'] for _ in res_atoms])
        # avg_curv = sum([_['curvature'] for _ in res_surfs['out'] + res_surfs['in']])/len(res_surfs)
        # Get the maximum curvature between residue and other atom
        # max_surf = max(res_surfs['out'], key=lambda x: x['curvature'])

        # Add the res info for analysis
        res_info = {'vol': vol, 'sa': sa}
                    # 'max_contact': [_ for _ in max_surf['atoms'] if _ not in res.atoms][0],
                    # 'max curv out': max_surf['curvature']}
        if res.name in amino_acids:
            amino_acids[res.name][res.seq] = res_info
        elif res.name in nucleic_acids:
            nucleic_acids[res.name][res.seq] = res_info
        elif res.name in ions:
            ions[res.name][res.seq] = res_info
        else:
            if res.name in other:
                other[res.name][res.seq] = res_info
            else:
                other[res.name] = {res.seq: res_info}
    return {'aminos': amino_acids, 'nucs': nucleic_acids, 'ions': ions, 'other': other}

