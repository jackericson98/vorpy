from System.system import System
from Data.Analyze.read_logs import read_logs


def compare_files(pdb_files, log_files, show_overview=True, show_group_data=True, compare_in_out=True, compare_residues=True):
    # Create the System
    # Read the files
    systems, logs = [], []
    print('Comparing: ')
    for i, file in enumerate(log_files):
        # Get the scheme

        systems.append(System(file=pdb_files[i]))
        logs.append(read_logs([log_files[i]]))
        print('File {}'.format(i + 1), systems[-1].name)

    # Compare build settings
    if show_overview:
        print("\nBuild Settings\n")
        print(*[_ for _ in logs[0][systems[0].name]['data']])
        for i, sys in enumerate(systems):
            print(*[logs[i][sys.name]['data'][_] for _ in logs[i][sys.name]['data']])

    # Compare full group values for the systems
    if show_group_data:
        print('\nGroup Data\n')
        print('\nTotal Volume, Surface Area\n')
        for i, sys in enumerate(systems):
            print(sys.name, logs[i][sys.name]['group data']['volume'], logs[i][sys.name]['group data']['sa'])

    # Compare inside vs outside:
    if compare_in_out:
        print('\nInside Vs. Outside Data')
        print('\nAverage Curvature Out, Average Curvature In, out vol H, out vol C, out vol ')
        for i, sys in enumerate(systems):
            print(in_out_data())


def in_out_data(sys, logs):
    # First we need to designate the atoms in each group that are in the outside and atoms that are on the inside
    in_atoms, out_atoms = [], []
    for atom in logs['atoms']:
        out = False
        for _ in atom['neighbors']:
            if _ not in sys.groups[0].atoms:
                out = True
        if out:
            out_atoms.append(atom['num'])
            atom['in'] = False
        else:
            in_atoms.append(atom['num'])
            atom['in'] = True
    in_vol, out_vol, in_sa, out_sa, in_vols, out_vols = 0, 0, 0, 0, [], []
    # Get the data for each atom
    for atom in logs['atoms']:
        if atom['in']:
            in_sa += atom['sa']
            in_vols.append(atom['vol'])
        else:
            out_vols.append(atom['vol'])
            out_sa += atom['sa']
    return {}


def residue_data(sys, logs):
    # We need to sort the atoms into their respective residues

    amino_acids = {'ALA': {}, 'ARB': {}, 'ASN': {}, 'ASP': {}, 'CYS': {}, 'GLN': {}, 'GLU': {}, 'HIS': {}, 'ILE': {},
                   'LEU': {}, 'LYS': {}, 'MET': {}, 'PHE': {}, 'PRO': {}, 'SER': {}, 'THR': {}, 'TRP': {}, 'TYR': {},
                   'VAL': {}, 'GLY': {}, 'ARG': {}}

    nucleic_acids = {'DT': {}, 'DA': {}, 'DG': {}, 'DC': {}, 'DU': {}, 'U': {}, 'G': {}, 'A': {}, 'T': {}, 'C': {}}

    ions = {}

    other = {}
    # We want to collect data from this: Residue volume, residue surface area, residue average curvature,
    # residue maximum curvature, inter-residue sa

    # Check if the system is coarse grained or not
    for res in sys.residues:
        # Get the logs atoms for the residue
        res_atoms, res_surfs = [], {'in': [], 'out': []}
        for atom_info_line in logs['atoms']:
            # Get the residue atoms information
            if atom_info_line['num'] in res.atoms:
                res_atoms.append(atom_info_line)
        # Get the residue surfaces
        for surf_info_line in logs['surfs']:
            # Check of one of the surfaces atoms is in the residue
            if surf_info_line['atoms'][0] in res.atoms:
                # Check for the other surface
                if surf_info_line['atoms'][1] in res.atoms:
                    res_surfs['in'].append(surf_info_line)
                else:
                    res_surfs['out'].append(surf_info_line)
        # Calculate the volume and surface area
        vol = sum([_['vol'] for _ in res_atoms])
        sa = sum([_['sa'] for _ in res_surfs['out']])
        # Get the curvatures
        max_curv = max([_['max curvature'] for _ in res_atoms])
        avg_curv = sum([_['curvature'] for _ in res_surfs['out'] + res_surfs['in']])/len(res_surfs)
        # Get the maximum curvature between residue and other atom
        max_surf = max(res_surfs['out'], key=lambda x: x['curvature'])

        # Add the res info for analysis
        res_info = {'vol': vol, 'sa': sa, 'max curv tot': max_curv, 'avg curv': avg_curv,
                    'max_contact': [_ for _ in max_surf['atoms'] if _ not in res.atoms][0],
                    'max curv out': max_surf['curvature']}
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
    return {'amines': amino_acids, 'nucs': nucleic_acids, 'ions': ions, 'other': other}





if __name__ == '__main__':
    prefix = 'C:/Users/i7-8700/Documents/test_files/'
    compare_files([prefix + 'atomistic/181L.pdb', prefix + 'avg_dist/181L_coarse_ad.pdb'],
                  [prefix + 'atomistic/181L_vor/sys/181L.csv',
                   prefix + 'avg_dist/181L_coarse_ad_vor/sys/181L_coarse_ad_logs.csv'])