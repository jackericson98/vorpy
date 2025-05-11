from vorpy.src.system.system import System
from Data.Analyze.tools.compare.read_logs import read_logs
from Data.Analyze.tools.compare.get_res_data import residue_data


def compare_files(pdb_files, log_files, build_data=False, totals=False, avg_distros=False, compare_in_out=False,
                  by_residues=False, vol=False, sa=False, curv=False, by_element=False, by_chain=False):
    """
    Takes in log files and pdbs and returns a dictionary of information

    """
    # Instantiate the return dictionary
    info = {'files': []}
    # Create the systems and logs lists
    systems, logs, names = [], [], []

    for i, file in enumerate(log_files):
        # Get the scheme
        my_sys = System(file=pdb_files[i])
        counter = 0
        while True:
            if my_sys.name in names:
                my_sys.name = my_sys.name + '_' + str(counter)
            else:
                break
        names.append(my_sys.name)
        systems.append(my_sys)
        try:
            my_logs = read_logs([log_files[i]])
        except ValueError:
            my_logs = read_logs([log_files[i]], new_logs=True)
        logs.append([my_logs[_] for _ in my_logs][0])
        line = {'num': (i + 1), 'name': systems[-1].name, 'net type': logs[-1]['data']['network_type']}
        info['files'].append(line)

    # Add build settings
    if build_data:
        info['build settings'] = {}
        for i, sys in enumerate(systems):
            info['build settings'][sys.name] = logs[i]['data']

    # Compare full group values for the systems
    if totals:
        info['totals'] = {}
        for i, sys in enumerate(systems):
            avg_curv = 0
            if curv:
                avg_curv = sum(logs[i]['surfs']['curvature']) / len(logs[i]['surfs'])
            line = {'name': sys.name, 'vol': logs[i]['group data']['volume'], 'sa': logs[i]['group data']['sa'],
                    'max curv': max(logs[i]['atoms']['max curv']), 'avg curv': avg_curv}
            info['totals'][sys.name] = line

    # Get distribution data for specific data_type
    if avg_distros:
        # Get the average max curvature, average volume, average outer sa for each residue type
        if by_residues:
            info['residues'] = {}
            for i, sys in enumerate(systems):
                line = residue_data(sys, logs[i], True, True, True)
                info['residues'][sys.name] = line

    # Get the average Curvature

        # Get the average max curvature, average volume, average outer sa for each atom type

    # # Compare inside vs outside:
    #     # if compare_in_out:
    #     #     print('\nInside Vs. Outside Data')
    #     #     csv_lines.append('Inside Vs. Outside Data')
    #     #     print('\nAverage Curvature Out, Average Curvature In, Average Volume Out, Average Volume In, Total Volume out, Total Volume in, ')
    #     #     for i, sys in enumerate(systems):
    #     #         print(in_out_data(sys, logs[i]))
    #     #
    #     # # Compare residues
    #     # if compare_residues:
    #     #     print('\nResidue Data:')
    #     #     print('\n')
    #     #     res_info_by_sys = []
    #     #     for i, sys in enumerate(systems):
    #     #         print(residue_data(sys, logs[i]))
    return info

