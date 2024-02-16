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
        print('\nAverage Curvature, ')





if __name__ == '__main__':
    prefix = 'C:/Users/i7-8700/Documents/test_files/'
    compare_files([prefix + 'atomistic/181L.pdb', prefix + 'avg_dist/181L_coarse_ad.pdb'],
                  [prefix + 'atomistic/181L_vor/sys/181L.csv',
                   prefix + 'avg_dist/181L_coarse_ad_vor/sys/181L_coarse_ad_logs.csv'])