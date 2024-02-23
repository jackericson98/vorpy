import csv
import numpy as np

from Data.Analyze.compare_files import compare_files
import matplotlib.pyplot as plt


if __name__ == '__main__':
    prefix = 'C:/Users/jacke/Documents/data1/'
    # my_info = compare_files(pdb_files=[prefix + '1BNA_coarse_ad.pdb'],
    #                         log_files=[prefix + '1BNA_coarse_ad_logs.csv'], avg_distros=True, by_residues=True)
    my_info = compare_files(pdb_files=[prefix + '1BNA.pdb',
                                       prefix + '1BNA.pdb',
                                       prefix + '1BNA_coarse_ad.pdb',
                                       prefix + '1BNA_coarse_ad.pdb',
                                       prefix + '1BNA_coarse_ncap.pdb',
                                       prefix + '1BNA_coarse_ncap.pdb',
                                       prefix + '1BNA_coarse_scbb_ad.pdb',
                                       prefix + '1BNA_coarse_scbb_ad.pdb',
                                       prefix + '1BNA_coarse_scbb_ncap.pdb',
                                       prefix + '1BNA_coarse_scbb_ncap.pdb'],
                            log_files=[prefix + '1BNA_atom_vor_logs.csv',
                                       prefix + '1BNA_atom_pow_logs.csv',
                                       prefix + '1BNA_coarse_ad_vor_logs.csv',
                                       prefix + '1BNA_coarse_ad_pow_logs.csv',
                                       prefix + '1BNA_coarse_ncap_vor_logs.csv',
                                       prefix + '1BNA_coarse_ncap_pow_logs.csv',
                                       prefix + '1BNA_coarse_scbb_ad_vor_logs.csv',
                                       prefix + '1BNA_coarse_scbb_ad_pow_logs.csv',
                                       prefix + '1BNA_coarse_scbb_ncap_vor_logs.csv',
                                       prefix + '1BNA_coarse_scbb_ncap_pow_logs.csv'], avg_distros=True,
                            by_residues=True)

    with open(prefix + '1BNA_residue_data.csv', 'w') as res_file:
        res_fl = csv.writer(res_file)
        res_fl.writerow(['file', 'residue type', 'name', 'residue', 'volume', 'surface area'])
        for file in my_info['residues']:
            for res_type in my_info['residues'][file]:
                for res_name in my_info['residues'][file][res_type]:
                    for res in my_info['residues'][file][res_type][res_name]:
                        res_fl.writerow([file, res_type, res_name, res] + [my_info['residues'][file][res_type][res_name][res][_] for _ in my_info['residues'][file][res_type][res_name][res]])

    vols, sas = {}, {}
    with open(prefix + '1BNA_residue_data.csv', 'r') as res_file:
        res_reader = csv.reader(res_file)
        for i, line in enumerate(res_reader):
            if i == 0:
                continue
            if len(line) == 0 or line[1] == 'other':
                continue
            # File
            if line[0] not in vols:
                vols[line[0]] = {}
                sas[line[0]] = {}
            # Residue Type
            if line[1] not in vols[line[0]]:
                vols[line[0]][line[1]] = {}
                sas[line[0]][line[1]] = {}
            # Residue Class
            if line[2] not in vols[line[0]][line[1]]:
                vols[line[0]][line[1]][line[2]] = []
                sas[line[0]][line[1]][line[2]] = []
            vols[line[0]][line[1]][line[2]].append(line[4])
            sas[line[0]][line[1]][line[2]].append(line[5])

    vol_res_data = {}
    sa_res_data = {}
    res_names = []
    for file in vols:
        vol_res_data[file] = {}
        sa_res_data[file] = {}
        for res_type in vols[file]:
            for res_name in vols[file][res_type]:
                res_names.append(res_name)
                vol_by_type = [float(_) for _ in vols[file][res_type][res_name]]
                vol_res_data[file][res_name] = {'avg': np.mean(vol_by_type), 'max': max(vol_by_type),
                                                'min': min(vol_by_type), 'sd': np.std(vol_by_type), 'data': vol_by_type}
                sa_by_type = [float(_) for _ in sas[file][res_type][res_name]]
                sa_res_data[file][res_name] = {'avg': np.mean(sa_by_type), 'max': max(sa_by_type),
                                               'min': min(sa_by_type), 'sd': np.std(sa_by_type), 'data': sa_by_type}
    titles = ['Atomistic', "Average Distance", 'Encapsulate', 'Sidechain/Backbone - Average Distance', 'Sidechain/Backbone - Encapsulate']

    for i, file in enumerate(sas):
        labels = [_ for _ in sas[file]['nucs']]
        my_data = [[float(__) for __ in sas[file]['nucs'][_]] for _ in labels]
        # Create box and whisker plot
        plt.figure(figsize=(8, 6))
        plt.boxplot(my_data, labels=labels, showmeans=True, vert=True, patch_artist=True, widths=0.6)

        # Add labels and title
        plt.xlabel('Residues')
        plt.ylabel('Surface Area \u212B\u00B2')
        if i % 2 == 0:
            pow_vor_desi = 'AW'
        else:
            pow_vor_desi = 'Pow'
        my_title = '1BNA ' + titles[i // 2] + ' - ' + pow_vor_desi
        plt.title(my_title)

        # Show the plot
        plt.tight_layout()
    plt.show()

