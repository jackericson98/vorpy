import csv
import numpy as np

from Data.Analyze.compare.compare_files import compare_files
import matplotlib.pyplot as plt


if __name__ == '__main__':
    prefix = 'C:/Users/jacke/Documents/data/'
    # my_info = compare_files(pdb_files=[prefix + '181L_coarse_ad.pdb'],
    #                         log_files=[prefix + '181L_coarse_ad_logs.csv'], avg_distros=True, by_residues=True)
    # my_info = compare_files(pdb_files=[prefix + '181L.pdb',
    #                                    prefix + '181L.pdb',
    #                                    prefix + '181L_coarse_ad.pdb',
    #                                    prefix + '181L_coarse_ad.pdb',
    #                                    prefix + '181L_coarse_ncap.pdb',
    #                                    prefix + '181L_coarse_ncap.pdb',
    #                                    prefix + '181L_coarse_scbb_ad.pdb',
    #                                    prefix + '181L_coarse_scbb_ad.pdb',
    #                                    prefix + '181L_coarse_scbb_ncap.pdb',
    #                                    prefix + '181L_coarse_scbb_ncap.pdb',
    #                                    prefix + '181L_martini.pdb',
    #                                    prefix + '181L_martini.pdb'],
    #                         log_files=[prefix + '181L_atom_vor_logs.csv',
    #                                    prefix + '181L_atom_pow_logs.csv',
    #                                    prefix + '181L_coarse_ad_logs.csv',
    #                                    prefix + '181L_coarse_ad_pow_logs.csv',
    #                                    prefix + '181L_coarse_ncap_vor_logs.csv',
    #                                    prefix + '181L_coarse_ncap_pow_logs.csv',
    #                                    prefix + '181L_coarse_scbb_ad_vor_logs.csv',
    #                                    prefix + '181L_coarse_scbb_ad_pow_logs.csv',
    #                                    prefix + '181L_coarse_scbb_ncap_vor_logs.csv',
    #                                    prefix + '181L_coarse_scbb_ncap_pow_logs.csv',
    #                                    prefix + '181L_martini_vor_logs.csv',
    #                                    prefix + '181L_martini_pow_logs.csv'], avg_distros=True, by_residues=True)

    # with open(prefix + '181L_residue_data.csv', 'w') as res_file:
    #     res_fl = csv.writer(res_file)
    #     res_fl.writerow(['file', 'residue type', 'name', 'residue', 'volume', 'surface area'])
    #     for file in my_info['residues']:
    #         for res_type in my_info['residues'][file]:
    #             for res_name in my_info['residues'][file][res_type]:
    #                 for res in my_info['residues'][file][res_type][res_name]:
    #                     res_fl.writerow([file, res_type, res_name, res] + [my_info['residues'][file][res_type][res_name][res][_] for _ in my_info['residues'][file][res_type][res_name][res]])

    vols, sas = {}, {}
    with open(prefix + '181L_residue_data.csv', 'r') as res_file:
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


    # Sample data
    labels = ['Atomistic', 'Average Distance', 'Encapsulate', 'Side-Chain/Backbone AD', 'Side-Chain/Backbone Encap', 'CG Martini']
    # Get the percent difference for each residue and average
    data = [0]
    for file in sas:
        if file == '181L':
            continue
        my_perc_diff = []

        for res_name in sas[file]['amines']:
            for i in range(len(sas[file]['amines'][res_name])):
                my_perc_diff.append(abs(float(sas['181L']['amines'][res_name][i]) - float(sas[file]['amines'][res_name][i]))/float(sas['181L']['amines'][res_name][i]))
        data.append(round(sum(my_perc_diff)/len(my_perc_diff), 4) * 100)


    data1 = data[::2]
    data2 = data[1::2]
    ymax = max(data)
    # Bar width
    bar_width = 0.35

    # Index for the x-axis
    x = range(len(labels))

    # Create the bar graph
    plt.bar(x, data1, width=bar_width, label='Additively Weighted')
    plt.bar([i + bar_width for i in x], data2, width=bar_width, label='Power')

    # Add labels and title
    plt.ylabel('% Difference Surface Area')
    plt.title('181L Average Residue % Difference (SA)')

    # Angle the labels and add values at the top of the bars
    plt.xticks([i + bar_width / 2 for i in x], labels, rotation=45, ha='right')
    for i, v in enumerate(data1):
        plt.text(i, v / 2, str(round(v, 2)) + ' %', ha='center', va='center', rotation=90)
    for i, v in enumerate(data2):
        plt.text(i + bar_width, v / 2, str(round(v, 2)) + ' %', ha='center', va='center', rotation=90)
    plt.ylim(0, 1.25 * ymax)
    # Add legend with appropriate layout
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 0.97), shadow=True, ncol=2)


    # Show the plot
    plt.tight_layout()
    plt.show()