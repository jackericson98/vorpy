from Data.Analyze.compare.compare_files import compare_files
import csv


if __name__ == '__main__':
    prefix1 = 'C:/Users/jacke/Documents/data/'
    prefix = 'C:/Users/i7-8700/Documents/logs_pdbs/'
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

