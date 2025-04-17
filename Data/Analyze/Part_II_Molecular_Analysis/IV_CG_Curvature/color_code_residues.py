from System.sys_funcs.output.atoms import make_pdb_line
from System.sys_funcs.input.balls.pdb import read_pdb_line
from os import path
import matplotlib.pyplot as plt
import csv
import numpy as np


"""
Looking to take in a pdb file and a log file and output a new pdb with a bfactor representing the % difference of
particular values
"""


# Values have to have a res and a res seq double dict
def color_pdb_by_res(pdb, values, output_pdb=None):
    pdb_dir = path.dirname(pdb)
    pdb_name = path.basename(pdb)
    if output_pdb is None:
        output_pdb = pdb_dir + pdb_name[:-4] +'_colorized.pdb'
    with open(pdb, 'r') as read_pdb, open(output_pdb, 'w') as write_pdb:
        for line in read_pdb:
            if line[:4] == 'ATOM':

                res = line[17:20].strip()
                if res == 'SOL' or res == 'CL':
                    write_pdb.write(line)
                    continue

                res_seq = line[22:26].strip()
                if res in values:
                    bfact = values[res][res_seq]*10
                else:
                    write_pdb.write(line)
                    continue
                new_line = line[:62] + '{:>1.2f}'.format(bfact) + line[66:]

                write_pdb.write(new_line)
            else:
                write_pdb.write(line)


if __name__ == '__main__':
    prefix = 'C:/Users/jacke/Documents/data1/'
    pdb_files = [prefix + '1BNA.pdb', prefix + '1BNA_coarse_ad.pdb', prefix + '1BNA_coarse_ad.pdb',
                 prefix + '1BNA_coarse_ncap.pdb', prefix + '1BNA_coarse_ncap.pdb', prefix + '1BNA_coarse_scbb_ad.pdb',
                 prefix + '1BNA_coarse_scbb_ad.pdb', prefix + '1BNA_coarse_scbb_ncap.pdb', prefix + '1BNA_coarse_scbb_ncap.pdb']

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
                vols[line[0]][line[1]][line[2]] = {}
                sas[line[0]][line[1]][line[2]] = {}
            # Specific residue
            vols[line[0]][line[1]][line[2]][line[3]] = line[4]
            sas[line[0]][line[1]][line[2]][line[3]] = line[5]

    vor_atom_vals_vol = {}
    vor_atom_vals_sa = {}
    vol_res_data = {}
    sa_res_data = {}
    res_names = []
    for file in vols:
        if file == '1BNA':
            vor_atom_vals_vol[file] = {}
            vor_atom_vals_sa[file] = {}
        else:
            vol_res_data[file] = {}
            sa_res_data[file] = {}
        for res_type in vols[file]:
            for res_name in vols[file][res_type]:
                if file == '1BNA':
                    vor_atom_vals_vol[file][res_name] = {}
                    vor_atom_vals_sa[file][res_name] = {}
                else:
                    vol_res_data[file][res_name] = {}
                    vol_res_data[file][res_name] = {}
                for res in vols[file][res_type][res_name]:
                    if file == '1BNA':
                        vor_atom_vals_vol[file][res_name][res] = vols[file][res_type][res_name][res]
                        vor_atom_vals_sa[file][res_name][res] = sas[file][res_type][res_name][res]
                    else:
                        vol_act = float(vor_atom_vals_vol['1BNA'][res_name][res])
                        vol_res_data[file][res_name][res] = abs(float(vols[file][res_type][res_name][res]) - vol_act) / vol_act
                        # sa_act = float(vor_atom_vals_sa['1BNA'][res_name][res])
                        # sa_res_data[file][res_name][res] = abs(float(sas[file][res_type][res_name][res]) - sa_act) / sa_act


    for i, file in enumerate(vols):
        if file == '1BNA':
            mypdb = pdb_files[i]
            continue

        color_pdb_by_res(mypdb, vol_res_data[file], path.dirname(mypdb) + '/' + file + '_colorized.pdb')

