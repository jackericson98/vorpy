from System.sys_funcs.output.atoms import make_pdb_line
from os import path


"""
Looking to take in a pdb file and a log file and output a new pdb with a bfactor representing the % difference of
particular values
"""

def color_pdb_by_res(pdb, values, output_pdb=None):
    pdb_dir = path.dirname(pdb)
    pdb_name = path.basename(pdb)
    if output_pdb is None:
        output_pdb = pdb_dir + pdb_name[:-3] +'_colorized.pdb'
    with open(pdb, 'r') as read_pdb, open(output_pdb, 'w') as write_pdb:
        for line in read_pdb:



