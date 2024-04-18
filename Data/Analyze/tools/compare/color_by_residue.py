from os import path



def color_pdb_by_res(pdb, values, output_pdb=None):
    """
    Takes in a pdb file and a set of values corresponding to each residue
    """
    pdb_dir = path.dirname(pdb)
    pdb_name = path.basename(pdb)
    if output_pdb is None:
        output_pdb = pdb_dir + pdb_name[:-4] + '_colorized.pdb'
    with open(pdb, 'r') as read_pdb, open(output_pdb, 'w') as write_pdb:
        for line in read_pdb:
            # Only read the lines with Atom in the front because that is all we care about
            if line[:4] == 'ATOM':
                # First get the residue name to clear out any
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
