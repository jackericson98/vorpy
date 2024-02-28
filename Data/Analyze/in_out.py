

def in_out_data(sys, logs, sa=False, vol=False, curv=False):
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
