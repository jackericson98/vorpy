import pandas as pd
from System.Network.network import Network


def read_verts(net, file):

    verts = []
    with open(file, 'r') as my_file:
        for i, line in enumerate(my_file.readlines()):
            line = line.split(' ')
            if i == 0:
                if len(line) == 15 and net.type.lower() != line[-1].lower()[:-1]:
                    print("\nWarning - Loaded Vertices do not match the set network type\n\n")
                try:
                    net.max_vert_rad = float(line[10][:-1])
                except ValueError:
                    pass

            if len(line) == 8:
                verts.append({'vatoms': [int(_) for _ in line[:4]], 'vloc': [float(_) for _ in line[4:7]], 'vrad': float(line[7])})
            elif len(line) == 9:
                verts.append({'vatoms': [int(_) for _ in line[:4]], 'vloc': [float(_) for _ in line[4:7]], 'vrad': float(line[7]), 'vdub': int(line[8])})
    verts = pd.DataFrame(verts)
    return verts

