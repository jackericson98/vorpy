import pandas as pd
from System.Network.network import Network


def read_verts(file):

    verts = []
    with open(file, 'r') as my_file:
        for line in my_file.readlines()[1:-1]:
            line = line.split(' ')
            verts.append({'vatoms': [int(_) for _ in line[:4]], 'vloc': [float(_) for _ in line[4:7]], 'vrad': float(line[7])})
    verts = pd.DataFrame(verts)
    return verts

