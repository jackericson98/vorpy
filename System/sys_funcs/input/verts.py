import pandas as pd


def read_verts(group, file):

    verts = []
    with open(file, 'r') as my_file:
        for i, line in enumerate(my_file.readlines()):
            line = line.split(' ')
            if i == 0:
                if len(line) == 15 and group.settings['net_type'].lower() != line[-1].lower()[:-1]:
                    print("\nWarning - Loaded Vertices do not match the set network type\n\n")
                try:
                    group.settings['max_vert'] = float(line[9][:-1])
                except ValueError:
                    pass
                continue
            if line[0] == 'END':
                continue
            if int(line[8]) == 1 and [int(_) for _ in line[:4]] == verts[-1]['balls']:
                verts[-1]['loc2'] = [float(_) for _ in line[4:7]]
                verts[-1]['rad2'] = float(line[7])
            else:
                verts.append(
                    {'balls': [int(_) for _ in line[:4]], 'loc': [float(_) for _ in line[4:7]], 'rad': float(line[7]),
                     'loc2': None, 'rad2': None})

    verts = pd.DataFrame(verts)
    return verts

