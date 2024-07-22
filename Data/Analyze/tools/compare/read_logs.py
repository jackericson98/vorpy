import csv

import pandas as pd


def read_atom(atom_line):
    try:
        atom = {'num': int(atom_line[0]), 'name': atom_line[1], 'volume': float(atom_line[2]), 'sa': float(atom_line[3]),
                'max curv': float(atom_line[4]), 'neighbors': [int(_) for _ in atom_line[5:] if _ != '']}
    except ValueError:
        atom = {'num': int(atom_line[0]), 'name': atom_line[1], 'volume': float(atom_line[2]),
                'sa': float(atom_line[3]), 'max curv': float(atom_line[4]), 'complete': atom_line[5],
                'neighbors': [int(_) for _ in atom_line[6:] if _ != '']}
    return atom


def read_surf(surf_line):
    surf = {'index': int(surf_line[0]), 'atoms': [int(_) for _ in surf_line[1:3]], 'sa': float(surf_line[3]),
            'curvature': float(surf_line[4]), 'atom vols': [float(_) for _ in surf_line[5:] if _ != '']}
    return surf


def read_edge(edge_line):
    edge = {'index': int(edge_line[0]), 'atoms': [int(_) for _ in edge_line[1:4]], 'length': float(edge_line[4])}
    return edge


def read_vert(vert_line):
    vert = {'index': int(vert_line[0]), 'atoms': [int(_) for _ in vert_line[1:5]],
            'loc': [float(_) for _ in vert_line[5:8]], 'rad': float(vert_line[8])}
    return vert


#
def read_logs(log_files, return_dict=False, no_sol=False):
    file_info = {}
    one_file = False
    if type(log_files) is str:
        one_file = True
        log_files = [log_files]
    for file in log_files:
        with open(file, 'r') as logs:
            log_reader = csv.reader(logs)
            data_type = 'data'
            atoms, surfs, edges, verts = [], [], [], []
            skip_next = False
            for i, line in enumerate(log_reader):
                if i in {0, 1, 3, 4}:
                    continue
                elif i == 2:
                    line = line + [0 for _ in range(11 - len(line))]
                    data = {'name': line[0], 'network_type': line[1], 'surface_resolution': float(line[2]),
                            'box_size': float(line[3]), 'max_vert': float(line[4]), 'Total_Time': float(line[5]),
                            'vert_time': float(line[6]), 'connect_time': float(line[7]), 'surf_time': float(line[8]),
                            'analysis_time': float(line[9]), 'max_vertex': float(line[10])}
                    continue
                elif i == 5:
                    group_data = {'index': int(line[0]), 'name': line[1], 'volume': float(line[2]), 'sa': float(line[3])}
                    continue
                if line[0] in {'build information', 'group information', 'Atoms', 'Edges', 'Surfaces', 'Vertices'}:
                    data_type = line[0]
                    skip_next = True
                    continue
                if skip_next:
                    skip_next = False
                    continue
                elif data_type == 'Atoms':
                    my_atom = read_atom(line)
                    if no_sol and my_atom['name'].strip().lower() in {'hw1', 'hw2', 'ow', 'h02', 'h01', 'na', 'cl', 'mg', 'k'}:
                        continue
                    else:
                        atoms.append(my_atom)
                elif data_type == 'Surfaces':
                    surfs.append(read_surf(line))
                elif data_type == 'Edges':
                    edges.append(read_edge(line))
                elif data_type == 'Vertices':
                    verts.append(read_vert(line))
                else:
                    continue
            file_name = data['name']
            index = 0
            while True:
                if file_name in file_info:
                    if index != 0:
                        file_name = file_name[:-len(str(index - 1))]
                    file_name = file_name + str(index)
                else:
                    break
                index += 1
            if return_dict:
                file_info[file_name] = {'data': data, 'group data': group_data, 'atoms': atoms, 'surfs': surfs,
                                        'edges': edges, 'verts': verts}
            else:
                file_info[file_name] = {'data': data, 'group data': group_data, 'atoms': pd.DataFrame(atoms),
                                        'surfs': pd.DataFrame(surfs), 'edges': pd.DataFrame(edges),
                                        'verts': pd.DataFrame(verts)}
    if one_file:
        my_file = [_ for _ in file_info][0]
        return file_info[my_file]
    return file_info
