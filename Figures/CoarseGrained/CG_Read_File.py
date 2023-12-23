import csv

with open("C:/Users/jacke/PycharmProjects/vorpy/Data/user_data/BSA_5ORF_CG_coarse/sys/BSA_5ORF_CG_coarse_logs_vor.csv", 'r') as cg_file:
    current_values = 'info'
    skip = False
    atoms, surfaces, edges, vertices = [], [], [], []
    for i, line in enumerate(cg_file.readlines()):
        # Skip the headers
        if skip:
            skip = False
            continue
        if line[0] in ['Atoms', 'Surfaces', 'Edges', 'Vertices']:
            current_values = {'a': 'Atoms', 's': 'Surfaces', 'e': 'Edges', 'v': 'Vertices'}[line[0]]
            skip = True
            continue
        if current_values == 'a':
            atoms.append({'index': line[0], 'name': line[1], 'volume': line[2], 'surface area': line[3],
                          'max curvature': line[4], 'neighbors': line[5:]})