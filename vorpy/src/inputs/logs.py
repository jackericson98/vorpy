import csv
import pandas as pd


def parse_string_lists_int(string_list, apply_type=int):
    # Test if it is just one single list
    if string_list[1] != '[':
        listy, current_number = [], ''
        for letter in string_list:
            if letter.isdigit() or letter == '.':
                current_number += letter
            elif letter == ',':
                listy.append(apply_type(current_number))
                current_number = ''
    else:
        listy, current_number = [[]], ''
        for letter in string_list:
            if letter.isdigit() or letter == '.':
                current_number += letter
            elif letter == ',' and len(current_number) > 0:
                listy[-1].append(apply_type(current_number))
                current_number = ''
            elif letter == ']':
                listy.append([])
                current_number = ''
    return listy


def parse_string_lists(string_list, apply_type=float):
    # Check if the string list contains np.float64 values
    if 'n' in string_list and 'p' in string_list:
        remove_float_vals = True
    else:
        remove_float_vals = False
    # Test if it is just one single list
    if string_list[1] != '[':
        listy, current_number = [], ''
        for letter in string_list:
            if letter.isdigit() or letter == '.':
                current_number += letter
            elif letter == ',':
                if remove_float_vals:
                    current_number = current_number[2:]
                listy.append(apply_type(current_number))
                current_number = ''
    else:
        listy, current_number = [[]], ''
        for letter in string_list:
            if letter.isdigit() or letter == '.':
                current_number += letter
            elif letter == ',' and len(current_number) > 0:
                if remove_float_vals:
                    current_number = current_number[2:]
                listy[-1].append(apply_type(current_number))
                current_number = ''
            elif letter == ']':
                listy.append([])
                current_number = ''
    return listy


def sort_bool(stringy):
    return True if stringy == 'True' else False


atom_vals = {'Index': int, 'Name': str, 'Residue': str, 'Residue Sequence': int, 'Chain': str, 'Mass': float,
             'X': float, 'Y': float, 'Z': float, 'Radius': float, 'Volume': float, 'Van Der Waals Volume': float,
             'Surface Area': float, 'Complete Cell?': sort_bool, 'Maximum Mean Curvature': float,
             'Average Mean Surface Curvature': float, 'Maximum Gaussian Curvature': float,
             'Average Gaussian Surface Curvature': float, 'Sphericity': float, 'Isometric Quotient': float,
             'Inner Ball?': sort_bool, 'Number of Neighbors': int, 'Closest Neighbor': int,
             'Closest Neighbor Distance': float, 'Layer Distance Average': parse_string_lists,
             'Layer Distance RMSD': parse_string_lists, 'Minimum Point Distance': float,
             'Maximum Point Distance': float, 'Number of Overlaps': int, 'Contact Area': float,
             'Non - Overlap Volume': float, 'Overlap Volume': float, 'Center of Mass': parse_string_lists,
             'Moment of Inertia Tensor': parse_string_lists, 'Bounding Box': parse_string_lists,
             'Neighbors': parse_string_lists_int}


atom_vals_old = {'Index': int, 'Name': str, 'Residue': str, 'Residue Sequence': int, 'Chain': str, 'Mass': float,
                 'X': float, 'Y': float, 'Z': float, 'Radius': float, 'Volume': float, 'Van Der Waals Volume': float,
                 'Surface Area': float, 'Complete Cell?': sort_bool, 'Maximum Curvature': float,
                 'Average Surface Curvature': float, 'Sphericity': float, 'Isometric Quotient': float,
                 'Inner Ball?': sort_bool, 'Number of Neighbors': int, 'Closest Neighbor': int,
                 'Closest Neighbor Distance': float, 'Layer Distance Average': parse_string_lists,
                 'Layer Distance RMSD': parse_string_lists, 'Minimum Point Distance': float,
                 'Maximum Point Distance': float, 'Number of Overlaps': int, 'Contact Area': float,
                 'Non - Overlap Volume': float, 'Overlap Volume': float, 'Center of Mass': parse_string_lists,
                 'Moment of Inertia Tensor': parse_string_lists, 'Bounding Box': parse_string_lists,
                 'Neighbors': parse_string_lists_int}

atom_vals_new = {'Index': int, 'Name': str, 'Residue': str, 'Residue Sequence': int, 'Chain': str, 'Mass': float,
                 'X': float, 'Y': float, 'Z': float, 'Radius': float, 'Volume': float, 'Van Der Waals Volume': float,
                 'Surface Area': float, 'Complete Cell?': sort_bool, 'Maximum Mean Curvature': float,
                 'Average Mean Surface Curvature': float, 'Maximum Gaussian Curvature': float,
                 'Average Gaussian Surface Curvature': float, 'Integrated Mean Curvature': float,
                 'Integrated Mean Curvature Squared': float, 'Integrated Gaussian Curvature': float,
                 'Sphericity': float, 'Isometric Quotient': float, 'Inner Ball?': sort_bool, 'Number of Neighbors': int,
                 'Closest Neighbor': int, 'Closest Neighbor Distance': float, 'Layer Distance Average': parse_string_lists,
                 'Layer Distance RMSD': parse_string_lists, 'Minimum Point Distance': float,
                 'Maximum Point Distance': float, 'Number of Overlaps': int, 'Contact Area': float,
                 'Non - Overlap Volume': float, 'Overlap Volume': float, 'Center of Mass': parse_string_lists,
                 'Moment of Inertia Tensor': parse_string_lists, 'Bounding Box': parse_string_lists,
                 'Neighbors': parse_string_lists_int}


def read_atom(atom_line):
    """
    Parse an atom row using its column count.

    Supported formats
    -----------------
    36 columns
        Existing VorPy atom log format.

    39 columns
        New VorPy atom log format containing:
        - Integrated Mean Curvature
        - Integrated Mean Curvature Squared
        - Integrated Gaussian Curvature
    """

    n = len(atom_line)

    # --------------------------------------------------------------
    # New format: 39 columns
    # --------------------------------------------------------------
    if n == 39:
        schema = atom_vals_new

    # --------------------------------------------------------------
    # Current/legacy format: 36 columns
    # --------------------------------------------------------------
    elif n == 36:
        schema = atom_vals

    # --------------------------------------------------------------
    # Older legacy format
    # --------------------------------------------------------------
    elif n == len(atom_vals_old):
        schema = atom_vals_old

    else:
        raise ValueError(
            f"Unrecognized atom log format: "
            f"expected 39, 36, or {len(atom_vals_old)} columns, "
            f"got {n}."
        )

    atom = {}

    for i, title in enumerate(schema):
        atom[title] = schema[title](atom_line[i])

    return atom


def read_surf(surf_line):
    """
    Parse a surface row based on its column count.

    Supported formats
    -----------------
    15 columns
        New VorPy format with integrated curvature values.

    12 columns
        Current VorPy format.

    10 columns
        Older mean/Gaussian curvature format.

    9 columns
        Older format without overlap.

    8 columns
        Legacy single-curvature format.
    """

    # Remove trailing empty CSV entries
    while surf_line and surf_line[-1] == "":
        surf_line = surf_line[:-1]

    n = len(surf_line)

    # --------------------------------------------------------------
    # New format: 15 columns
    # --------------------------------------------------------------
    if n == 15:
        return {
            "Index": int(surf_line[0]),
            "Balls": [
                int(surf_line[1]),
                int(surf_line[2])
            ],

            "Surface Area": float(surf_line[3]),

            # Keep existing names unchanged
            "Mean Curvature": float(surf_line[4]),
            "Average Mean Curvature": float(surf_line[5]),
            "Gauss Curvature": float(surf_line[6]),
            "Average Gauss Curvature": float(surf_line[7]),

            # New values
            "Integrated Mean Curvature": float(surf_line[8]),
            "Integrated Mean Curvature Squared": float(surf_line[9]),
            "Integrated Gaussian Curvature": float(surf_line[10]),

            "Ball Volumes": [
                float(surf_line[11]),
                float(surf_line[12])
            ],

            "Contact Area": float(surf_line[13]),
            "Overlap": float(surf_line[14]),
        }

    # --------------------------------------------------------------
    # Current format: 12 columns
    # --------------------------------------------------------------
    elif n == 12:
        return {
            "Index": int(surf_line[0]),
            "Balls": [
                int(surf_line[1]),
                int(surf_line[2])
            ],

            "Surface Area": float(surf_line[3]),

            "Mean Curvature": float(surf_line[4]),
            "Average Mean Curvature": float(surf_line[5]),
            "Gauss Curvature": float(surf_line[6]),
            "Average Gauss Curvature": float(surf_line[7]),

            "Ball Volumes": [
                float(surf_line[8]),
                float(surf_line[9])
            ],

            "Contact Area": float(surf_line[10]),
            "Overlap": float(surf_line[11]),
        }

    # --------------------------------------------------------------
    # Older format: 10 columns
    # --------------------------------------------------------------
    elif n == 10:
        return {
            "Index": int(surf_line[0]),
            "Balls": [
                int(surf_line[1]),
                int(surf_line[2])
            ],

            "Surface Area": float(surf_line[3]),

            "Mean Curvature": float(surf_line[4]),
            "Gauss Curvature": float(surf_line[5]),

            "Ball Volumes": [
                float(surf_line[6]),
                float(surf_line[7])
            ],

            "Contact Area": float(surf_line[8]),
            "Overlap": float(surf_line[9]),
        }

    # --------------------------------------------------------------
    # Older format: 9 columns
    # --------------------------------------------------------------
    elif n == 9:
        return {
            "Index": int(surf_line[0]),
            "Balls": [
                int(surf_line[1]),
                int(surf_line[2])
            ],

            "Surface Area": float(surf_line[3]),

            "Mean Curvature": float(surf_line[4]),
            "Gauss Curvature": float(surf_line[5]),

            "Ball Volumes": [
                float(surf_line[6]),
                float(surf_line[7])
            ],

            "Contact Area": float(surf_line[8]),
            "Overlap": 0.0,
        }

    # --------------------------------------------------------------
    # Legacy format: 8 columns
    # --------------------------------------------------------------
    elif n == 8:
        return {
            "Index": int(surf_line[0]),
            "Balls": [
                int(surf_line[1]),
                int(surf_line[2])
            ],

            "Surface Area": float(surf_line[3]),
            "Curvature": float(surf_line[4]),

            "Ball Volumes": [
                float(surf_line[5]),
                float(surf_line[6])
            ],

            "Contact Area": float(surf_line[7]),
            "Overlap": 0.0,
        }

    raise ValueError(
        f"Unrecognized surface log format: "
        f"expected 15, 12, 10, 9, or 8 columns, got {n}."
    )


def read_edge(edge_line):
    edge = {'Index': int(edge_line[0]), 'Balls': [int(_) for _ in edge_line[1:4]], 'Length': float(edge_line[4])}
    return edge


def read_vert(vert_line):
    vert = {'Index': int(vert_line[0]), 'Balls': [int(_) for _ in vert_line[1:5]],
            'loc': [float(_) for _ in vert_line[5:8]], 'rad': float(vert_line[8])}
    return vert


def read_logs(log_files, return_dict=False, no_sol=False, all_=True, balls=False, surfs=False, edges=False, verts=False):
    # Set up the dictionary to store the data
    file_info = {}
    # Create the one_file variable to track this. 
    one_file = False
    # Figure out if the log_files is a single file or a list of files and change the variable accordingly
    if type(log_files) is str:
        one_file = True
        log_files = [log_files]
    # Loop through the files
    for file in log_files:
        # Open the file
        with open(file, 'r') as logs:
            # Create the log reader
            log_reader = csv.reader(logs)
            # Set up the data type
            data_type = 'data'
            # Set up the lists to store the data
            atoms, surf_list, edge_list, vert_list = [], [], [], []
            # Set up the skip_next variable to track if the next line should be skipped
            skip_next = False
            # Loop through the lines
            for i, line in enumerate(log_reader):
                # Skip the first, the second, the fourth, and the fifth lines
                if i in {0, 1, 3, 4}:
                    continue
                # Get the main data from the logs file. 
                elif i == 2:
                    line = line + [0 for _ in range(11 - len(line))]
                    # Try to get the data from the line
                    try:
                        data = {'name': line[0], 'network_type': line[1], 'surface_resolution': float(line[2]),
                                'box_size': float(line[3]), 'max_vert': float(line[4]), 'Total_Time': float(line[5]),
                                'vert_time': float(line[6]), 'connect_time': float(line[7]), 'surf_time': float(line[8]),
                                'analysis_time': float(line[9]), 'max_vertex': float(line[10])}
                        continue
                    # If the data is not found, get the data from the new logs type
                    except ValueError:
                        data = {'name': line[0], 'location': line[1], 'time': line[2], 'network_type': line[3],
                                'surface_resolution': float(line[4]), 'box_size': float(line[5]),
                                'max_vert': float(line[6]), 'Total_Time': float(line[7]),
                                'vert_time': float(line[8]), 'connect_time': float(line[9]),
                                'surf_time': float(line[10]), 'analysis_time': float(line[11]),
                                'max_vertex': float(line[12]), 'version': '< 3.2.0' if len(line) <= 13 else line[13]}
                        continue
                # Get the group data
                elif i == 5:
                    group_data = {'Name': line[0], 'Volume': float(line[1]), 'Surface Area': float(line[2]),
                                  'Mass': float(line[3]), 'Density': float(line[4]),
                                  'Center of Mass': parse_string_lists(line[5]), 'VDW Volume': float(line[6]),
                                  'VDW Center of Mass': parse_string_lists(line[7]),
                                  'Moment of Inertia': parse_string_lists(line[8]),
                                  'Spatial Moment of Inertia': parse_string_lists(line[9])}
                    continue

                # If the line is a build information, group information, Atoms, Edges, Surfaces, or Vertices, set the
                # data type and skip the next line
                if line[0] in {'build information', 'group information', 'Atoms', 'Edges', 'Surfaces', 'Vertices'}:
                    data_type = line[0]
                    skip_next = True
                    continue

                # If the skip_next variable is True, skip the next line
                if skip_next:
                    skip_next = False
                    continue
                # If the data type is Atoms and the all_ or balls variable is True, read the atom data
                elif data_type == 'Atoms' and (all_ or balls):
                    my_atom = read_atom(line)
                    my_atom['rad'], my_atom['loc'] = my_atom['Radius'], [my_atom['X'], my_atom['Y'], my_atom['Z']]
                    if no_sol and my_atom['name'].strip().lower() in {'hw1', 'hw2', 'ow', 'h02', 'h01', 'na', 'cl', 'mg', 'k'}:
                        continue
                    else:
                        atoms.append(my_atom)
                # If the data type is Surfaces and the all_ or surfs variable is True, read the surface data
                elif data_type == 'Surfaces' and (all_ or surfs):
                    surf_list.append(read_surf(line))
                # If the data type is Edges and the all_ or edges variable is True, read the edge data
                elif data_type == 'Edges' and (all_ or edges):
                    edge_list.append(read_edge(line))
                # If the data type is Vertices and the all_ or verts variable is True, read the vertex data
                elif data_type == 'Vertices' and (all_ or verts):
                    vert_list.append(read_vert(line))
                # If the data type is not one of the above, skip the line
                else:
                    continue
            # Get the file name
            file_name = data['name']
            # Set up the index to track the number of files with the same name
            index = 0
            # Loop through the file name
            while True:
                # If the file name is already in the dictionary, add a number to the end of the file name
                if file_name in file_info:
                    if index != 0:
                        file_name = file_name[:-len(str(index - 1))]
                    file_name = file_name + str(index)
                else:
                    break
                index += 1
            # If the return_dict variable is True, add the data to the dictionary
            if return_dict:
                file_info[file_name] = {'data': data, 'group data': group_data, 'atoms': atoms, 'surfs': surf_list,
                                        'edges': edge_list, 'verts': vert_list}
            else:
                file_info[file_name] = {'data': data, 'group data': group_data, 'atoms': pd.DataFrame(atoms),
                                        'surfs': pd.DataFrame(surf_list), 'edges': pd.DataFrame(edge_list),
                                        'verts': pd.DataFrame(vert_list)}
    # If the one_file variable is True, return the first file in the dictionary
    if one_file:
        # Get the first file in the dictionary
        my_file = [_ for _ in file_info][0]
        return file_info[my_file]
    # If the one_file variable is False, return the dictionary
    return file_info
