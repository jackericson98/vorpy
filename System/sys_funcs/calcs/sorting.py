from numba import jit
import numpy as np


def global_vars(sub_boxes, my_box_verts, my_num_splits, my_max_ball_rad, my_sub_box_size):
    global balls_matrix, box_verts, num_splits, max_ball_rad, sub_box_size
    balls_matrix = sub_boxes
    box_verts = my_box_verts
    num_splits = my_num_splits
    max_ball_rad = my_max_ball_rad
    sub_box_size = my_sub_box_size


@jit(nopython=True)
def box_search_numba(loc, num_splits, box_verts):
    # Calculate the size of the sub boxes
    sub_box_size = [round((box_verts[1][i] - box_verts[0][i]) / num_splits, 3) for i in range(3)]
    # Find the sub box for the ball
    box_ndxs = [int((loc[j] - box_verts[0][j]) / sub_box_size[j]) for j in range(3)]
    if box_ndxs[0] >= num_splits or box_ndxs[1] >= num_splits or box_ndxs[2] >= num_splits:
        return
    # Return the box indices
    return box_ndxs


def box_search(loc):
    """
    Locates the sub box indices for a given location
    """
    return box_search_numba(np.array(loc), num_splits, np.array(box_verts))


def get_balls(cells, dist=0, cell_reach=0, my_balls_matrix=None, my_sub_box_size=None, my_max_ball_rad=None):
    """
    Takes in the cells and the number of additional cells to search and returns an ball list
    :param cells: The initial boxes in the network to stem from
    :param dist: The number of cells out from the initial set of cells to search
    """
    # Get the universal variables
    global balls_matrix, sub_box_size, max_ball_rad
    # If the three variables are not specified set them equal to the globals
    if my_balls_matrix is not None:
        balls_matrix, sub_box_size, max_ball_rad = my_balls_matrix, my_sub_box_size, my_max_ball_rad
    # Get the reach around the box to grab balls from
    reach = int(dist / min(sub_box_size)) + 2
    # Grab the number of cells in the grid
    n = balls_matrix[-1, -1, -1][0]
    # If a single cell is entered
    if type(cells[0]) is int:
        cells = [cells]
    # Get the min and max of the cells
    ndx_min = [np.inf, np.inf, np.inf]
    ndx_max = [-np.inf, -np.inf, -np.inf]
    # Go through the cells and set the minimum and maximum indexes for xyz for a rectangle containing the balls
    for cell in cells:
        # Check each xyz index to see if they are larger or smaller than the max or min
        for i in range(3):
            if cell[i] < ndx_min[i]:
                ndx_min[i] = cell[i]
            if cell[i] > ndx_max[i]:
                ndx_max[i] = cell[i]
    xs = [x for x in range(max(0, -reach + ndx_min[0] - cell_reach), reach + ndx_max[0] + cell_reach)]
    ys = [y for y in range(max(0, -reach + ndx_min[1] - cell_reach), reach + ndx_max[1] + cell_reach)]
    zs = [z for z in range(max(0, -reach + ndx_min[2] - cell_reach), reach + ndx_max[2] + cell_reach)]
    balls = []
    # Get balls
    for i in xs:
        if 0 <= i < n:
            for j in ys:
                if 0 <= j < n:
                    for k in zs:
                        if 0 <= k < n:
                            try:
                                balls += balls_matrix[i, j, k]
                            except KeyError:
                                pass
    return balls


def ndx_search(ndxs_list, ndxs):
    """
     Searches a list of indices of balls sorted by smallest ball and where the vertex would be
    :param ndxs_list: The index for checking
    :param ndxs: The indices to check against
    :return: The vertex index of the vertex or where the vertex should be inserted
    """
    # If the length of the test list is equal to 0 return the next index
    if len(ndxs_list) <= 1:
        # If there exists one vertex already and the new vertex is less than the old vertex return 1
        if len(ndxs_list) > 0 and ndxs > ndxs_list[0]:
            return 1
        # Otherwise, return 0
        return 0
    # Get the middle of the list of vertices
    mid_list_ndx = len(ndxs_list) // 2
    # If the search element (my_list) is greater than the test element (test_lol) search the lower half of test_lol
    if ndxs > ndxs_list[mid_list_ndx]:
        ndxs_ndx = ndx_search(ndxs_list[mid_list_ndx:], ndxs)
        return ndxs_ndx + mid_list_ndx
    # If the search element (my_list) is less than the test element (test_lol) search the upper half of test_lol
    elif ndxs < ndxs_list[mid_list_ndx]:
        ndxs_ndx = ndx_search(ndxs_list[:mid_list_ndx], ndxs)
        return ndxs_ndx
    # If the search element (my_list) is greater than the test element (test_lol) search the lower half of test_lol
    elif ndxs == ndxs_list[mid_list_ndx]:
        return mid_list_ndx


def get_radius(ball):
    """
    Finds the radius of the ball from the symbol or vice versa
    :return: The radius of the ball from the symbol or vice versa
    """
    radii, special_radii = ball['sys'].radii, ball['sys'].special_radii
    # Get the radius and the element from the name of the ball
    if ball['res'] is not None and ball['res'].name in special_radii:
        # Check if no ball name exists or its empty
        if ball['name'] is not None and ball['name'] != '':
            for i in range(len(ball['name'])):
                name = ball['name'][:-i]
                # Check the residue name
                if name in special_radii[ball['res'].name]:
                    ball['rad'] = special_radii[ball['res'].name][name]
    # If we have the type and just want the radius, keep scanning until we find the radius
    if ball['rad'] is None and ball['element'].lower() in radii:
        ball['rad'] = radii[ball['element'].lower()]
    # If indicated we return the symbol of ball that the radius indicates
    if ball['rad'] is None or ball['rad'] == 0:
        # Check to see if the radius is in the system
        if ball['rad'] in {radii[_] for _ in radii[1]}:
            ball['element'] = radii[ball['rad']]
        else:
            # Get the closest ball to it
            min_diff = np.inf
            # Go through the radii in the system looking for the smallest difference
            for radius in radii:
                if radii[radius] - ball['rad'] < min_diff:
                    ball['element'] = radii[radius]
    return ball['rad']


def divide_box(net_box, divisions, c=0):
    # Convert the divisions to two_pow
    two_pow = 0
    while True:
        def poly(x):
            return 0.03704228 * x ** 3 + 0.33267327 * x ** 2 + 0.94711614 * x + 0.65148515
        my_divs = poly(two_pow)
        if my_divs >= divisions:
            break
        two_pow += 1

    # Find the order of dimensional subdivisions
    dims = [abs(net_box[0][i] - net_box[1][i]) for i in range(3)]
    sorted_dims, sorted_dim_ndxs = zip(*sorted(zip(dims, [0, 1, 2]), key=lambda x: x[0], reverse=True))

    # Determines the number of divisions per dimension
    num_divs = [two_pow // 3 + (1 if two_pow % 3 > i else 0) for i in range(3)]

    # Create the list of sub boxes
    my_sub_boxes = []

    # Get the divisions
    _, xyz_divs = zip(*sorted(zip(sorted_dim_ndxs, num_divs), key=lambda x: x[0]))

    # If one division
    if two_pow == 1:
        if xyz_divs[0] == 1:
            my_sub_boxes = [[[net_box[0][0] - c, net_box[0][1] - c, net_box[0][2] - c],
                             [net_box[0][0] + dims[0] / 2 + c, net_box[1][1] + c, net_box[1][2] + c]],
                            [[net_box[0][0] + dims[0] / 2 - c, net_box[0][1] - c, net_box[0][2] - c],
                             [net_box[1][0] + c, net_box[1][1] + c, net_box[1][2] + c]]]
        elif xyz_divs[1] == 1:
            my_sub_boxes = [[[net_box[0][0] - c, net_box[0][1] - c, net_box[0][2] - c],
                             [net_box[1][0] + c, net_box[0][1] + dims[1] / 2 + c, net_box[1][2] + c]],
                            [[net_box[0][0] - c, net_box[0][1] + dims[1] / 2 - c, net_box[0][2] - c],
                             [net_box[1][0] + c, net_box[1][1] + c, net_box[1][2] + c]]]
        elif xyz_divs[2] == 1:
            my_sub_boxes = [[[net_box[0][0] - c, net_box[0][1] - c, net_box[0][2] - c],
                             [net_box[1][0] + c, net_box[1][1] + c, net_box[0][2] + dims[2] / 2 + c]],
                            [[net_box[0][0] - c, net_box[0][1] - c, net_box[0][2] + dims[2] / 2 - c],
                             [net_box[1][0] + c, net_box[1][1] + c, net_box[1][2] + c]]]
        return my_sub_boxes

    # If two divisions
    elif two_pow == 2:
        xs, ys, zs = net_box[0]
        xm, ym, zm = [net_box[0][i] + dims[i] / 2 for i in range(3)]
        xe, ye, ze = net_box[1]
        if xyz_divs[0] == 0:
            my_sub_boxes = [[[xs - c, ys - c, zs - c], [xe + c, ym + c, zm + c]],
                            [[xs - c, ym - c, zs - c], [xe + c, ye + c, zm + c]],
                            [[xs - c, ys - c, zm - c], [xe + c, ym + c, ze + c]],
                            [[xs - c, ym - c, zm - c], [xe + c, ye + c, ze + c]]]
        elif xyz_divs[1] == 0:
            my_sub_boxes = [[[xs - c, ys - c, zs - c], [xm + c, ye + c, zm + c]],
                            [[xm - c, ys - c, zs - c], [xe + c, ye + c, zm + c]],
                            [[xs - c, ys - c, zm - c], [xm + c, ye + c, ze + c]],
                            [[xm - c, ys - c, zm - c], [xe + c, ye + c, ze + c]]]
        elif xyz_divs[2] == 0:
            my_sub_boxes = [[[xs - c, ys - c, zs - c], [xm + c, ym + c, ze + c]],
                            [[xm - c, ys - c, zs - c], [xe + c, ym + c, ze + c]],
                            [[xs - c, ym - c, zs - c], [xm + c, ye + c, ze + c]],
                            [[xm - c, ym - c, zs - c], [xe + c, ye + c, ze + c]]]
        return my_sub_boxes
    # Create the subnets
    for i in range(xyz_divs[0] + 1):
        for j in range(xyz_divs[1] + 1):
            for k in range(xyz_divs[2] + 1):
                # Create the vertices for the sub net
                my_sub_boxes.append([[net_box[0][0] + i * dims[0] / (xyz_divs[0] + 1) - c,
                                      net_box[0][1] + j * dims[1] / (xyz_divs[1] + 1) - c,
                                      net_box[0][2] + k * dims[2] / (xyz_divs[2] + 1) - c],
                                     [net_box[0][0] + (i + 1) * dims[0] / (xyz_divs[0] + 1) + c,
                                      net_box[0][1] + (j + 1) * dims[1] / (xyz_divs[1] + 1) + c,
                                      net_box[0][2] + (k + 1) * dims[2] / (xyz_divs[2] + 1) + c]])
    return my_sub_boxes


def get_sys_type(my_sys):
    # Get the type of system it is (poly, nucleic, both, other)
    sys_type = 'Molecule'
    nucs = {'T', 'DT', 'G', 'DG', 'A', 'DA', 'C', 'DC', 'U', 'DU'}
    if len(my_sys.residues) > 0:
        for res in my_sys.residues:
            if res.name in my_sys.special_radii:
                if sys_type == 'Nucleic':
                    sys_type = 'Complex'
                    break
                sys_type = 'Protein'
            elif res.name in nucs:
                if sys_type == 'Protein':
                    sys_type = 'Complex'
                    break
                sys_type = 'Nucleic'
    return sys_type
