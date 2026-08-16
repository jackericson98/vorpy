from numba import jit
import numpy as np
from vorpy.src.chemistry import special_radii


def global_vars(sub_boxes, my_box_verts, my_num_splits, my_max_ball_rad, my_sub_box_size):
    """Set global variables for box searching and ball sorting operations.

    This function initializes the global variables used throughout the sorting module
    for efficient spatial partitioning and ball retrieval operations. These variables
    define the spatial grid structure and parameters used for organizing and accessing
    balls in 3D space.

    Parameters
    ----------
    sub_boxes : numpy.ndarray
        Matrix containing the spatial partitioning of balls into sub-boxes
    my_box_verts : list of numpy.ndarray
        Vertices defining the bounding box of the system
    my_num_splits : int
        Number of divisions along each axis for the spatial grid
    my_max_ball_rad : float
        Maximum radius of any ball in the system
    my_sub_box_size : list of float
        Size of each sub-box along each axis
    """
    # Initialize global variables
    global balls_matrix, box_verts, num_splits, max_ball_rad, sub_box_size
    
    # Assign the input values to the global variables
    balls_matrix = sub_boxes
    box_verts = my_box_verts
    num_splits = my_num_splits
    max_ball_rad = my_max_ball_rad
    sub_box_size = my_sub_box_size


@jit(nopython=True)
def box_search_numba(loc, num_splits, box_verts):
    """Find the sub box indices for a given location in 3D space.

    This function calculates which sub-box a given point belongs to within a larger
    bounding box that has been divided into a grid of smaller sub-boxes. The function
    is optimized with Numba for performance.

    Parameters
    ----------
    loc : numpy.ndarray
        The 3D coordinates [x, y, z] of the point to locate
    num_splits : int
        Number of divisions along each axis of the bounding box
    box_verts : numpy.ndarray
        The vertices of the bounding box, where box_verts[0] is the minimum point
        and box_verts[1] is the maximum point

    Returns
    -------
    list or None
        A list of three integers [i, j, k] representing the sub-box indices,
        or None if the point lies outside the bounding box
    """
    # Calculate the size of the sub boxes
    sub_box_size = [round((box_verts[1][i] - box_verts[0][i]) / num_splits, 3) for i in range(3)]
    # Find the sub box for the ball
    box_ndxs = [int((loc[j] - box_verts[0][j]) / sub_box_size[j]) for j in range(3)]
    if box_ndxs[0] >= num_splits or box_ndxs[1] >= num_splits or box_ndxs[2] >= num_splits:
        return
    # Return the box indices
    return box_ndxs


def box_search(loc):
    """Locate the sub box indices for a given location in 3D space.

    This function serves as a wrapper for the Numba-optimized box_search_numba function,
    converting the input location to a numpy array and using the global variables
    num_splits and box_verts to determine which sub-box a point belongs to.

    Parameters
    ----------
    loc : array-like
        The 3D coordinates [x, y, z] of the point to locate

    Returns
    -------
    list or None
        A list of three integers [i, j, k] representing the sub-box indices,
        or None if the point lies outside the bounding box
    """
    return box_search_numba(np.array(loc), num_splits, np.array(box_verts))


def get_balls(cells, dist=0, cell_reach=0, my_balls_matrix=None, my_sub_box_size=None, my_max_ball_rad=None):
    """Retrieves a list of balls from a 3D grid of cells based on specified search parameters.

    This function searches for balls within a specified distance of given cells in a 3D grid.
    It expands the search area by a configurable number of cells and returns all balls found
    within the expanded search region.

    Parameters
    ----------
    cells : list of list of int or list of int
        The initial set of cells to search from. Can be a single cell [i,j,k] or a list of cells.
    dist : float, optional
        The distance to expand the search region from the initial cells (default: 0)
    cell_reach : int, optional
        Additional number of cells to expand the search region (default: 0)
    my_balls_matrix : numpy.ndarray, optional
        Custom ball matrix to search through (default: None)
    my_sub_box_size : float, optional
        Custom sub-box size for distance calculations (default: None)
    my_max_ball_rad : float, optional
        Custom maximum ball radius for search optimization (default: None)

    Returns
    -------
    list or None
        A list of balls found within the search region, or None if the input cells are invalid
    """
    # Get the universal variables
    global balls_matrix, sub_box_size, max_ball_rad
    # If cells is none there is an issue with the box search and we need to return none
    if cells is None:
        return
    # If the three variables are not specified set them equal to the globals
    if my_balls_matrix is not None:
        balls_matrix, sub_box_size, max_ball_rad = my_balls_matrix, my_sub_box_size, my_max_ball_rad
    # Get the reach around the box to grab balls from
    reach = int(dist / min(sub_box_size)) + 3
    # Grab the number of cells in the grid
    n = balls_matrix[-1, -1, -1][0]
    # If a single cell is entered
    if type(cells[0]) is int:
        cells = [cells]
    # Get the min and max of the cells
    ndx_min = list(cells[0])
    ndx_max = list(cells[0])

    for cell in cells[1:]:
        for i in range(3):
            if cell[i] < ndx_min[i]:
                ndx_min[i] = cell[i]
            if cell[i] > ndx_max[i]:
                ndx_max[i] = cell[i]

    # Get the bounded range of cells to search
    xs = range(max(0, -reach + ndx_min[0] - cell_reach), min(n, reach + ndx_max[0] + cell_reach))
    ys = range(max(0, -reach + ndx_min[1] - cell_reach), min(n, reach + ndx_max[1] + cell_reach))
    zs = range(max(0, -reach + ndx_min[2] - cell_reach), min(n, reach + ndx_max[2] + cell_reach))

    # Get the balls
    balls = []
    for i in xs:
        for j in ys:
            for k in zs:
                try:
                    balls += balls_matrix[i, j, k]
                except KeyError:
                    pass

    return balls


def ndx_search(ndxs_list, ndxs):
    """Performs a binary search on a sorted list of ball indices to find the insertion point for a new vertex.

    This function implements a binary search algorithm to efficiently locate where a new vertex index should be inserted
    into a sorted list of ball indices. The list is maintained in ascending order based on ball size.

    Parameters
    ----------
    ndxs_list : list
        A sorted list of ball indices to search through
    ndxs : int
        The new vertex index to find the insertion point for

    Returns
    -------
    int
        The index position where the new vertex should be inserted to maintain sorted order
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


def divide_box(net_box, divisions, c=0):
    """Divides a bounding box into smaller sub-boxes based on the number of divisions.

    This function takes a bounding box and recursively divides it into smaller sub-boxes
    based on the specified number of divisions. The division process prioritizes splitting
    along the longest dimension first to maintain balanced sub-box sizes.

    Parameters
    ----------
    net_box : list of lists
        The bounding box to divide, represented as [[x_min, y_min, z_min], [x_max, y_max, z_max]]
    divisions : int
        The desired number of divisions to create
    c : float, optional
        A small constant used to adjust box boundaries to prevent edge cases (default: 0)

    Returns
    -------
    list of lists
        A list of sub-boxes, where each sub-box is represented as [[x_min, y_min, z_min], [x_max, y_max, z_max]]
    """
    # Convert the divisions to two_pow
    two_pow = 0
    # Loop until the number of divisions is reached
    while True:
        # Define the polynomial function
        def poly(x):
            return 0.03704228 * x ** 3 + 0.33267327 * x ** 2 + 0.94711614 * x + 0.65148515
        # Get the number of divisions
        my_divs = poly(two_pow)
        # If the number of divisions is reached break
        if my_divs >= divisions:
            break
        # Increment the number of divisions
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
        # If the first dimension is not divided
        if xyz_divs[0] == 0:
            my_sub_boxes = [[[xs - c, ys - c, zs - c], [xe + c, ym + c, zm + c]],
                            [[xs - c, ym - c, zs - c], [xe + c, ye + c, zm + c]],
                            [[xs - c, ys - c, zm - c], [xe + c, ym + c, ze + c]],
                            [[xs - c, ym - c, zm - c], [xe + c, ye + c, ze + c]]]
        # If the second dimension is not divided
        elif xyz_divs[1] == 0:
            my_sub_boxes = [[[xs - c, ys - c, zs - c], [xm + c, ye + c, zm + c]],
                            [[xm - c, ys - c, zs - c], [xe + c, ye + c, zm + c]],
                            [[xs - c, ys - c, zm - c], [xm + c, ye + c, ze + c]],
                            [[xm - c, ys - c, zm - c], [xe + c, ye + c, ze + c]]]
        # If the third dimension is not divided
        elif xyz_divs[2] == 0:
            my_sub_boxes = [[[xs - c, ys - c, zs - c], [xm + c, ym + c, ze + c]],
                            [[xm - c, ys - c, zs - c], [xe + c, ym + c, ze + c]],
                            [[xs - c, ym - c, zs - c], [xm + c, ye + c, ze + c]],
                            [[xm - c, ym - c, zs - c], [xe + c, ye + c, ze + c]]]
        return my_sub_boxes
    # Create the subnets
    for i in range(xyz_divs[0] + 1):
        # If the first dimension is not divided
        for j in range(xyz_divs[1] + 1):
            # If the second dimension is not divided
            for k in range(xyz_divs[2] + 1):
                # Create the vertices for the sub net
                my_sub_boxes.append([[net_box[0][0] + i * dims[0] / (xyz_divs[0] + 1) - c,
                                      net_box[0][1] + j * dims[1] / (xyz_divs[1] + 1) - c,
                                      net_box[0][2] + k * dims[2] / (xyz_divs[2] + 1) - c],
                                     [net_box[0][0] + (i + 1) * dims[0] / (xyz_divs[0] + 1) + c,
                                      net_box[0][1] + (j + 1) * dims[1] / (xyz_divs[1] + 1) + c,
                                      net_box[0][2] + (k + 1) * dims[2] / (xyz_divs[2] + 1) + c]])
    # Return the sub boxes
    return my_sub_boxes


def get_sys_type(my_sys):
    """Determines the type of molecular system based on its composition.

    This function analyzes the residues in a molecular system to classify it as one of:
    - 'Protein': Contains only protein residues
    - 'Nucleic': Contains only nucleic acid residues
    - 'Complex': Contains both protein and nucleic acid residues
    - 'Molecule': Default type if no residues are present or if residues don't match known types

    Parameters
    ----------
    my_sys : object
        A molecular system object containing a residues attribute

    Returns
    -------
    str
        A string indicating the system type: 'Protein', 'Nucleic', 'Complex', or 'Molecule'
    """
    # Get the type of system it is (poly, nucleic, both, other)
    sys_type = 'Molecule'
    nucs = {'T', 'DT', 'G', 'DG', 'A', 'DA', 'C', 'DC', 'U', 'DU'}
    if len(my_sys.residues) > 0:
        for res in my_sys.residues:
            if res.name in special_radii:
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


def sort_lists(*lists, reverse=False):
    """Sorts multiple lists based on the values in the first list.

    This function takes multiple lists and sorts them all based on the values in the first list.
    The sorting maintains the relative order between elements across all lists. For example,
    if the first list is [3, 1, 2] and the second list is ['c', 'a', 'b'], after sorting
    the first list to [1, 2, 3], the second list will be sorted to ['a', 'b', 'c'].

    Parameters
    ----------
    *lists : list
        Variable number of lists to be sorted. All lists must have the same length.
    reverse : bool, optional
        If True, sorts in descending order. Default is False.

    Returns
    -------
    list of lists
        The sorted lists in the same order as input.

    Raises
    ------
    ValueError
        If the input lists have different lengths.
    """
    if not lists:
        return []

    # Ensure that all lists have the same length as the first list
    if not all(len(lst) == len(lists[0]) for lst in lists):
        raise ValueError("All lists must have the same length")

    # Combine the lists into a list of tuples for sorting
    combined = list(zip(*lists))

    # Sort based on the elements of the first list
    sorted_combined = sorted(combined, key=lambda x: x[0], reverse=reverse)

    # Unzip the sorted tuples back into separate lists
    return list(map(list, zip(*sorted_combined)))
