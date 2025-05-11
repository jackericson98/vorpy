from pandas import DataFrame
from vorpy.src.network import Network


# Add Voronota data method. Takes in voronota data and adds it to the System
def read_vta(grp, ball_file, vert_file):
    """
    Reads Voronota data files and populates a group with ball and vertex information.

    This function processes two Voronota data files:
    - ball_file: Contains atom information including indices and properties
    - vert_file: Contains vertex information connecting atoms

    Parameters:
    -----------
    grp : Group
        The group object to populate with ball and vertex data
    ball_file : str
        Path to the Voronota ball data file
    vert_file : str
        Path to the Voronota vertex data file

    Returns:
    --------
    None
        Modifies the group object in place by:
        - Creating a Network object if one doesn't exist
        - Populating the network with vertex data from the files
        - Linking vertices to the appropriate balls in the system

    Notes:
    ------
    - Progress is displayed during file loading
    - Vertex data is sorted and deduplicated
    - Ball indices are adjusted to be 0-based
    """

    # Create the System and load the files
    with open(ball_file, 'r') as b, open(vert_file, 'r') as v:
        b_file, v_file = b.readlines(), v.readlines()
    # Create the ball and vert lists
    verts, balls = [], []
    # Loop through the ball file
    for i in range(len(b_file)):
        # Print the progress
        print("\rLoading Balls - {:.2f}%".format(100 * i/len(b_file)), end='')
        # Split the data
        data = b_file[i].split(" ")
        # Grab the data reference for the atoms
        balls.append(int(data[5]) - 1)
    # Interpret the vertices
    for i in range(len(v_file)):
        print("\rLoading verts - {:.2f}%".format(100 * i/len(v_file)), end='')
        # Split the data
        data = v_file[i].split(" ")
        # Add the vertex data
        loc, rad = [float(data[4]), float(data[5]), float(data[6])], float(data[7])
        # Create the atoms list
        atoms = [balls[int(data[0])], balls[int(data[1])], balls[int(data[2])], balls[int(data[3])]]
        # Sort the atoms
        atoms.sort()
        # Create the dub flag
        dub = 0
        # If the atoms are the same as the last vertex, set the dub flag to 1
        if i > 0 and atoms == verts[-1]['vatoms']:
            dub = 1
        # Add the vertex to the list
        verts.append({'vatoms': atoms, 'vloc': loc, 'vrad': rad, 'vdub': dub})
    # Check to see if there is anetwork associated with the group
    if grp.net is None:
        # Create the network
        grp.net = Network(locs=grp.sys.balls['loc'], rads=grp.sys.balls['rad'], group=grp.ball_ndxs,
                          settings=grp.settings, masses=grp.sys.balls['mass'])
    # Add the vertices to the network
    grp.net.verts = DataFrame(verts)
