from pandas import DataFrame
from System.Network.network import Network


# Add Voronota data method. Takes in voronota data and adds it to the System
def read_vta_data(grp, ball_file, vert_file):

    # Create the System and load the files
    with open(ball_file, 'r') as b, open(vert_file, 'r') as v:
        b_file, v_file = b.readlines(), v.readlines()
    # Create the ball and vert lists
    verts, balls = [], []
    for i in range(len(b_file)):
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
        atoms = [balls[int(data[0])], balls[int(data[1])], balls[int(data[2])], balls[int(data[3])]]
        atoms.sort()
        dub = 0
        if i > 0 and atoms == verts[-1]['vatoms']:
            dub = 1
        verts.append({'vatoms': atoms, 'vloc': loc, 'vrad': rad, 'vdub': dub})
    # Check to see if there is anetwork associated with the group
    if grp.net is None:
        grp.net = Network(locs=grp.sys.balls['loc'], rads=grp.sys.balls['rad'], group=grp.ball_ndxs,
                          settings=grp.settings)
    grp.net.verts = DataFrame(verts)
