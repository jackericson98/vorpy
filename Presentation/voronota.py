from objects import Vertex


# Read verts function. Used to interpret Voronota data
def read_verts(file):
    file = open(file).readlines()
    verts = []

    for i in range(len(file)):
        data = file[i].split(" ")
        verts.append(Vertex(data[3:6], data[6]))
        return verts