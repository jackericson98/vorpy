from System.Network.connect_network import *


# Build network function. Takes in a system and returns a fully connected network
def build_network(sys):
    # Find the vertices of the system
    find_vertices(sys)
    # Connect the network of vertices
    connect_network(sys)
    # Return the system
    return sys.net
