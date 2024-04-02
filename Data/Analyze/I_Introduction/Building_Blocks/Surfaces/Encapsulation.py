from numpy import array as ar
import matplotlib.pyplot as plt
from Visualize.mpl_visualize import plot_atoms, plot_verts, plot_edges, plot_surfs
from System.Network.verts.calc_vert import calc_vert, calc_flat_vert
from System.Network.edges.build_edge import build_edge
from System.Network.surfs.build_surf import build_surf


plot_atoms([[0, 0, -0.5], [0.1, 0.0, 0.2]], [1.0, 3.0], Show=True, alpha=0.3)
