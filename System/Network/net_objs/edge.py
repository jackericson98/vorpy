from System.sys_funcs.calcs import calc_circ, calc_angle
from System.Network.net_objs.surface import Surface
import numpy as np


class Edge:
    """Edge object. Used to build the network and calculate the surfaces"""
    def __init__(self, atoms=None, net=None, verts=None, surfs=None, doublet=False, points=None, center=None, rad=None,
                 normal=None, dist=None, pv0=None, pv1=None, ndx=None, load_ndxs=None, point_refs=None, straight=False):

        # If no network was given have a catch
        if net is not None and net.atoms is not None and atoms is not None:
            ndx = [atom.num for atom in atoms]
            ndx.sort()
        self.ndx = ndx                   # Index         :   Indices of the atoms of the surface
        self.net = net                   # Network       :   Network of the System
        self.atoms = atoms               # Atoms         :   List of Atom type objects for the edge
        self.verts = verts               # Vertices      :   List of Vertex type objects
        self.surfs = surfs               # Surfaces      :   List of 2 surfaces attached to the edge
        self.load_ndxs = load_ndxs       # Load indices  :   List of object load indices
        self.point_refs = point_refs     # Point refs    :   List of the surface and indices of the points for the edge

        self.loc = center                # Location      :   Location of the center of the 3 atoms that make up the edge
        self.rad = rad                   # Radius        :   Radius of the inscribed circle of the three atoms
        self.norm = normal
        self.dist = dist
        self.points = points             # Points        :   List of points along the
        self.pv0 = pv0                   # Vertex pt 0   :   The points on the ends of the edges
        self.pv1 = pv1                   # Vertex pt 1   :   The points on the ends of the edges
        self.pa = None                   # Projection pt :   The projection point from which the edge is built
        self.doublet = doublet           # Doublet       :   Boolean for if the edge is part of a doublet or not
        self.loc2 = None                 # Loc2          :   Allows edges to be checked like vertices
        self.draw_points = None
        self.draw_tris = None
        self.ref = None                  # Reference     :   Tuple holding a surface and a range for efficient storage
        self.straight = straight         # Straight edge :   Straight edge or not
