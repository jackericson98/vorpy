from System.sys_funcs.calcs import *


class Atom:
    def __init__(self, system=None, location=None, radius=None, index='', name='', residue='', chain='', res_seq="",
                 ocp="", t_fact="", seg_id="", element="", charge="", bonds=None, chn=None, res=None):

        # System groups
        self.sys = system           # System       :   Main system object
        self.res = res              # Residue      :   Residue object of which the atom is a part
        self.chn = chn              # Chain        :   Chain object of which the atom is a part

        self.loc = location         # Location     :   Set the location of the center of the sphere
        self.rad = radius           # Radius       :   Set the radius for the sphere object. Default is 1

        # Calculated Traits
        self.vol = 0                # Cell Volume  :   Volume of the voronoi cell for the atom
        self.sa = 0                 # Surface Area :   Surface area of the atom's cell
        self.box = []               # Box          :   The grid location of the atom

        # Network objects
        self.verts = []             # Vertices     :   List of Vertex type objects
        self.surfs = []             # Surfaces     :   List of Surface type objects
        self.edges = []             # Edges        :   List of Edge type objects

        # Input traits
        self.num = index            # Number       :   The index from the initial atom file
        self.name = name            # Name         :   Name retrieved from pdb file
        self.chain = chain          # Chain        :   Molecule chain the atom is a part of
        self.residue = residue      # Residue      :   Class of molecule that the atom is a part of
        self.res_seq = res_seq      # Sequence     :   Sequence of the residue that the atom is a part of
        self.occupancy = ocp        # Occupancy    :   Occupancy of the atom
        self.t_fact = t_fact        # Temp Factor  :   Temperature factor for the atom
        self.seg_id = seg_id        # Segment ID   :   Segment identifier for the atom
        self.element = element      # Symbol       :   Element of the atom
        self.charge = charge        # Charge       :   Charge of the atom
        self.bonds = bonds          # Bonds        :   Bonds to other atoms

        self.get_radius()

    def calc_vol(self):
        # Create the volume variable
        vol = 0
        # Go through each surface on the atom
        for surf in self.surfs:
            # If the surface hasn't been constructed yet, construct it
            if surf.points is None or surf.tris is None:
                surf.build()
            self.sa += surf.sa
            # Check to see if the surface's volume has been calculated already
            if surf.vols[surf.ndx.index(self.num)] != 0:
                vol += surf.vols[surf.ndx.index(self.num)]
            else:
                # Calculate the volume of the
                for tri in surf.tris:
                    p0, p1, p2, p3 = self.loc, surf.points[tri[0]], surf.points[tri[1]], surf.points[tri[2]]
                    my_vol = calc_tetra_vol(p0, p1, p2, p3)
                    surf.vols[surf.ndx.index(self.num)] = my_vol
                    vol += my_vol
        # Return the volume
        self.vol = vol
        return vol

    def get_radius(self):
        """
            Finds the radius of the atom from the symbol or vice versa

        :return: The radius of the atom from the symbol or vice versa
        """
        radii = self.sys.radii
        # If indicated we return the symbol of atom that the radius indicates
        if self.element is None:
            # Check to see if the radius is in the system
            if self.rad in {radii[_] for _ in radii[1]}:
                self.element = radii[self.rad]
            else:
                # Get the closest atom to it
                min_diff = np.inf
                # Go through the radii in the system looking for the smallest difference
                for radius in radii:
                    if radii[radius] - self.rad < min_diff:
                        self.element = radii[radius]
        # If we have the type and just want the radius, keep scanning until we find the radius
        elif self.rad is None:
            self.rad = radii[self.element.lower()]
