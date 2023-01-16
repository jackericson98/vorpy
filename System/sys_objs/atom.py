from System.sys_funcs.calcs import *


class Atom:
    """
    Atom object class used to represent loaded atoms

    location: list
        set the location of the center of the sphere
    radius : float
        set the radius for the sphere object. Default is 1
    system : System object
        set the atom's system attribute
    element : str
        element of the atom
    chain : str
        molecule chain the atom is a part of
    residue : str
        residue of the molecule that the atom is a part of
    res_seq : int
        sequence of the residue that the atom is a part of
    name : str
        name retrieved from pdb file
    ocp : str
        Occupancy of the atom
    t_fact : str
        Temperature factor for the atom
    seg_id : str
        Segment identifier for the atom
    charge : float
        Charge of the atom
    verts : list
        Vertex objects connected to the atom
    surfs : list
        Surface objects connected to the atom
    edges : list
        Edge objects connected to the atom
    load_ndxs : list
        Holds the object indices for when the system is loaded back in
    cell_vol : float
        Volume of the voronoi cell for the atom
    box : list
        The grid location of the atom

    """
    def __init__(self, location=None, radius=None, system=None, element=None, chain=None, mol_class=None, residue=None,
                 molecule=None, res_seq=None, name=None, ocp=None, t_fact=None, seg_id=None, charge=None, load_ndxs=None, index=None):

        # Calculated Traits
        self.loc = location         # Location     :   Set the location of the center of the sphere
        self.rad = radius           # Radius       :   Set the radius for the sphere object. Default is 1
        self.vol = 0                # Cell Volume  :   Volume of the voronoi cell for the atom
        self.sa = 0                 # Surface Area :   Surface area of the atom's cell
        self.box = []               # Box          :   The grid location of the atom

        # Network connections
        self.sys = system           # System       :   Set the atom's system attribute
        self.resid = residue        # Residue      :   The residue of the atom
        self.mol = molecule         # Molecule     :   The molecule that the atom is a part of
        self.verts = []             # Vertices     :   List of Vertex type objects
        self.surfs = []             # Surfaces     :   List of Surface type objects
        self.edges = []             # Edges        :   List of Edge type objects
        self.load_ndxs = load_ndxs  # Load indices :   Holds the object indices for when the system is loaded back in

        # Inherent traits
        self.num = index
        self.element = element      # Symbol       :   Element of the atom
        self.chain = chain          # Chain        :   Molecule chain the atom is a part of
        self.mol_class = mol_class  # Mol Class    :   Class of molecule that the atom is a part of
        self.res_seq = res_seq      # Sequence     :   Sequence of the residue that the atom is a part of
        self.name = name            # Name         :   Name retrieved from pdb file
        self.occupancy = ocp        # Occupancy    :   Occupancy of the atom
        self.t_fact = t_fact        # Temp Factor  :   Temperature factor for the atom
        self.seg_id = seg_id        # Segment ID   :   Segment identifier for the atom
        self.charge = charge        # Charge       :   Charge of the atom

        self.sort_()

        # Calculated traits
    def sort_(self):
        """
        Puts the atom in the correct spot in the system
        :return:
        """
        # If no system exists, there is no place to be sorted to
        if self.sys is None:
            return
        # Find the molecule

    def calc_vol(self):
        # Create the volume variable
        vol = 0
        # Go through each surface on the atom
        for surf in self.surfs:
            self.sa += surf.sa
            if surf.tris is None or len(surf.tris) == 0:
                if surf.file is not None:
                    try:
                        surf.read_file(surf.file)
                    except FileNotFoundError:
                        surf.build()
                    except PermissionError:
                        surf.build()
                else:
                    surf.build()

            for tri in surf.tris:
                if tri is None:
                    print(surf.tris)
                p0, p1, p2, p3 = self.loc, surf.points[tri[0]], surf.points[tri[1]], surf.points[tri[2]]
                vol += calc_tetra_vol(p0, p1, p2, p3)
        # Return the volume
        self.vol = vol
        return vol



def get_radius(radius, system, return_symbol=False):
    """
        Finds the radius of the atom from the symbol or vice versa

    :param radius: Either the elemental symbol for the atom or it's radius
    :param system: System to reference radii from
    :param return_symbol: Boolean for whether to return the symbol or not
    :return: The radius of the atom from the symbol or vice versa
    """
    radii = system.radii
    # If indicated we return the symbol of atom that the radius indicates
    if return_symbol:
        # Set the atom type to nothing
        atom_type = ""
        # Check to see if the radius is in the system
        if radius in radii[1]:
            return radii[0][radii[1].index(radius)]
        else:
            # Get the closest atom to it
            min_diff = np.inf
            # Go through the radii in the system looking for the smallest difference
            for i in range(len(radii[1])):
                if radii[1][i] is not None and radii[1][i] - radius < min_diff:
                    atom_type = radii[0][i]
    # If we have the type and just want the radius, keep scanning until we find the radius
    else:
        radius = radius.strip()
        return radii[1][radii[0].index(radius.lower())]
    # If nothing is found to be exact return the closest atom type
    return atom_type
