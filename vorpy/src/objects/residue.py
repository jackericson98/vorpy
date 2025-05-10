

class Residue:
    """
    A class representing a molecular residue in a system.
    
    A residue is a fundamental building block of biomolecules (proteins, nucleic acids, etc.)
    that contains a collection of atoms. This class manages residue-specific information
    and provides methods for atom management.

    Attributes:
        atoms (list): List of Atom objects belonging to this residue
        name (str): Name of the residue (e.g., 'ALA' for alanine)
        sys: Reference to the parent system object
        mol: Reference to the parent molecule object
        seq (int): Sequence number/position of the residue in the chain
        id (str): Segment identifier for the residue
        chain: Reference to the parent chain object
        print_name (str): Formatted name for display/output purposes
    """
    def __init__(self, atoms=None, name=None, sys=None, mol=None, sequence=None, seg_id=None, chain=None):
        """
        Initialize a new Residue object.

        Args:
            atoms (list, optional): List of Atom objects to initialize the residue with
            name (str, optional): Name of the residue (e.g., 'ALA', 'GLY')
            sys (optional): Reference to the parent system object
            mol (optional): Reference to the parent molecule object
            sequence (int, optional): Sequence number/position in the chain
            seg_id (str, optional): Segment identifier
            chain (optional): Reference to the parent chain object
        """
        self.atoms = atoms if atoms is not None else []
        self.name = name
        self.sys = sys
        self.mol = mol
        self.seq = sequence
        self.id = seg_id
        self.chain = chain
        self.print_name = None

    def add_atom(self, atom):
        """
        Add an atom to the residue's atom list.

        Args:
            atom: The Atom object to add to the residue
        """
        self.atoms.append(atom)

nucleic_acids = {'DT', 'DA', 'DG', 'DC', 'DU', 'U', 'G', 'A', 'T', 'C', 'GDP', 'OMC'}
amino_acids = {'ALA', 'ARB', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
               'THR', 'TRP', 'TYR', 'VAL', 'GLY', 'ARG'}