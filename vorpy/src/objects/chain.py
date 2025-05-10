

class Chain:
    """
    A class representing a molecular chain in a system.
    
    Attributes:
        name (str): The name of the chain
        atoms (list): List of atoms belonging to this chain
        residues (list): List of residues in this chain
        vol (float): Volume of the chain (default: 0)
        sa (float): Surface area of the chain (default: 0)
    """
    def __init__(self, sys=None, atoms=None, residues=None, name=None):
        # Initialize chain attributes
        self.name = name  # Name of the chain
        self.atoms = atoms  # List of atoms in the chain
        self.residues = residues  # List of residues in the chain
        self.vol = 0  # Volume of the chain
        self.sa = 0  # Surface area of the chain

    def add_atom(self, atom):
        """
        Add an atom to the chain's atom list.
        
        Args:
            atom: The atom object to add to the chain
        """
        self.atoms.append(atom)


class Sol(Chain):
    """
    A specialized Chain class representing a solvent molecule (default: water).
    Inherits all attributes and methods from the Chain class.
    
    Attributes:
        name (str): Name of the solvent (default: "H2O")
        atoms (list): List of atoms in the solvent molecule
        residues (list): List of residues in the solvent molecule
        vol (float): Volume of the solvent (default: 0)
        sa (float): Surface area of the solvent (default: 0)
    """
    def __init__(self, sys=None, atoms=None, residues=None, name="H2O"):
        # Initialize the parent Chain class
        super().__init__()
        # Initialize solvent-specific attributes
        self.atoms = atoms  # List of atoms in the solvent
        self.residues = residues  # List of residues in the solvent
        self.name = name  # Name of the solvent (defaults to "H2O")
        self.vol = 0  # Volume of the solvent
        self.sa = 0  # Surface area of the solvent
