

class Chain:
    def __init__(self, atoms=None, residues=None, name=None):
        self.name = name
        self.atoms = atoms
        self.residues = residues

    def add_atom(self, atom):
        self.atoms.append(atom)


class Sol(Chain):
    def __init__(self, atoms=None, residues=None, name="H2O"):
        super().__init__()
        self.atoms = atoms
        self.residues = residues
        self.name = name
