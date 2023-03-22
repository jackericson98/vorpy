

class Residue:
    def __init__(self, atoms=None, name=None, sys=None, mol=None, sequence=None, seg_id=None, chain=None):
        """
        Residue Object for holding specific residue information
        :param atoms:
        :param name:
        :param sys:
        :param mol:
        :param sequence:
        :param seg_id:
        """
        self.atoms = atoms
        self.name = name
        self.sys = sys
        self.mol = mol
        self.seq = sequence
        self.id = seg_id
        self.chain = chain
        self.print_name = None
        h_count, o_count = 0, 0
        for atom in atoms:
            if atom.element.lower() == 'h':
                h_count += 1
            elif atom.element.lower() == 'o':
                o_count += 1
        if h_count == 2 and o_count == 1:
            self.print_name = 'water'
        else:
            self.print_name = "I Dont Know Gosh"

    def add_atom(self, atom):
        self.atoms.append(atom)
