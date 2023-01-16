from System.sys_funcs.calcs import *


class Residue:
    def __init__(self, atoms=None, sys=None, mol=None, sequence=None, seg_id=None):
        self.atoms = atoms
        self.sys = sys
        self.mol = mol
        self.seq = sequence
        self.id = seg_id
