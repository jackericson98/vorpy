from System.sys_funcs.calcs import *


class Molecule:
    def __init__(self, atoms=None, sys=None, name=None, residues=None):
        self.atoms = atoms
        self.resids = residues
        self.sys = sys
        self.name = name
