from System.sys_funcs.calcs import *


class Residue:
    def __init__(self, atoms=None, name=None, sys=None, mol=None, sequence=None, seg_id=None):
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
