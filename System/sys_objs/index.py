from System.sys_funcs.input import read_ndx


class Index:
    """
    Index object used to hold gromacs style index files and their lists
    """
    def __init__(self, sys, atom_range=None, name=None, file=None, file_ndx=None):
        self.sys = sys
        self.name = name
        self.atoms = None
        self.atom_rng = atom_range
        self.file = file
        self.file_ndx = file_ndx

        if self.file is not None:
            self.load()

    # Load method. Used to load the details of an index file
    def load(self, file=None, file_ndx=None):
        if file is not None:
            self.file = file
        if file_ndx is not None:
            self.file_ndx = file_ndx
        read_ndx(sys=self.sys, file=self.file)
