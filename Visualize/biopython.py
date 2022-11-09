from Bio.PDB import *
import nglview as nv
import ipywidgets

pdb_parser = PDBParser()
structure = pdb_parser.get_structure("PHA-L", "../Data/test_data/Na5.pdb")
view = nv.show_biopython(structure)
