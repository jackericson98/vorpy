"""Residue categories shared by loading and visualization."""
from vorpy.src.boundary import ION_RESIDUES, SOLVENT_RESIDUES, WATER_RESIDUES

def is_solvent(atom):
    return str(atom.residue_name).strip().upper() in SOLVENT_RESIDUES
