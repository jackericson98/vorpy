"""Residue categories shared by loading and visualization."""
WATER_RESIDUES = {"HOH", "WAT", "SOL", "TIP3", "TIP3P", "SPC", "SPCE"}
ION_RESIDUES = {
    "LI",
    "NA",
    "K",
    "RB",
    "CS",
    "MG",
    "CA",
    "SR",
    "BA",
    "ZN",
    "CD",
    "FE",
    "FE2",
    "FE3",
    "MN",
    "CU",
    "CU1",
    "CU2",
    "CO",
    "NI",
    "AL",
    "F",
    "CL",
    "BR",
    "I",
    "IOD",
    "SO4",
    "PO4",
    "NH4",
}


SOLVENT_RESIDUES = WATER_RESIDUES | ION_RESIDUES

def is_solvent(atom):
    return atom.residue_name.strip().upper() in SOLVENT_RESIDUES
