from Analysis.anal_calcs import *
from Analysis.find_simps import find_simps


# Analyze system function. Finds the surfaces and volumes of the system
def analyze(sys):
    # Go through each surface in the system and find the simplices and the surface area
    for surf in sys.net.surf:
        # Get the surfaces simplices
        surf.simps = find_simps(surf)
        # Get the surface area of the surface
        surf.sa = calc_sa(surf)

    # Go through each atom in the system
    for atom in sys.atoms:
        atom.cell_vol = calc_vol(atom)
