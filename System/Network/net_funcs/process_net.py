from System.sys_funcs.calcs import *

# The point of this file is to pre-process all the information for the system





def calc_sas(net):
    # Go through every surface and calculate it's surface area
    for surf in net.surfs:
        surf.sa = calc_sa(surf)


def calc_vols(net):
    # Go through each atom in the system and calculate its volume
    for atom in net.atoms:
        calc_vol(atom)


def calc_curve(net):
    # Go through each surface and calculate the curvature
    for surf in net.surfs:
        pass
