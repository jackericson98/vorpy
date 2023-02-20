from Visualize.commands.interpret import *
from Visualize.commands.group import group


def export(sys, usr_npt, my_group=None):
    """
    Takes in input strings and exports them based on their option choices
    :param sys:
    :param usr_npt:
    :return:
    """
    if my_group is None:
        my_group = group(sys, usr_npt)

    print("choose one to export: 1. Shell, 2. Surfaces, 3. Layers, 4. Atoms, 5. Filled Body, 6. Info File, 7. Vertices, 8. All")
    while True:
        # Export the group exports
        xpt_npt = input("export >>>   ")
        xpt_npt = xpt_npt.strip()
        # Check for a quit
        if xpt_npt.lower() in quits:
            return
        # Check for help request
        elif xpt_npt.lower() in helps:
            help_()
            continue
        # Export the shell:
        elif xpt_npt.lower() in ['1', '1.', 'shell', 'sh']:
            my_group.exports(shell=True)
            print("\r{} shell exported to {}".format(my_group.name, my_group.dir))
        # Export the Surfaces
        elif xpt_npt.lower() in ['2', '2.', 'surfs', 'surfaces']:
            my_group.exports(surfaces=True)
            print("\r{} surfaces exported to {}".format(my_group.name, my_group.dir + "/surfaces"))
        # Export the layers
        elif xpt_npt.lower() in ['3', '3.', 'layers', 'lyrs', 'l']:
            my_group.exports(layers=True)
            print("\r{} layers exported to {}".format(my_group.name, my_group.dir + "/layers"))
        # Export the atoms
        elif xpt_npt.lower() in ['4', '4.', 'atoms', 'a', 'atms']:
            my_group.exports(atoms=True)
            print("\r{} atoms exported to {}".format(my_group.name, my_group.dir))
        # Export the filled body
        elif xpt_npt.lower() in ['5', '5.', 'filled body', 'fb', 'f', 'filled_body', 'f_b', 'fld_bdy']:
            my_group.exports(fill=True)
            print("\r{} filled body exported to {}".format(my_group.name, my_group.dir))
        # Export the Info file
        elif xpt_npt.lower() in ['6', '6.', 'info', 'info_file', 'info file']:
            my_group.exports(info=True)
            print("\r{} info file exported to {}".format(my_group.name, my_group.dir))
        # Export vertices
        elif xpt_npt.lower() in ['7', '7.', 'verts', 'vertices']:
            my_group.exports(verts=True)
            print("\r {} exported to {}".format(my_group.name, my_group.dir))
            # Export all
        elif xpt_npt.lower() in ['8', '8.', 'eight', 'ocho', '8 ', 'all', 'a']:
            my_group.exports(all_=True)
            print("\r{} shell, surfaces, layers, atoms, filled body and info file exported to {}".format(group.name,                                                                                                   group.dir))
        else:
            invalid_input(xpt_npt)
