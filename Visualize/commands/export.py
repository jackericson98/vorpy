from Visualize.commands.interpret import *
from Visualize.commands.group import group


def export_surfs(sys, my_group=None):
    # Main prompt
    print("{} Surface Exports: 1. Shell 2. Filled Body 3. Separate Surfaces 4. All".format(my_group.name))
    while True:
        # Input
        surf_npt = input("surfacee (1-4) >>>   ")
        surf_npt.strip()
        # Export the shell
        alls = ['4', '4.', 'four', 'all', 'a']
        if surf_npt.lower() in ['1', '1.', 'one', 'shell', 'sh'] + alls:
            my_group.exports(shell=True)
            print("\r{} shell surfaces exported to {}".format(my_group.name, my_group.dir))
        if surf_npt.lower() in ['2', '2.', 'two', 'filled', 'fb'] + alls:
            my_group.exports(fill=True)
            print("\r{} filled surfaces exported to {}".format(my_group.name, my_group.dir))
        if surf_npt.lower() in ['3', '3.', 'three', 'surfs', 'surfaces'] + alls:
            my_group.exports(surfaces=True)
            print("\r{} surfaces exported to {}".format(my_group.name, my_group.dir))
        if surf_npt.lower() in quits:
            break



def export_atoms(sys, my_group=None):
    # Main prompt
    print("{} Atoms Exports 1. All Atoms 2. Surrounding Atoms 3. Exterior Atoms 4. All".format(my_group.name))
    while True:
        # Input
        surf_npt = input("atoms (1-4) >>>   ")
        surf_npt.strip()
        # Export the shell
        alls = ['4', '4.', 'four', 'all', 'a']
        if surf_npt.lower() in ['1', '1.', 'one', 'atoms'] + alls:
            my_group.exports(atoms=True)
            print("\r{} atoms exported to {}".format(my_group.name, my_group.dir + "/surfaces"))
        if surf_npt.lower() in ['2', '2.', 'two', 'sur', 'sa'] + alls:
            my_group.exports(surr_atoms=True)
            print("\r{} surrounding atoms exported to {}".format(my_group.name, my_group.dir + "/surfaces"))
        if surf_npt.lower() in ['3', '3.', 'three', 'ext', 'ea'] + alls:
            my_group.exports(ext_atoms=True)
            print("\r{} exterior atoms exported to {}".format(my_group.name, my_group.dir + "/surfaces"))
        if surf_npt.lower() in quits:
            break


def export_info(sys, my_group=None):
    # Main prompt
    print("{} Information 1. Show Info 2. Info Export".format(my_group.name))
    while True:
        # Input Prompt
        nfo_npt = input("info (1-2) >>>   ")
        # Export the shell
        alls = ['4', '4.', 'four', 'all', 'a']
        if nfo_npt.lower() in ['1', '1.', 'one', 'show', 'sh'] + alls:
            my_group.get_info()
            print("Volume:", my_group.vol, "Surface Area:", my_group.sa)
        if nfo_npt.lower() in ['2', '2.', 'two'] + alls:
            my_group.exports(info=True)
            print("\r{} information file exported to {}".format(my_group.name, my_group.dir))
        if nfo_npt.lower() in quits:
            break


def export_verts(sys, my_group=None):
    # Main prompt
    print("{} Network Exports 1. Vertices 2. Edges 3. Shell Vertices 4. Shell Edges".format(my_group.name))
    while True:
        # Input
        verts_npt = input("network (1-4) >>>   ")
        if verts_npt.lower() in ['1', '1.', 'one', 'verts', 'vertices']:
            my_group.exports(verts=True)
            print("\r{} vertices exported to {}".format(my_group.name, my_group.dir))
        if verts_npt.lower() in ['2', '2.', 'two', "edges", "e"]:
            my_group.exports(edges=True)
            print("\r{} edges exported to {}".format(my_group.name, my_group.dir))
        if verts_npt.lower() in ['3', '3.', 'three', 'sv', 'shell_verts']:
            my_group.exports(shell_verts=True)
            print("\r{} shell vertices exported to {}".format(my_group.name, my_group.dir))
        if verts_npt.lower() in ['4', '4.', 'four', "se", "shell_edges"]:
            my_group.exports(shell_edges=True)
            print("\r{} shell edges exported to {}".format(my_group.name, my_group.dir))
        if verts_npt.lower() in quits:
            break


def export(sys, usr_npt, my_group=None):
    """
    Takes in input strings and exports them based on their option choices
    :param my_group:
    :param sys:
    :param usr_npt:
    :return:
    """
    if my_group is None:
        my_group = group(sys, usr_npt)
    if sys.dir is None:
        sys.set_output_directory()
    # Main prompt
    print("Exports 1. Surfaces and interfaces 2. Atoms 3. Information 4. Vertices and Edges 5. All")
    while True:
        # Export the group exports
        xpt_npt = input("export type (1-5) >>>   ")
        xpt_npt = xpt_npt.strip()
        # Check for a quit
        if xpt_npt.lower() in quits + ns:
            return
        # Check for help request
        elif xpt_npt.lower() in helps:
            help_()
            continue
        alls = ['5', '5.', 'five', 'all', 'a']
        # Export the shell:
        if xpt_npt.lower() in ['1', '1.', 'surfaces'] + alls:
            # Call the Surfaces function
            export_surfs(sys=sys, my_group=my_group)
        # Export the Surfaces
        if xpt_npt.lower() in ['2', '2.', 'atoms', 'a'] + alls:
            export_atoms(sys=sys, my_group=my_group)
        # Export the layers
        if xpt_npt.lower() in ['3', '3.', 'info', 'information', 'i'] + alls:
            export_info(sys=sys, my_group=my_group)
        # Export the atoms
        if xpt_npt.lower() in ['4', '4.', 'vertices', 'verts', 'v'] + alls:
            export_verts(sys=sys, my_group=my_group)
