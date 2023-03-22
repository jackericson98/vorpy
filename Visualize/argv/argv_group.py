from Visualize.cmnd.group import group
from Visualize.cmnd.commands import *
from System.Group.group import Group


def argv_group(my_sys, usr_npt, add_more=False, bff=None):
    # Check if a bff was specified
    if bff is not None and len(bff) >= 1:
        bff = group(sys=my_sys, usr_npt=bff[0])

    if len(usr_npt) == 0 or usr_npt[0][0] == 'ns':
        my_sys.groups.append(Group(sys=my_sys, chains=my_sys.chains, name="no_SOL"))
    else:
        for grouping in usr_npt:
            my_sys.groups.append(group(sys=my_sys, usr_npt=grouping, bff=bff))


def group_argv(my_sys, usr_npt):
    # Create a group variable
    my_group = None
    # Go through the user inputs loading files
    while usr_npt:
        # Pop the file descriptor
        descriptor = usr_npt.pop(0)
        # Check to see that it is a descriptor
        if descriptor.lower() not in my_objects:
            return
        elif descriptor.lower() == 'ns':
            my_group = Group(my_sys, chains=my_sys.chains, name=my_sys.name + "_no_SOL")
            continue
        elif descriptor.lower in full_objs:
            return Group(my_sys, atoms=my_sys.atoms, name=my_sys.name + "_full")
        # Load the file
        my_group = group(my_sys, [descriptor, usr_npt.pop(0)], my_group)
        # If the next value is && go again
        if len(usr_npt) > 0 and usr_npt[0] == '+':
            usr_npt.pop(0)
    # Return the group object
    return my_group
