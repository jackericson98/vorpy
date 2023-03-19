from Visualize.cmnd.set import sett
from Visualize.cmnd.commands import my_settings


def argv_sett(my_sys, usr_npt):
    # Go through the user inputs loading files
    while usr_npt:
        my_npt = usr_npt.pop(0)
        # Pop the file descriptor
        descriptor = my_npt[0]
        # Check to see that it is a descriptor
        if descriptor.lower() not in my_settings or len(my_npt) == 0:
            return
        # Load the file
        sett(my_sys, my_npt, vorpy2_set=False)
