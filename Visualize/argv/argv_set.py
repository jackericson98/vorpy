from Visualize.cmnd.set import sett, my_settings


def argv_sett(my_sys, usr_npt):
    print(usr_npt)
    # Go through the user inputs loading files
    while usr_npt:
        my_npt = usr_npt.pop(0)
        # Pop the file descriptor
        descriptor = my_npt[0]
        # Check to see that it is a descriptor
        if descriptor.lower() not in my_settings or len(usr_npt) == 0:
            return
        # Load the file
        sett(my_sys, [descriptor, my_npt[1]], vorpy2_set=False)
