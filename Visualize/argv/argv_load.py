from Visualize.cmnd.load import load


def argv_load(my_sys, usr_npt):
    if len(usr_npt) > 0:
        my_npts = []
        for npt in usr_npt:
            my_npts += ['', npt[0]]
        load(my_sys, my_npts)


def argv_load_atoms(my_sys, usr_npt):
    if len(usr_npt) > 0:
        load(my_sys, usr_npt)
