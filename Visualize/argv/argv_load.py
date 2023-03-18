from Visualize.cmnd.load import load


def argv_load(my_sys, usr_npt):
    if len(usr_npt) > 0:
        load(my_sys, usr_npt)
