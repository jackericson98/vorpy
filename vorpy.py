# Outside Imports
import os
import multiprocessing as mp
from Visualize.vorpy_gui import Vorpy, ErrorBox, LoadingBox
from System.system import System


if __name__ == '__main__':
    os.chdir("..")
    # Launch the gui
    myVorpy = Vorpy()
    # Create the System
    if myVorpy.sys_file_address is not None:
        mySys = System(myVorpy.sys_file_address, box_size=myVorpy.sys_box_x_flt.get())
        mySys.name = myVorpy.sys_file_name.get()
        if myVorpy.vert_file_address is not None:
            mySys.add_verts(myVorpy.vert_file_address)
    elif len(myVorpy.sys_atom_list) > 0:
        mySys = System(myVorpy.sys_atom_list, box_size=myVorpy.sys_box_x_flt.get())
        mySys.name = "User_Data"
    else:
        ErrorBox("Please select a file or enter atoms")
        mySys = None

    # If no outer directory was specified use the current one
    if myVorpy.sys_output_directory is None:
        myVorpy.sys_output_directory = os.getcwd()
    # Catch for existing directories
    i = 0
    while True:
        # Try creating the directory
        try:
            i_str = str(i)
            if i == 0:
                i_str = ""
            os.mkdir(os.getcwd() + "/" + mySys.name + i_str)
            break
        except FileExistsError:
            i += 1
    # Set the output directory for the system
    mySys.output_directory = os.getcwd() + "/" + mySys.name + i_str
    mySys.net.find_verts()
    mySys.export_verts()
    mySys.net.connect()

    surfs = mySys.net.surfs
    for i in range(len(surfs)):
        # Calculate and print the running percentage for mesh calculations
        percentage = int((i + 1) / len(surfs) * 100)
        print("\rBuilding Surfaces: ",
              '#' * (percentage // 10) + ' ' * (10 - (percentage // 10)), percentage, "%", end='')
        surfs[i].build(simps=True, min_dist=myVorpy.sys_res_flt.get())
    print("\rBuilding Surfaces:  ########## 100 %")
    print("\rSurfaces Built")

    # Analyze the network
    if myVorpy.output_analysis.get() or myVorpy.output_all.get():
        mySys.analyze()

    # Export the system
    mySys.export(directory=mySys.output_directory, export_all=myVorpy.output_all.get(),
                       export_sys=myVorpy.output_sys.get(), export_mols=myVorpy.output_mols.get(),
                       export_atoms=myVorpy.output_atoms.get(), export_analysis=myVorpy.output_analysis.get(),
                       export_surfs=myVorpy.output_surfs.get())
    os.chdir(myVorpy.vorpy_directory)
