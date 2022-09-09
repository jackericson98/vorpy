# Outside Imports
import os
import multiprocessing as mp
from Visualize.vorpy_gui import Vorpy, ErrorBox, LoadingBox


if __name__ == '__main__':
    # Launch the gui
    myVorpy = Vorpy()

    myLoading = LoadingBox()
    # Create the System
    if myVorpy.sys is None:
        ErrorBox("Please select a system")
    # Create a new directory of one has not been indicated
    if myVorpy.sys.output_directory is None:
        # If the system doesn't have a name
        if myVorpy.sys.name == '':
            myVorpy.sys.name = "User_Data"
        i = 0
        while True:
            try:
                i_str = str(i)
                if i == 0:
                    i_str = ""
                os.mkdir(os.getcwd() + "/" + myVorpy.sys.name + i_str)
                break
            except FileExistsError:
                i += 1
        myVorpy.sys.output_directory = os.getcwd() + "/" + myVorpy.sys.name + i_str


    # Define the process that we want to split up into different cores of the cpu.
    # We have find vertices (Hardest), Connect Network (Similarly as hard), Build surfaces
    # Build the network
    # myVorpy.sys.build_network(get_verts=not myVorpy.verts_loaded,
    #                           export_verts=(myVorpy.output_verts.get() or myVorpy.output_all.get()),
    #                           directory=myVorpy.output_directory, box_size=float(myVorpy.box_size.get()),
    #                           min_dist=float(myVorpy.resolution.get()))
    myVorpy.sys.net.find_verts()
    myVorpy.sys.export_verts()
    myVorpy.sys.net.connect()

    surfs = myVorpy.sys.net.surfs
    for i in range(len(surfs)):
        # Calculate and print the running percentage for mesh calculations
        percentage = int((i + 1) / len(surfs) * 100)
        print("\rBuilding Surfaces: ",
              '#' * (percentage // 10) + ' ' * (10 - (percentage // 10)), percentage, "%", end='')
        surfs[i].build(simps=True, min_dist=float(myVorpy.sys_res_str.get()))
    print("\rBuilding Surfaces:  ########## 100 %")
    print("\rSurfaces Built")

    # Analyze the network
    if myVorpy.output_analysis.get() or myVorpy.output_all.get():
        myVorpy.sys.analyze()

    # Export the system
    myVorpy.sys.export(directory=myVorpy.sys.output_directory, export_all=myVorpy.output_all.get(),
                       export_sys=myVorpy.output_sys.get(), export_mols=myVorpy.output_mols.get(),
                       export_atoms=myVorpy.output_atoms.get(), export_analysis=myVorpy.output_analysis.get(),
                       export_surfs=myVorpy.output_surfs.get())
    os.chdir(myVorpy.vorpy_directory)
