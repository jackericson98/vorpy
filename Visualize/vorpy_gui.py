import os
import tkinter as tk
import tkinter.ttk as ttk
from tkinter import filedialog
from System.system import System
from System.group import Group


# Loading gui class. Holds the settings for the load/settings gui
class Vorpy:

    def __init__(self, width=750, height=650):

        self.main = tk.Tk()
        self.main.geometry(str(width) + "x" + str(height))
        self.main.title('vorpy')
        self.running = False
        self.exporting = False
        self.index_file = None
        self.index_list = []
        self.output_directory = os.getcwd()
        self.net_file = None
        self.sys = System(atoms=[], mols=[], residues=[])
        self.g1 = None
        self.g2 = None

        self.mol_list, self.res_list, self.atom_list = [], [], []

        self.sys_name = tk.StringVar(self.main, "No System Selected")
        self.net_name = tk.StringVar(self.main, "No Network Selected")
        self.sys_file = None
        self.net_file_address = None
        self.frame_files = []
        self.user_atoms = []

        """_______________________________________________Load Frame_________________________________________________"""

        # Load subframe attributes
        self.load_frame = tk.Frame(self.main, name="load")
        self.load_frame.pack()

        # Load frame attributes
        self.sys_pros = tk.StringVar(self.main, "   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n")
        self.net_sets = tk.StringVar(self.main, "   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n")

        # Build subframe attributes
        self.sys_res_flt = tk.DoubleVar(self.main, 0.1)
        self.sys_alpha_value = tk.DoubleVar(self.main, 5)
        self.sys_box_x_flt = tk.DoubleVar(self.main, 1.4)

        ############################################## Load Header subframe ############################################

        # Set up the frame, place it and add the separator
        tk.Label(self.load_frame, text="VorPy", font=("Times New Roman bold", 40))\
            .grid(row=0, column=0, columnspan=3, padx=10, pady=10)
        ttk.Separator(self.load_frame).grid(row=1, columnspan=3, sticky='ew')

        ############################################## Load system subframe ############################################

        # Set up the frame and place it
        load_system_frame = tk.Frame(self.load_frame, width=375)
        load_system_frame.grid(row=2, column=0, padx=10, pady=10)
        # Write and place the frame's title
        tk.Label(load_system_frame, text="Load System", font=("Times New Roman Bold", 10))\
            .grid(row=0, column=0, columnspan=3, sticky='nw')
        # Create the system's name with a default letting the user know nothing has been selected, then place it
        tk.Label(load_system_frame, textvariable=self.sys_name, font=("Times New Roman bold", 15))\
            .grid(row=1, column=0, columnspan=3)
        # Create the "Add Atoms", "Load System" and "Load Frames" buttons and place them
        tk.Button(load_system_frame, text="Add Atoms ", command=self.add_atoms_button).grid(row=2, column=0, sticky='w')
        tk.Button(load_system_frame, text="Load System ", command=self.load_sys_button).grid(row=2, column=1)
        tk.Button(load_system_frame, text="Load Frames ", command=self.load_frames_button)\
            .grid(row=2, column=2, sticky='e')
        # Add the system information header and labels
        tk.Label(load_system_frame, text="System Information").grid(row=3, column=1, sticky='n')
        tk.Label(load_system_frame, text="Atoms:\nMolecules:\nSolution:\n   ~   \n   ~   \n   ~   \n")\
            .grid(row=4, column=0, sticky='w')
        tk.Label(load_system_frame, text=":\n:\n:\n:\n:\n:").grid(row=4, column=1)
        # Add the system information to match the labels
        tk.Label(load_system_frame, textvariable=self.sys_pros).grid(row=4, column=2)

        # Add a seperator for the load system and load network frames
        ttk.Separator(self.load_frame, orient=tk.VERTICAL).grid(row=1, column=1, rowspan=3, sticky='ns')

        ############################################# Load network subframe ############################################

        # Set up and place the load network frame
        load_network_frame = tk.Frame(self.load_frame, width=375)
        load_network_frame.grid(row=2, column=2, padx=10, pady=10)
        # Add the header for the load network frame
        tk.Label(load_network_frame, text="Load Network", font=("Times New Roman bold", 10))\
            .grid(row=0, column=0, columnspan=3, sticky='nw')
        # Set the network's name at the top of the frame
        tk.Label(load_network_frame, textvariable=self.net_name, font=("Times New Roman bold", 15))\
            .grid(row=1, column=0, columnspan=3)
        # Add the buttons for loading the network and analyzing it
        tk.Button(load_network_frame, text="Load Network", command=self.load_network).grid(row=2, column=1)
        tk.Button(load_network_frame, text="Analyze", command=self.load_analyze_button).grid(row=2, column=2)
        # Report the loaded network's information
        tk.Label(load_network_frame, text="Network Information").grid(row=3, column=1, sticky='n')
        tk.Label(load_network_frame, text="Box Size:\nResolution:\nMax Vertex:\n    ~    \n    ~    \n    ~    \n")\
            .grid(row=4, column=0)
        tk.Label(load_network_frame, text=":\n:\n:\n:\n:\n:").grid(row=4, column=1)
        tk.Label(load_network_frame, textvariable=self.net_sets).grid(row=4, column=2)

        # Separate the system and network frame from the build frame
        ttk.Separator(self.load_frame).grid(row=3, columnspan=3, sticky='ew')

        ################################################ Build Subframe ################################################

        # Instantiate the frame
        bld_sbfrm = tk.Frame(self.load_frame)
        bld_sbfrm.grid(row=4, column=0, columnspan=3, pady=10, padx=10)

        # Set the title for the subframe
        tk.Label(bld_sbfrm, text="Build Network", font=("Times New Roman bold", 15))\
            .grid(row=0, column=0, sticky='n', columnspan=6)

        # System resolution value
        tk.Label(bld_sbfrm, text="\nResolution: ", font=("Times New Roman bold", 12))\
            .grid(row=1, column=0,  columnspan=3)
        # Slider for setting the network's resolution and then a label for the units
        tk.Scale(bld_sbfrm, from_=0.005, to=0.5, orient=tk.HORIZONTAL, variable=self.sys_res_flt,
                 resolution=0.005).grid(row=2, column=0, columnspan=3, sticky='ew')
        tk.Label(bld_sbfrm, text=u'\u212B', font=('Times New Roman', 15)).grid(row=2, column=3, sticky='s')

        # Maximum vertex radius value for the system
        tk.Label(bld_sbfrm, text="\nMax Vertex:", font=("Times New Roman bold", 12)).grid(row=3, column=0, columnspan=3)
        # Slider for setting the network's resolution and then a label for the units
        tk.Scale(bld_sbfrm, from_=2, to=10, orient=tk.HORIZONTAL, variable=self.sys_alpha_value,
                 resolution=0.05).grid(row=4, column=0, columnspan=3, sticky='ew')
        tk.Label(bld_sbfrm, text=u'\u212B', font=('Times New Roman', 15)).grid(row=4, column=3, sticky='s')

        # Box size
        tk.Label(bld_sbfrm, text="Box Size: ", font=("Times New Roman bold", 12)).grid(row=5, column=0, columnspan=3)
        tk.Scale(bld_sbfrm, from_=1.05, to=3, orient=tk.HORIZONTAL, variable=self.sys_box_x_flt,
                 resolution=.05, length=400).grid(row=6, column=0, columnspan=3, sticky='ew')
        tk.Label(bld_sbfrm, text="x", font=('Times New Roman', 20)).grid(row=6, column=3, sticky='s')

        # Parallelize check
        self.parallelize = tk.BooleanVar(self.main)
        tk.Checkbutton(bld_sbfrm, text="Parallelize ", variable=self.parallelize, onvalue=True,
                       offvalue=False).grid(row=2, column=4, sticky='w', padx=10)

        # Find solution vertices check
        self.sol_verts = tk.BooleanVar(self.main, True)
        tk.Checkbutton(bld_sbfrm, text="Find SOL Verts", variable=self.sol_verts, onvalue=True,
                       offvalue=False).grid(row=3, column=4, sticky='w', padx=10)

        # Curved faces check
        self.curved_faces = tk.BooleanVar(self.main, True)
        tk.Checkbutton(bld_sbfrm, text="Curved Faces", variable=self.curved_faces, onvalue=True, offvalue=False) \
            .grid(row=4, column=4, sticky='w', padx=10)

        # Parallelize check
        self.flat_faces = tk.BooleanVar(self.main, False)
        tk.Checkbutton(bld_sbfrm, text="Flat Faces", variable=self.flat_faces, onvalue=True, offvalue=False) \
            .grid(row=2, column=5, sticky='w')

        # Use loaded vertices check
        self.use_loaded_verts = tk.BooleanVar(self.main, False)
        tk.Checkbutton(bld_sbfrm, text="Use Loaded Vertices", variable=self.use_loaded_verts, onvalue=True,
                       offvalue=False).grid(row=3, column=5, sticky='w')

        # Build network button
        tk.Button(bld_sbfrm, text="Build Network", font=("Times New Roman bold", 20),
                  command=self.build_network_button).grid(row=5, rowspan=2, column=4, columnspan=2)

        """########################################### Build Frame ##################################################"""

        self.build_frame = tk.Frame(self.main, name="build")
        # self.build_frame.pack()

        # Header
        tk.Label(self.build_frame, text="VorPy", font=("Times New Roman bold", 40)).grid(row=0, column=0,
                                                                                         columnspan=3, padx=10, pady=10)
        ttk.Separator(self.build_frame).grid(row=1, columnspan=3, sticky='ew')

        # System Name
        tk.Label(self.build_frame, textvariable=self.net_name, font=("Times New Roman bold", 20)).grid(row=2, padx=10,
                                                                                                       pady=10)

        # Network information
        build_net_info_frame = tk.Frame(self.build_frame)
        build_net_info_frame.grid(row=3, padx=10, pady=10)
        tk.Label(build_net_info_frame, text="Network Information", font=("Times New Roman bold", 15)).grid(row=0,
                                                                                                           column=1)
        tk.Label(build_net_info_frame, text="Vertices").grid(row=1, column=1)
        tk.Label(build_net_info_frame, text="Surfaces").grid(row=1, column=2)
        tk.Label(build_net_info_frame, text="Number:\n\nmy Time:\n\nCPU Time:\n\nSettings:\n").grid(row=2, column=0,
                                                                                                    sticky='w')

        # Progress
        build_net_progress_frame = tk.Frame(self.build_frame)
        build_net_progress_frame.grid(row=4, padx=10, pady=20)
        tk.Label(build_net_progress_frame, text="Progress", font=("Times New Roman bold", 20)).grid(row=0, columnspan=3)
        self.current_process = tk.IntVar(self.main, 0)
        self.canvas_width, self.canvas_height = 400, 40
        self.progress_canvas = tk.Canvas(build_net_progress_frame, bg='white', width=self.canvas_width,
                                         height=self.canvas_height)
        self.progress_canvas.grid(row=1, column=0, columnspan=3, pady=10)
        self.update_progress_canvas(0)

        # Current process progress
        self.processes = ["Loading System", "Finding Vertices", "Connecting Network", "Building Surfaces",
                          "Analyzing System"]
        self.current_process_str = tk.StringVar(self.main, self.processes[self.current_process.get()])
        self.percentage = tk.StringVar(self.main, "0%")
        tk.Label(build_net_progress_frame, textvariable=self.current_process_str).grid(row=2, column=0, sticky='w')
        tk.Label(build_net_progress_frame, textvariable=self.percentage).grid(row=2, column=2, sticky='e')
        self.loading_bar = ttk.Progressbar(build_net_progress_frame, length=400, mode='determinate')
        self.loading_bar.grid(row=3, columnspan=3)

        # Cancel button
        tk.Button(self.build_frame, text='Quit').grid(row=5)

        """######################################### Analysis Frame #################################################"""

        # Set up the analysis frame
        self.analysis_frame = tk.Frame(self.main, name="analysis")

        # Header
        tk.Label(self.analysis_frame, text="VorPy", font=("Times New Roman bold", 40)).grid(row=0, columnspan=3,
                                                                                            padx=10, pady=10,)
        # Separator
        ttk.Separator(self.analysis_frame).grid(row=1, columnspan=3, sticky='ew')

        # Network name
        tk.Label(self.analysis_frame, textvariable=self.net_name, font=("Times New Roman bold", 20))\
            .grid(row=2, columnspan=3, padx=10, pady=10)

        # Separator
        ttk.Separator(self.analysis_frame).grid(row=3, columnspan=3, sticky='ew')

        # Network information List
        anal_net_frame = tk.Frame(self.analysis_frame)
        anal_net_frame.grid(row=4, column=0, padx=10, pady=10)

        tk.Label(anal_net_frame, text="Network Information:", font=("Times New Roman bold", 20))\
            .grid(row=0, column=0, columnspan=2, sticky='nw')

        ttk.Separator(anal_net_frame, orient=tk.VERTICAL).grid(row=1, rowspan=3, column=1, sticky='ns')
        ttk.Separator(anal_net_frame, orient=tk.HORIZONTAL).grid(row=2, columnspan=3, column=0, sticky='ew')

        # Create the system information sub frame
        anal_sys_info_subframe = tk.Frame(anal_net_frame)
        anal_sys_info_subframe.grid(row=1, column=0, padx=5, pady=5)

        tk.Label(anal_sys_info_subframe, text="System Information").grid(row=3, columnspan=3, sticky='n')
        tk.Label(anal_sys_info_subframe, text="Atoms:\nMolecules:\nResidues:\n    ~    \n    ~    \n    ~    \n") \
            .grid(row=4, column=0)
        tk.Label(anal_sys_info_subframe, text=":\n:\n:\n:\n:\n:").grid(row=4, column=1)
        tk.Label(anal_sys_info_subframe, textvariable=self.sys_pros).grid(row=4, column=2)

        # Create the system information sub frame
        anal_cpu_info_subframe = tk.Frame(anal_net_frame)
        anal_cpu_info_subframe.grid(row=1, column=2, padx=5, pady=5)

        tk.Label(anal_cpu_info_subframe, text="CPU Information").grid(row=3, columnspan=3, sticky='n')
        tk.Label(anal_cpu_info_subframe, text="My Time:\nCPU Time:\nComputations:\n    ~    \n    ~    \n    ~    \n") \
            .grid(row=4, column=0)
        tk.Label(anal_cpu_info_subframe, text=":\n:\n:\n:\n:\n:").grid(row=4, column=1)
        self.cpu_info = tk.StringVar(self.main, "   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n")
        tk.Label(anal_cpu_info_subframe, textvariable=self.cpu_info).grid(row=4, column=2)

        # Create the system information sub frame
        anal_settings_subframe = tk.Frame(anal_net_frame)
        anal_settings_subframe.grid(row=3, column=0, padx=5, pady=5)

        tk.Label(anal_settings_subframe, text="Settings Information").grid(row=3, columnspan=3, sticky='n')
        tk.Label(anal_settings_subframe,
                 text="Box Size:\nResolution:\nMax Vertex:\n    ~    \n    ~    \n    ~    \n").grid(row=4, column=0)
        tk.Label(anal_settings_subframe, text=":\n:\n:\n:\n:\n:").grid(row=4, column=1)
        tk.Label(anal_settings_subframe, textvariable=self.net_sets).grid(row=4, column=2)

        # Create the system information sub frame
        anal_outputs_info_subframe = tk.Frame(anal_net_frame)
        anal_outputs_info_subframe.grid(row=3, column=2, padx=5, pady=5)

        tk.Label(anal_outputs_info_subframe, text="Output Information").grid(row=3, columnspan=3, sticky='n')
        tk.Label(anal_outputs_info_subframe, text="Vertices:\nSurfaces:\nEdges:\nFull Cells\n   ~   \n   ~   ") \
            .grid(row=4, column=0)
        tk.Label(anal_outputs_info_subframe, text=":\n:\n:\n:\n:\n:").grid(row=4, column=1)
        self.output_info = tk.StringVar(self.main, "   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n")
        tk.Label(anal_outputs_info_subframe, textvariable=self.output_info).grid(row=4, column=2)

        # Seperator
        ttk.Separator(self.analysis_frame, orient=tk.VERTICAL).grid(row=3, column=1, rowspan=3, sticky='ns')

        # Export information List
        group_info_frame = tk.Frame(self.analysis_frame)
        group_info_frame.grid(row=4, column=2, sticky='ns', padx=10, pady=10)

        tk.Label(group_info_frame, text="Group Information:", font=("Times New Roman bold", 18))\
            .grid(row=0, column=0, sticky='nw', columnspan=5)
        tk.Label(group_info_frame, text="Group 1").grid(row=1, column=2)
        ttk.Separator(group_info_frame, orient=tk.VERTICAL).grid(row=1, rowspan=3, column=1, sticky='ns')
        ttk.Separator(group_info_frame, orient=tk.VERTICAL).grid(row=1, rowspan=3, column=3, sticky='ns')

        tk.Label(group_info_frame, text="Group 2").grid(row=1, column=4)
        ttk.Separator(group_info_frame).grid(row=2, columnspan=5, sticky='ew')
        tk.Label(group_info_frame, text="Surface Area:\n\nVolume:\n\nBoundary atoms:\n\nOuter Atoms:\n")\
            .grid(row=3, column=0, sticky='w')
        self.g1_info = tk.StringVar(self.main, "   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n   ~   ")
        tk.Label(group_info_frame, textvariable=self.g1_info).grid(row=3, column=2)
        self.g2_info = tk.StringVar(self.main, "   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n   ~   ")
        tk.Label(group_info_frame, textvariable=self.g2_info).grid(row=3, column=4)
        ttk.Separator(group_info_frame).grid(row=4, column=0, columnspan=5, sticky='ew')
        tk.Label(group_info_frame, text="Interface:").grid(row=6, column=1)
        tk.StringVar()
        tk.Label(group_info_frame, text="   ~   \n   ~   \n   ~   ")

        self.interface_info = tk.StringVar(self.main, )

        ttk.Separator(self.analysis_frame).grid(row=5, columnspan=5, sticky='ew')

        # Main frame Network object gathering
        export_obj_frame = tk.Frame(self.analysis_frame)
        export_obj_frame.grid(row=6, column=0, columnspan=3, padx=10, pady=10)

        # Group sub frame
        self.using_group1 = tk.BooleanVar(self.main, True)
        self.using_group2 = tk.BooleanVar(self.main, not self.using_group1.get())
        anal_check_subframe = tk.Frame(export_obj_frame)
        anal_check_subframe.grid(column=0, row=0, padx=10, pady=10)

        tk.Label(anal_check_subframe, text="Groups:", font=("Times New Roman bold", 20)).grid(row=0)
        self.group1_selections = tk.StringVar(self.main, "")
        self.group2_selections = tk.StringVar(self.main, "")
        tk.Checkbutton(anal_check_subframe, text="Group 1", font=("Times New Roman", 15), variable=self.using_group1,
                       command=self.flip_g2).grid(row=1, column=0)
        tk.Label(anal_check_subframe, text="Selections").grid(row=2, column=0)
        tk.Label(anal_check_subframe, textvariable=self.group1_selections).grid(row=3, column=0)
        tk.Checkbutton(anal_check_subframe, text="Group 2", font=("Times New Roman", 15), variable=self.using_group2,
                       command=self.flip_g1).grid(row=4, column=0)
        tk.Label(anal_check_subframe, text="Selections").grid(row=5, column=0)
        tk.Label(anal_check_subframe, textvariable=self.group2_selections).grid(row=6, column=0)

        # Choose index

        # Seperator
        ttk.Separator(export_obj_frame, orient=tk.VERTICAL).grid(column=1, row=0, sticky='ns', rowspan=3)
        # Choose index sub frame
        self.choose_index_subframe = tk.Frame(export_obj_frame)
        self.choose_index_subframe.grid(column=2, row=0, sticky='ns', padx=10, pady=10)

        # Set the label header for the "choose index" sub frame
        tk.Label(self.choose_index_subframe, text="Get Atoms:", font=("Times New Roman bold", 20))\
            .grid(row=0, columnspan=3, sticky='nw')
        tk.Label(self.choose_index_subframe, text="Choose Index:", font=("Times New Roman bold", 15))\
            .grid(row=1, column=0)
        self.current_ndx = tk.StringVar(self.main, "[None]")
        tk.OptionMenu(self.choose_index_subframe, variable=self.current_ndx, value=self.index_list)\
            .grid(row=1, column=1)
        tk.Button(self.choose_index_subframe, text="Browse", command=self.load_index).grid(row=1, column=2)

        # Create the index
        # Atoms
        tk.Label(self.choose_index_subframe, text="Create Index", font=("Times New Roman bold", 12))\
            .grid(row=3, column=0, columnspan=3)
        self.current_atom_selection = tk.StringVar(self.main, "")
        tk.Label(self.choose_index_subframe, text="Atom").grid(row=4, column=0)
        self.atom_options = tk.OptionMenu(self.choose_index_subframe, self.current_atom_selection, "", *self.sys.atoms)
        self.atom_options.grid(row=5, column=0)

        # Molecules
        self.current_mol_selection = tk.StringVar(self.main)
        tk.Label(self.choose_index_subframe, text="Molecule").grid(row=4, column=1)
        self.mol_options = tk.OptionMenu(self.choose_index_subframe, self.current_mol_selection, "", *self.sys.mols)
        self.mol_options.grid(row=5, column=1)

        # Residues
        self.current_res_selection = tk.StringVar(self.main, "")
        tk.Label(self.choose_index_subframe, text="Residue").grid(row=4, column=2)
        self.res_options = tk.OptionMenu(self.choose_index_subframe, self.current_res_selection, "", *self.sys.residues)
        self.res_options.grid(row=5, column=2)

        # Selection Button
        tk.Button(self.choose_index_subframe, text="Reset Group", command=self.reset_group).grid(row=6, column=0)
        tk.Button(self.choose_index_subframe, text="Undo Selection", command=self.undo_selection).grid(row=6, column=1)
        tk.Button(self.choose_index_subframe, text="Add Selection", command=self.add_selection).grid(row=6, column=2)

        # Seperator
        ttk.Separator(export_obj_frame, orient=tk.VERTICAL).grid(column=3, row=0, sticky='ns', rowspan=3)

        # Export object sub frame
        export_obj_subframe = tk.Frame(export_obj_frame)
        export_obj_subframe.grid(column=4, row=0, padx=10, pady=10, sticky='ns')

        tk.Label(export_obj_subframe, text="Export:", font=("Times New Roman bold", 20))\
            .grid(column=0, row=0, columnspan=2, sticky='nw')
        tk.Label(export_obj_subframe, text="Change Output Directory:").grid(row=1, column=0)
        tk.Button(export_obj_subframe, text="Browse", command=self.change_output_directory).grid(row=1, column=1)
        self.output_dir_str = tk.StringVar(self.main, os.getcwd()[:12] + "..." + "/Data/User_data")
        tk.Label(export_obj_subframe, textvariable=self.output_dir_str).grid(row=2, column=0, columnspan=2)
        self.export_info = tk.BooleanVar(self.main, False)
        tk.Checkbutton(export_obj_subframe, text="Export Info", variable=self.export_info)\
            .grid(row=3, column=0, columnspan=2)

        tk.Button(export_obj_subframe, text="Export", font=("Times New Roman bold", 20),
                  command=self.export_selections).grid(row=4, column=0, columnspan=2, rowspan=2)

        self.main.mainloop()

    """############################################# Functions  #####################################################"""

    # Load system button function. Calls the file browser and sets the system
    def load_sys_button(self):
        # Reset the system
        self.sys = System()
        # Re-connect the gui and the system
        self.sys.gui = self
        # File grabber pop up
        file_path = filedialog.askopenfilename()
        # Set the file path
        if file_path:
            self.sys_file = file_path
        # Get the file name
        filename = ""
        i = -1
        while self.sys_file[i] != "/":
            filename = filename + self.sys_file[i]
            i -= 1
        # Set the file name
        self.sys.name = filename[::-1][:-4]
        self.sys_name.set(self.sys.name)
        # Load the system
        self.sys.load_sys(file_path)
        # Set the system information
        # We want to get the number of atoms, the number of molecules, etc
        myStr = str(len(self.sys.atoms)) + '\n' + str(len(self.sys.mols)) + '\n' + str(len(self.sys.residues)) + \
                "\n   ~   \n   ~   \n   ~   "
        # Set the variables
        self.sys_pros.set(myStr)
        # Set the networks atoms
        self.sys.net.atoms = self.sys.atoms
        # Set the molecule names
        for mol in self.sys.mols:
            self.mol_list.append(mol[0].chain)
        for res in self.sys.residues:
            self.res_list.append(res[0].res + " " + res[0].res_seq)
        for atom in self.sys.atoms:

            self.atom_list.append(str(self.sys.atoms.index(atom)) + " " + atom.element)
        self.mol_options.destroy()
        self.mol_options = tk.OptionMenu(self.choose_index_subframe, self.current_mol_selection, *self.sys.mol_names)
        self.mol_options.grid(row=5, column=1)

        self.res_options.destroy()
        self.res_options = tk.OptionMenu(self.choose_index_subframe, self.current_res_selection, *self.sys.res_names)
        self.res_options.grid(row=5, column=2)

        self.atom_options.destroy()
        self.atom_options = tk.OptionMenu(self.choose_index_subframe, self.current_atom_selection, *self.sys.atom_names)
        self.atom_options.grid(row=5, column=0)

    # Load frames function.
    def load_frames_button(self):
        # File grabber pop up
        file_path = filedialog.askopenfilename()
        self.frame_files.append(file_path)

    # Load network button function. Pulls up the file browser and lets the user select their vorpy saved system
    def load_network(self):
        # File grabber pop up
        file_path = filedialog.askopenfilename()
        self.net_file = file_path
        # Check to see if there is a system file
        if len(self.sys.atoms) < 1:
            return
        else:
            self.sys.load_net(self.net_file)
        self.set_analyze_info()
        self.sys.set_output_directory()

    # Set analyze frame info. When the analysis frame is pulled up this is called and sets the info
    def set_analyze_info(self):
        # We want to get the number of atoms, the number of molecules, etc
        myStr = str(self.sys.net.box_size) + '\n' + str(self.sys.net.min_dist) + '\n' + str(self.sys.net.beta_val) + \
                "\n   ~   \n   ~   \n   ~   "
        # Set the variables
        self.net_sets.set(myStr)
        # Set the network name
        self.net_name.set(self.sys_name.get() + "_network")
        # Set the output info and the cpu info
        myStr = str(self.sys.net.my_time) + '\n' + str(self.cpu_info) + '\n' + "   ~   \n   ~   \n   ~   \n   ~   "
        self.cpu_info.set(myStr)
        myStr = str(len(self.sys.net.verts)) + '\n' + str(len(self.sys.net.edges)) + '\n' + \
                str(len(self.sys.net.surfs)) + "\n   ~   \n   ~   \n   ~   "

        self.output_info.set(myStr)

    # Load analysis button method. Moves the screen to the analysis screen once a network has been loaded
    def load_analyze_button(self):

        if self.net_file is not None:
            self.load_frame.pack_forget()
            self.main.geometry("800x900")
            self.analysis_frame.pack()

    # Load index file method. Allows the user to load
    def load_index(self):
        # File grabber pop up
        file_path = filedialog.askopenfilename()
        self.index_file = file_path

    # Add atoms button method. Adds user atoms to the system for analysis
    def add_atoms_button(self):
        pass

    # Build network button method. Cements the settings and destroys the gui
    def build_network_button(self):
        if self.sys is None:
            return

        self.sys.build_network()
        self.load_frame.destroy()
        self.main.geometry("800x900")
        self.analysis_frame.pack()
        self.set_analyze_info()

    # Change output directory method. Updates the location of the output directory
    def change_output_directory(self):
        # File grabber pop up
        file_path = filedialog.askdirectory()
        # Create the System
        if file_path:
            self.output_directory = file_path
            self.output_dir_str.set(file_path[:12] + ' ... ' + file_path[-12:])

    # Update progressbar function. Updates the color of the progressbar depending on the stage that vorpy is in
    def update_progress_canvas(self, task_ndx=None):

        # Set up the canvas variable
        canvas = self.progress_canvas
        # Set up the dimensions
        w, h, = int(self.canvas_width/4), int(self.canvas_height/2)

        # Check to see we aren't moving to a specific task
        if task_ndx is None:
            # Move to the next task and set the string for the task
            self.current_process.set(self.current_process.get() + 1)
            self.current_process_str.set(self.processes[self.current_process.get()])
            task_ndx = self.current_process.get()

        # Determine the colors of the arrows
        colors = ['green' for _ in range(task_ndx)]
        colors += ['yellow']
        colors += ['red' for _ in range(4 - task_ndx - 1)]

        # Create the polygons
        canvas.create_polygon(0, 0, w-h/2, 0, w+h/2, h, w-h/2, 2*h, 0, 2*h, fill=colors[0], outline='black')
        canvas.create_polygon(w-h/2, 0, 2*w-h/2, 0, 2*w+h/2, h, 2*w-h/2, 2*h, w-h/2, 2*h, w+h/2, h, fill=colors[1],
                              outline='black')
        canvas.create_polygon(2*w-h/2, 0, 3*w-h/2, 0, 3*w+h/2, h, 3*w-h/2, 2*h, 2*w-h/2, 2*h, 2*w+h/2, h,
                              fill=colors[2], outline='black')
        canvas.create_polygon(3*w-h/2, 0, 4*w, 0, 4*w, 2*h, 3*w-h/2, 2*h, 3*w+h/2, h, fill=colors[3], outline='black')
        # Add the text on the polygons
        canvas.create_text(0.5 * w, h, text="Finding\nVertices", font=("Times New Roman bold", 10))
        canvas.create_text(1.5 * w, h, text="Connecting\nNetwork", font=("Times New Roman bold", 10))
        canvas.create_text(2.5 * w, h, text="Building\nSurfaces", font=("Times New Roman bold", 10))
        canvas.create_text(3.5 * w, h, text="Analyzing\nNetwork", font=("Times New Roman bold", 10))

    def flip_g1(self):
        self.using_group1.set(not self.using_group1.get())

    def flip_g2(self):
        self.using_group2.set(not self.using_group2.get())

    def reset_group(self):
        if self.using_group1.get():
            self.g1 = None
            self.group2_selections.set("")
            if self.g2 is not None:
                self.g2 = Group(self.sys.net, self.g2.atoms)
                self.g2.get_info()
        if self.using_group2.get():
            self.g2 = None
            self.group2_selections.set("")
            if self.g1 is not None:
                self.g1 = Group(self.sys.net, self.g1.atoms)
                self.g1.get_info()

    def undo_selection(self):
        if self.using_group1.get() and len(self.g1.atoms) >= 1:
            self.g1.atoms = self.g1.atoms[:-min(self.g1.prev_sele, len(self.g1.atoms))]
            self.group1_selections.set(self.group1_selections.get()[:-len(self.g1.prev_str)])

        elif self.using_group2.get() and len(self.g2.atoms) >= 1:
            self.g2.atoms = self.g2.atoms[:-self.g2.prev_sele]
            self.group2_selections.set(self.group2_selections.get()[:-min(self.g2.prev_sele, len(self.g2.atoms))])

    def add_selection(self):
        atoms = []
        # Get the atoms from the three lists
        if self.current_mol_selection.get() != '':
            atoms += self.sys.mols[self.sys.mol_names.index(self.current_mol_selection.get()) - 1]
        if self.current_res_selection.get() != '':
            atoms += self.sys.residues[self.sys.res_names.index(self.current_res_selection.get()) - 1]
        if self.current_atom_selection.get() != '':
            atoms += [self.sys.atoms[self.sys.atom_names.index(self.current_atom_selection.get()) - 1]]
        if len(atoms) == 0:
            return

        # Sort out the atoms that are already in the list
        for atom in atoms:
            if self.using_group1 and self.g1 is not None and atom in self.g1.atoms:
                atoms.pop(atoms.index(atom))
            elif self.using_group2 and self.g2 is not None and atom in self.g2.atoms:
                atoms.pop(atoms.index(atom))

        current_string = self.current_mol_selection.get() + self.current_res_selection.get() + \
                         self.current_atom_selection.get()
        # Add the atoms to the current group selection
        if self.using_group1.get():
            if self.g1 is None:
                self.g1 = Group(self.sys.net, atoms)
                self.g1.name = "Group1_" + "".join([str(self.sys.atoms.index(_)) for _ in
                                                    self.g1.atoms[:min(12, len(self.g1.atoms) - 1)]])
            else:
                self.g1.atoms += atoms
            self.g1.prev_str = current_string
            self.g1.prev_sele = len(atoms)
            self.group1_selections.set(self.group1_selections.get() + "\n" + current_string)
            self.g1.get_info()
            self.g1_info.set(str(round(self.g1.body_sa, 3)) + "\n\n" + str(round(self.g1.body_vol, 3)) + "\n\n" +
                             str(len(self.g1.outer_body_atoms)) + "\n\n" + str(len(self.g1.surr_body_atoms)))

            if self.g2 is not None:
                self.g2.bff, self.g1.bff = self.g1, self.g2
        elif self.using_group2:
            if self.g2 is None:
                self.g2 = Group(self.sys.net, atoms)
                self.g2.name = "Group2_" + "".join([str(self.sys.atoms.index(_))
                                                    for _ in self.g2.atoms[:min(12, len(self.g2.atoms) - 1)]])
            else:
                self.g2.atoms += atoms
            self.g2.prev_str = current_string
            self.g2.prev_sele = len(atoms)
            self.group2_selections.set(self.group2_selections.get() + "\n" + current_string)
            self.g2.get_info()
            self.g2_info.set(str(round(self.g2.body_sa, 3)) + "\n" + str(round(self.g2.body_vol, 3)) + "\n" +
                             str(len(self.g2.outer_body_atoms)) + "\n" + str(len(self.g2.surr_body_atoms)))
            if self.g2 is not None:
                self.g2.bff, self.g1.bff = self.g1, self.g2

    def export_selections(self):
        if self.g1 is not None:
            self.sys.export_selection(self.g1, group2=self.g2, info=self.export_info.get())
