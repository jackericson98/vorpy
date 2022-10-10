import os
import tkinter as tk
import tkinter.ttk as ttk
from tkinter import filedialog
from System.system import System, Network


# Loading gui class. Holds the settings for the load/settings gui
class Vorpy:

    def __init__(self, sys=None, width=750, height=650):

        self.main = tk.Tk()
        self.main.geometry(str(width) + "x" + str(height))
        self.main.title('vorpy')
        self.running = False
        self.exporting = False
        self.atom_list = ["Molecullllle"]
        self.mol_lis = ["Molecullllle"]
        self.res_list = ["Molecullllle"]
        self.index_list = []
        self.output_directory = os.getcwd()
        self.net_file = None
        self.sys = None

        #############################################  Load Frame Attributes ###########################################

        self.load_frame = tk.Frame(self.main, name="load")
        self.load_frame.pack()

        # Set up the loading attributes
        self.sys_file = None
        self.frame_files = []
        self.user_atoms = []


        # Set up the load network
        self.net_file_address = None

        # Header frame
        tk.Label(self.load_frame, text="VorPy", font=("Times New Roman bold", 40)).grid(row=0, column=0, columnspan=3, padx=10, pady=10)
        ttk.Separator(self.load_frame).grid(row=1, columnspan=3, sticky='ew')

        # File stuff
        load_system_frame = tk.Frame(self.load_frame, width=375)
        load_system_frame.grid(row=2, column=0, padx=10, pady=10)
        tk.Label(load_system_frame, text="Load System", font=("Times New Roman Bold", 10)).grid(row=0, column=0, columnspan=3, sticky='nw')
        self.sys_name = tk.StringVar(self.main, "No System Selected")
        tk.Label(load_system_frame, textvariable=self.sys_name, font=("Times New Roman bold", 15)).grid(row=1, column=0, columnspan=3)
        tk.Button(load_system_frame, text="Add Atoms ", command=self.add_atoms_button).grid(row=2, column=0, sticky='w')
        tk.Button(load_system_frame, text="Load System ", command=self.load_sys_button).grid(row=2, column=1)
        tk.Button(load_system_frame, text="Load Frames ", command=self.load_frames_button).grid(row=2, column=2, sticky='e')
        tk.Label(load_system_frame, text="System Information").grid(row=3, column=1, sticky='n')
        tk.Label(load_system_frame, text="Atoms:\nMolecules:\nSolution:\n   ~   \n   ~   \n   ~   \n")\
            .grid(row=4, column=0, sticky='w')
        tk.Label(load_system_frame, text=":\n:\n:\n:\n:\n:").grid(row=4, column=1)
        self.sys_pros = tk.StringVar(self.main, "   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n")
        tk.Label(load_system_frame, textvariable=self.sys_pros).grid(row=4, column=2)

        ttk.Separator(self.load_frame, orient=tk.VERTICAL).grid(row=1, column=1, rowspan=3, sticky='ns')

        # Load network stuff
        load_network_frame = tk.Frame(self.load_frame, width=375)
        load_network_frame.grid(row=2, column=2, padx=10, pady=10)
        tk.Label(load_network_frame, text="Load Network", font=("Times New Roman bold", 10)).grid(row=0, column=0, columnspan=3, sticky='nw')
        self.net_name = tk.StringVar(self.main, "No Network Selected")
        tk.Label(load_network_frame, textvariable=self.net_name, font=("Times New Roman bold", 15)).grid(row=1, column=0, columnspan=3)
        tk.Button(load_network_frame, text="Load Network", command=self.load_network).grid(row=2, column=1)
        tk.Label(load_network_frame, text="Network Information").grid(row=3, column=1, sticky='n')
        tk.Label(load_network_frame, text="Box Size:\nResolution:\nMax Vertex:\n    ~    \n    ~    \n    ~    \n")\
            .grid(row=4, column=0)
        tk.Label(load_network_frame, text=":\n:\n:\n:\n:\n:").grid(row=4, column=1)
        self.net_sets = tk.StringVar(self.main, "   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n")
        tk.Label(load_network_frame, textvariable=self.net_sets).grid(row=4, column=2)

        ttk.Separator(self.load_frame).grid(row=3, columnspan=3, sticky='ew')

        # Settings frame
        load_settings_frame = tk.Frame(self.load_frame)
        load_settings_frame.grid(row=4, column=0, columnspan=3, pady=10, padx=10)

        tk.Label(load_settings_frame, text="Build Network", font=("Times New Roman bold", 15)).grid(row=0, column=0, sticky='n', columnspan=6)

        # System resolution value
        self.sys_res_flt = tk.DoubleVar(self.main, 0.1)
        tk.Label(load_settings_frame, text="\nResolution: ", font=("Times New Roman bold", 12)).grid(row=1, column=0,  columnspan=3)
        tk.Scale(load_settings_frame, from_=0.005, to=0.5, orient=tk.HORIZONTAL, variable=self.sys_res_flt, resolution=0.005) \
            .grid(row=2, column=0, columnspan=3, sticky='ew')
        tk.Label(load_settings_frame, text=u'\u212B', font=('Times New Roman', 15)).grid(row=2, column=3, sticky='s')

        # Maximum vertex radius value for the system
        self.sys_alpha_value = tk.DoubleVar(self.main, 5)
        tk.Label(load_settings_frame, text="\nMax Vertex: ", font=("Times New Roman bold", 12)).grid(row=3, column=0,  columnspan=3)
        tk.Scale(load_settings_frame, from_=2, to=10, orient=tk.HORIZONTAL, variable=self.sys_alpha_value, resolution=0.05) \
            .grid(row=4, column=0, columnspan=3, sticky='ew')
        tk.Label(load_settings_frame, text=u'\u212B', font=('Times New Roman', 15)).grid(row=4, column=3, sticky='s')

        # System box size multiplier
        self.sys_box_x_flt = tk.DoubleVar(self.main, 1.4)
        tk.Label(load_settings_frame, text="Box Size: ", font=("Times New Roman bold", 12)).grid(row=5, column=0,  columnspan=3)
        tk.Scale(load_settings_frame, from_=1.05, to=3, orient=tk.HORIZONTAL, variable=self.sys_box_x_flt, resolution=.05, length=400) \
            .grid(row=6, column=0, columnspan=3, sticky='ew')
        tk.Label(load_settings_frame, text="x", font=('Times New Roman', 20)).grid(row=6, column=3, sticky='s')


        # Parallelize check
        self.parallelize = tk.BooleanVar(self.main)
        tk.Checkbutton(load_settings_frame, text="Parallelize ", variable=self.parallelize, onvalue=True, offvalue=False) \
            .grid(row=2, column=4, sticky='w', padx=10)

        # Find solution vertices check
        self.measure_sol = tk.BooleanVar(self.main, True)
        tk.Checkbutton(load_settings_frame, text="Find SOL Verts", variable=self.measure_sol, onvalue=True, offvalue=False) \
            .grid(row=3, column=4, sticky='w', padx=10)

        # Curved faces check
        self.curved_faces = tk.BooleanVar(self.main, True)
        tk.Checkbutton(load_settings_frame, text="Curved Faces", variable=self.curved_faces, onvalue=True, offvalue=False) \
            .grid(row=4, column=4, sticky='w', padx=10)

        # Parallelize check
        self.flat_faces = tk.BooleanVar(self.main, False)
        tk.Checkbutton(load_settings_frame, text="Flat Faces", variable=self.flat_faces, onvalue=True, offvalue=False) \
            .grid(row=2, column=5, sticky='w')

        # Use loaded vertices check
        self.use_loaded_verts = tk.BooleanVar(self.main, False)
        tk.Checkbutton(load_settings_frame, text="Use Loaded Vertices", variable=self.use_loaded_verts, onvalue=True, offvalue=False) \
            .grid(row=3, column=5, sticky='w')


        tk.Button(load_settings_frame, text="Build Network", font=("Times New Roman bold", 20), command=self.build_network_button) \
            .grid(row=5, rowspan=2, column=4, columnspan=2)



        ########################################### Build Frame ########################################################

        self.build_frame = tk.Frame(self.main, name="build")

        # Header
        tk.Label(self.build_frame, text="VorPy", font=("Times New Roman bold", 40)).grid(row=0, column=0, columnspan=3, padx=10, pady=10)
        ttk.Separator(self.build_frame).grid(row=1, columnspan=3, sticky='ew')

        # System Name
        tk.Label(self.build_frame, textvariable=self.net_name, font=("Times New Roman bold", 20)).grid(row=2, padx=10, pady=10)

        # Network information
        build_net_info_frame = tk.Frame(self.build_frame)
        build_net_info_frame.grid(row=3, padx=10, pady=10)
        tk.Label(build_net_info_frame, text="Network Information", font=("Times New Roman bold", 15)).grid(row=0, column=1)
        tk.Label(build_net_info_frame, text="Vertices").grid(row=1, column=1)
        tk.Label(build_net_info_frame, text="Surfaces").grid(row=1, column=2)
        tk.Label(build_net_info_frame, text="Number:\n\nmy Time:\n\nCPU Time:\n\nSettings:\n").grid(row=2, column=0, sticky='w')

        # Progress
        build_net_progress_frame = tk.Frame(self.build_frame)
        build_net_progress_frame.grid(row=4, padx=10, pady=20)
        tk.Label(build_net_progress_frame, text="Progress", font=("Times New Roman bold", 20)).grid(row=0, columnspan=3)
        self.current_process = tk.IntVar(self.main, 0)
        self.canvas_width, self.canvas_height = 400, 40
        self.progress_canvas = tk.Canvas(build_net_progress_frame, bg='white', width=self.canvas_width, height=self.canvas_height)
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



        ######################################### Analysis Frame #######################################################

        # Set up the analysis frame
        self.analysis_frame = tk.Frame(self.main, name="analysis")

        # Header
        tk.Label(self.analysis_frame, text="VorPy", font=("Times New Roman bold", 40)).grid(row=0, columnspan=3, padx=10, pady=10,)

        # Seperator
        ttk.Separator(self.analysis_frame).grid(row=1, columnspan=3, sticky='ew')

        # Network name
        tk.Label(self.analysis_frame, textvariable=self.net_name, font=("Times New Roman bold", 20))\
            .grid(row=2, columnspan=3, padx=10, pady=10)

        # Seperator
        ttk.Separator(self.analysis_frame).grid(row=3, columnspan=3, sticky='ew')

        # Network information List
        anal_net_frame = tk.Frame(self.analysis_frame)
        anal_net_frame.grid(row=4, column=0, padx=10, pady=10)

        tk.Label(anal_net_frame, text="Network Information:", font=("Times New Roman bold", 20)).grid(row=0, column=0, columnspan=2, sticky='nw')

        ttk.Separator(anal_net_frame, orient=tk.VERTICAL).grid(row=1, rowspan=3, column=1, sticky='ns')
        ttk.Separator(anal_net_frame, orient=tk.HORIZONTAL).grid(row=2, columnspan=3, column=0, sticky='ew')

        # Create the system information sub frame
        anal_sys_info_subframe = tk.Frame(anal_net_frame)
        anal_sys_info_subframe.grid(row=1, column=0, padx=5, pady=5)

        tk.Label(anal_sys_info_subframe, text="System Information").grid(row=3, columnspan=3, sticky='n')
        tk.Label(anal_sys_info_subframe, text="Atoms:\nMolecules:\nResidues:\n    ~    \n    ~    \n    ~    \n") \
            .grid(row=4, column=0)
        tk.Label(anal_sys_info_subframe, text=":\n:\n:\n:\n:\n:").grid(row=4, column=1)
        self.sys_info = tk.StringVar(self.main, "   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n")
        tk.Label(anal_sys_info_subframe, textvariable=self.sys_info).grid(row=4, column=2)

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
        tk.Label(anal_settings_subframe, text="Box Size:\nResolution:\nMax Vertex:\n    ~    \n    ~    \n    ~    \n") \
            .grid(row=4, column=0)
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

        tk.Label(group_info_frame, text="Export Information:", font=("Times New Roman bold", 20)).grid(row=0, column=0, sticky='nw', columnspan=2)
        tk.Label(group_info_frame, text="Group 1").grid(row=1, column=0)
        ttk.Separator(group_info_frame, orient=tk.VERTICAL).grid(row=1, rowspan=3, column=1, sticky='ns')
        tk.Label(group_info_frame, text="Group 2").grid(row=1, column=2)
        ttk.Separator(group_info_frame).grid(row=2, columnspan=3, sticky='ew')
        tk.Label(group_info_frame, text="Number:\n\nmy Time:\n\nCPU Time:\n\nSettings:\n").grid(row=3, column=0,
                                                                                                   sticky='w')
        ttk.Separator(self.analysis_frame).grid(row=5, columnspan=3, sticky='ew')



        # Main frame Network object gathering
        export_obj_frame = tk.Frame(self.analysis_frame)
        export_obj_frame.grid(row=6, column=0, columnspan=3, padx=10, pady=10)



        # Group sub frame
        self.using_group1 = tk.BooleanVar(self.main, True)
        self.using_group2 = tk.BooleanVar(self.main, not self.using_group1.get())
        anal_check_subframe = tk.Frame(export_obj_frame)
        anal_check_subframe.grid(column=0, row=0, padx=10, pady=10)


        tk.Label(anal_check_subframe, text="Groups:", font=("Times New Roman bold", 20)).grid(row=0)
        self.group1_selections = tk.StringVar(self.main, "No sub-groups selected")
        self.group2_selections = tk.StringVar(self.main, "No sub-groups selected")
        tk.Checkbutton(anal_check_subframe, text="Group 1", font=("Times New Roman", 15), variable=self.using_group1, command=self.flip_g2).grid(row=1, column=0)
        tk.Label(anal_check_subframe, text="Selections").grid(row=2, column=0)
        tk.Label(anal_check_subframe, textvariable=self.group1_selections).grid(row=3, column=0)
        tk.Checkbutton(anal_check_subframe, text="Group 2", font=("Times New Roman", 15), variable=self.using_group2, command=self.flip_g1).grid(row=4, column=0)
        tk.Label(anal_check_subframe, text="Selections").grid(row=5, column=0)
        tk.Label(anal_check_subframe, textvariable=self.group2_selections).grid(row=3, column=0)


        # Choose index

        # Seperator
        ttk.Separator(export_obj_frame, orient=tk.VERTICAL).grid(column=1, row=0, sticky='ns', rowspan=3)
        # Choose index sub frame
        choose_index_subframe = tk.Frame(export_obj_frame)
        choose_index_subframe.grid(column=2, row=0, sticky='ns', padx=10, pady=10)

        # Set the label header for the "choose index" sub frame
        tk.Label(choose_index_subframe, text="Get Atoms:", font=("Times New Roman bold", 20)).grid(row=0, columnspan=3, sticky='nw')
        tk.Label(choose_index_subframe, text="Choose Index:", font=("Times New Roman bold", 15)).grid(row=1, column=0)
        self.current_ndx = tk.StringVar(self.main, "[None]")
        tk.OptionMenu(choose_index_subframe, variable=self.current_ndx, value=self.index_list).grid(row=1, column=1)
        tk.Button(choose_index_subframe, text="Browse", command=self.load_index).grid(row=1, column=2)

        # Create the index
        # Atoms
        tk.Label(choose_index_subframe, text="Create Index", font=("Times New Roman bold", 15)).grid(row=3, column=0)
        self.current_atom_selection = tk.StringVar(self.main, "None")
        tk.Label(choose_index_subframe, text="Atom").grid(row=4, column=0)
        tk.OptionMenu(choose_index_subframe, variable=self.current_atom_selection, value=self.atom_list).grid(row=5, column=0)

        # Molecules
        self.current_mol_selection = tk.StringVar(self.main, "None")
        tk.Label(choose_index_subframe, text="Molecule").grid(row=4, column=1)
        tk.OptionMenu(choose_index_subframe, variable=self.current_mol_selection, value=self.atom_list).grid(row=5, column=1)

        # Residues
        self.current_res_selection = tk.StringVar(self.main, "None")
        tk.Label(choose_index_subframe, text="Residue").grid(row=4, column=2)
        tk.OptionMenu(choose_index_subframe, variable=self.current_res_selection, value=self.atom_list).grid(row=5, column=2)


        # Selection Button
        tk.Button(choose_index_subframe, text="Reset Group", command=self.reset_group).grid(row=6, column=0)
        tk.Button(choose_index_subframe, text="Undo Selection", command=self.undo_selection).grid(row=6, column=1)
        tk.Button(choose_index_subframe, text="Add Selection", command=self.add_selection).grid(row=6, column=2)

        # Seperator
        ttk.Separator(export_obj_frame, orient=tk.VERTICAL).grid(column=3, row=0, sticky='ns', rowspan=3)

        # Export object sub frame
        export_obj_subframe = tk.Frame(export_obj_frame)
        export_obj_subframe.grid(column=4, row=0, padx=10, pady=10, sticky='ns')


        tk.Label(export_obj_subframe, text="Export:", font=("Times New Roman bold", 20)).grid(column=0, row=0, columnspan=2, sticky='nw')
        tk.Label(export_obj_subframe, text="Change Output Directory:").grid(row=1, column=0)
        tk.Button(export_obj_subframe, text="Browse", command=self.change_output_directory).grid(row=1, column=1)
        self.output_dir_str = tk.StringVar(self.main, os.getcwd()[:12] + "..." + "/Data/User_data")
        tk.Label(export_obj_subframe, textvariable=self.output_dir_str).grid(row=2, column=0, columnspan=2)
        self.export_info = tk.BooleanVar(self.main, False)
        tk.Checkbutton(export_obj_subframe, text="Export Info", variable=self.export_info.get()).grid(row=3, column=0, columnspan=2)

        tk.Button(export_obj_subframe, text="Export", font=("Times New Roman bold", 20), command=self.export_selections)\
            .grid(row=4, column=0, columnspan=2, rowspan=2)

        self.main.mainloop()

    ############################################# Functions  ###########################################################

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
        self.set_sys_info()
        # Set the networks atoms
        self.sys.net.atoms = self.sys.atoms



    def set_sys_info(self):
        # We want to get the number of atoms, the number of molecules, etc
        mystr = str(len(self.sys.atoms)) + '\n' + str(len(self.sys.mols)) + '\n' + str(len(self.sys.residues)) + \
                "\n   ~   \n   ~   \n   ~   "
        # Set the variables
        self.sys_pros.set(mystr)

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

    def load_index(self):
        # File grabber pop up
        file_path = filedialog.askopenfilename()
        self.index_file = file_path

    def add_atoms_button(self):
        pass

    # Build network button method. Cements the settings and destroys the gui
    def build_network_button(self):
        if self.sys is not None:
            self.sys.build_network()
            self.load_frame.destroy()
            self.main.geometry("800x900")
            self.analysis_frame.pack()

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
        pass

    def undo_selection(self):
        pass

    def add_selection(self):
        pass

    def export_selections(self):
        self.exporting = True