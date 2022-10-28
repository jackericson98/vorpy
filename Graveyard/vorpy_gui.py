import os
import tkinter as tk
from tkinter import filedialog
from tkinter.messagebox import showinfo
from tkinter.ttk import Progressbar


class Vorpy:
    """Vorpy GUI class. When instantiated the Gui will launch"""
    def __init__(self, width=450, height=700):

        # Set up the window
        self.sys_atom_list = []
        self.vp_main = tk.Tk()
        self.width = width
        self.height = height
        self.vp_main.geometry(str(width) + "x" + str(height))
        self.vp_main.title('vorpy')
        self.vorpy_directory = os.getcwd()
        self.sys_output_directory = ""
        self.sys_file_address = None
        self.vert_file_address = None
        self.building = False

        # Set up the main frame
        f = tk.Frame(self.vp_main)
        f.place(relx=0.5, rely=0.5, anchor='c')

        # Header
        tk.Label(f, text="Vorpy", font=("Times New Roman bold", 40)).grid(row=0, column=0, columnspan=10)
        tk.Label(f, text="").grid(row=1)

        # Load file variables
        # Load header
        tk.Label(f, text="Load Atoms", font=('Times New Roman bold', 20))\
            .grid(row=2, column=0, columnspan=9, sticky='w')

        # System file information
        self.sys_file_name = tk.StringVar(self.vp_main, "")
        tk.Label(f, text="System: ", font=('Times New Roman', 12)).grid(row=3, column=0, columnspan=2, sticky='w')
        tk.Label(f, textvariable=self.sys_file_name, font=('Times New Roman', 15)) \
            .grid(row=3, column=1, columnspan=2, sticky='w')
        tk.Button(f, text="Load System ", command=self.load_sys_button).grid(row=3, column=9, sticky='e')

        # Vertices file information
        self.vert_file_name = tk.StringVar(self.vp_main, "")
        tk.Label(f, text="Vertices: ", font=('Times New Roman', 12)).grid(row=4, column=0, columnspan=3, sticky='w')
        tk.Label(f, textvariable=self.vert_file_name, font=('Times New Roman', 15)) \
            .grid(row=4, column=3, columnspan=9, sticky='w')
        tk.Button(f, text="Load Vertices", command=self.add_vertices).grid(row=4, column=9, sticky='e')

        # Build system information
        self.user_atoms_x, self.user_atoms_y = tk.StringVar(self.vp_main), tk.StringVar(self.vp_main)
        self.user_atoms_z, self.user_atoms_rad = tk.StringVar(self.vp_main), tk.StringVar(self.vp_main)
        self.user_atoms_list = tk.StringVar(self.vp_main)

        tk.Label(f, text="Add Atoms: ", font=("Times New Roman", 15)).grid(row=6, column=0, columnspan=2, sticky='w')
        tk.Label(f, textvariable=self.user_atoms_list).grid(row=6, column=1, columnspan=9)
        tk.Label(f, text="X:").grid(row=7, column=0, sticky='e')
        tk.Entry(f, textvariable=self.user_atoms_x, width=5).grid(row=7, column=1, sticky='w')
        tk.Label(f, text="Y:").grid(row=7, column=2, sticky='e')
        tk.Entry(f, textvariable=self.user_atoms_y, width=5).grid(row=7, column=3, sticky='w')
        tk.Label(f, text="Z:").grid(row=7, column=4, sticky='e')
        tk.Entry(f, textvariable=self.user_atoms_z, width=5).grid(row=7, column=5, sticky='w')
        tk.Label(f, text="Rad:").grid(row=7, column=6, sticky='e')
        tk.Entry(f, textvariable=self.user_atoms_rad, width=5).grid(row=7, column=7, sticky='w')
        tk.Button(f, text="Add Atom", command=self.add_atoms).grid(row=7, column=9, sticky='e')
        tk.Label(f, text="").grid(row=8)

        # Setting variables:
        # Settings header
        tk.Label(f, text='Settings', font=('Times New Roman bold', 20)).grid(row=9, column=0, columnspan=9, sticky='w')

        # System box size multiplier
        self.sys_box_x_flt = tk.DoubleVar(self.vp_main, 1.4)
        tk.Label(f, text="Box Size: ").grid(row=10, column=0, columnspan=3, sticky='ws')
        tk.Scale(f, from_=1.05, to=3, orient=tk.HORIZONTAL, variable=self.sys_box_x_flt, resolution=.05)\
            .grid(row=10, column=2, columnspan=5, sticky='ewn')
        tk.Label(f, text="x", font=('Times New Roman', 15)).grid(row=10, column=7, columnspan=2, sticky='s')

        # System resolution value
        self.sys_res_flt = tk.DoubleVar(self.vp_main, 0.1)
        tk.Label(f, text="Resolution: ").grid(row=11, column=0, columnspan=3, sticky='ws')
        tk.Scale(f, from_=0.05, to=1, orient=tk.HORIZONTAL, variable=self.sys_res_flt, resolution=0.05)\
            .grid(row=11, column=2, columnspan=5, sticky='ewn')
        tk.Label(f, text=u'\u212B', font=('Times New Roman', 15)).grid(row=11, column=7, columnspan=2, sticky='s')

        # Beta value for the system
        self.sys_alpha_value = tk.DoubleVar(self.vp_main, 5)
        tk.Label(f, text="Max Vertex: ").grid(row=12, column=0, columnspan=3, sticky='ws')
        tk.Scale(f, from_=2, to=10, orient=tk.HORIZONTAL, variable=self.sys_alpha_value, resolution=0.05) \
            .grid(row=12, column=2, columnspan=5, sticky='ewn')
        tk.Label(f, text=u'\u212B', font=('Times New Roman', 15)).grid(row=12, column=7, columnspan=2, sticky='s')

        # Voronota system check
        self.vta = tk.BooleanVar(self.vp_main)
        tk.Checkbutton(f, text="Flat Faces ", variable=self.vta, onvalue=True, offvalue=False)\
            .grid(row=10, column=9, sticky='ws')

        # Parallelize check
        self.parallelize = tk.BooleanVar(self.vp_main)
        tk.Checkbutton(f, text="Parallelize ", variable=self.parallelize, onvalue=True, offvalue=False)\
            .grid(row=11, column=9, sticky='ws')

        # Output variables
        # Output Header
        tk.Label(f, text="Outputs", font=("Times New Roman bold", 20)).grid(row=14, column=0, columnspan=9, sticky='w')

        # Output directory information
        self.output_dir_str = tk.StringVar(self.vp_main, os.getcwd()[:min(60, len(os.getcwd()))]
                                           + "/Data/User_data" + self.sys_file_name.get())
        tk.Label(f, text="Output Directory: ").grid(row=15, column=0, columnspan=3, sticky='w')
        tk.Label(f, textvariable=self.output_dir_str).grid(row=16, column=0, columnspan=10, sticky='w')
        tk.Button(f, text="     Change    ", command=self.change_output_directory)\
            .grid(row=15, column=8, columnspan=2, sticky='en')

        # Check boxes
        self.output_all = tk.BooleanVar(self.vp_main, True)
        tk.Checkbutton(f, text="All", variable=self.output_all, onvalue=True, offvalue=False)\
            .grid(row=17, column=0, columnspan=3, sticky='w')
        self.output_sys = tk.BooleanVar(self.vp_main, False)
        tk.Checkbutton(f, text="System", variable=self.output_sys, onvalue=True, offvalue=False) \
            .grid(row=17, column=4, columnspan=3, sticky='w')
        self.output_pdb = tk.BooleanVar(self.vp_main, False)
        tk.Checkbutton(f, text="System pdb ", variable=self.output_pdb, onvalue=True, offvalue=False) \
            .grid(row=17, column=8, columnspan=3, sticky='w')
        self.output_atoms = tk.BooleanVar(self.vp_main, False)
        tk.Checkbutton(f, text="Atoms", variable=self.output_atoms, onvalue=True, offvalue=False) \
            .grid(row=18, column=0, columnspan=3, sticky='w')
        self.output_surfs = tk.BooleanVar(self.vp_main, False)
        tk.Checkbutton(f, text="Surfaces", variable=self.output_surfs, onvalue=True, offvalue=False) \
            .grid(row=18, column=4, columnspan=3, sticky='w')
        self.output_mols = tk.BooleanVar(self.vp_main, False)
        tk.Checkbutton(f, text="Molecule", variable=self.output_mols, onvalue=True, offvalue=False) \
            .grid(row=18, column=8, columnspan=3, sticky='w')
        self.output_residues = tk.BooleanVar(self.vp_main, False)
        tk.Checkbutton(f, text="Residues ", variable=self.output_residues, onvalue=True, offvalue=False) \
            .grid(row=19, column=0, columnspan=3, sticky='w')
        self.output_verts = tk.BooleanVar(self.vp_main, False)
        tk.Checkbutton(f, text="Vertices ", variable=self.output_verts, onvalue=True, offvalue=False) \
            .grid(row=19, column=4, columnspan=3, sticky='w')
        tk.Label(f, text="").grid(row=20)

        # Run program
        tk.Button(f, text="Build System", font=("Times New Roman bold", 20), command=self.build_network_button) \
            .grid(row=21, column=0, columnspan=10)

        # End the loop
        self.vp_main.mainloop()

    # Load system button function. Calls the file browser and sets the system
    def load_sys_button(self):
        # File grabber pop up
        file_path = filedialog.askopenfilename()
        # Set the file path
        if file_path:
            self.sys_file_address = file_path
        # Get the file name
        filename = ""
        i = -1
        if len(self.sys_file_address) > 0:
            while self.sys_file_address[i] != "/":
                filename = filename + self.sys_file_address[i]
                i -= 1
        else:
            filename = "    atad_resU"
        # Set the file name
        self.sys_file_name.set(filename[::-1][:-4])

    # Add atoms method. When this command runs it updates the system atoms list
    def add_atoms(self):
        # Add the values as a list of atoms
        self.sys_atom_list.append([[float(self.user_atoms_x.get()), float(self.user_atoms_y.get()),
                                    float(self.user_atoms_z.get())], float(self.user_atoms_rad.get())])
        # Print the first and last two atoms in the list
        if len(self.sys_atom_list) <= 4:
            self.user_atoms_list.set(self.sys_atom_list)
        else:
            self.user_atoms_list.set(str(self.sys_atom_list[:2]) + " ... " + str(self.sys_atom_list[-2:]))


    # Add vertices method. Adds the vertices in the file to the system
    def add_vertices(self):
        # File grabber pop up
        file_path = filedialog.askopenfilename()
        if file_path is None:
            return
        self.vert_file_address = file_path
        # Get the file name
        filename = ""
        i = -1
        if self.sys_file_address and len(self.sys_file_address) > 0:
            while self.sys_file_address[i] != "/":
                filename = filename + self.sys_file_address[i]
                i -= 1
        else:
            # Written in reverse
            filename = "    strev_atad_resU"
        if file_path:
            self.vert_file_name.set(filename[::-1][:-4] + "_verts")
        else:
            self.vert_file_address = None

    # Change output directory method. Updates the location of the output directory
    def change_output_directory(self):
        # File grabber pop up
        file_path = filedialog.askdirectory()
        # Create the System
        if file_path:
            self.sys_output_directory = file_path
            self.output_dir_str.set(file_path[:12] + ' ... ' + file_path[-12:])
        else:
            ErrorBox("Directory not changed")

    # Build network button method. Cements the settings and destroys the gui
    def build_network_button(self):
        self.vp_main.destroy()


class ErrorBox:
    """Error box class. Used to indicate if an error has occurred and takes the message as an input"""
    def __init__(self, error_message):
        self.root = tk.Tk()
        self.root.geometry("300x75")

        self.error_message = tk.StringVar(self.root, error_message)
        self.root.title("vorpy")
        error_frame = tk.Frame(self.root)
        error_frame.place(relx=0.5, rely=0.45, anchor='c')

        tk.Label(error_frame, text="Error!", font=("Times New Roman bold", 15)).grid(row=0, column=0, sticky='w')
        tk.Label(error_frame, textvariable=self.error_message).grid(row=1, column=0, sticky='e')
        tk.Button(error_frame, text="OK", command=self.root.quit()).grid(row=2, column=0, sticky='e')

        # End the loop
        self.root.mainloop()


class LoadingBox:
    """Loading box gui"""
    def __init__(self):
        # root window
        self.root = tk.Tk()
        self.root.geometry('300x120')
        self.root.title('Processing System')
        # progressbar
        self.pb = Progressbar(self.root, orient='horizontal', mode='determinate', length=280)
        # place the progressbar
        self.pb.grid(column=0, row=0, columnspan=2, padx=10, pady=20)

        # label
        self.value_label = tk.Label(self.root, text=self.update_progress_label())
        self.value_label.grid(column=0, row=1, columnspan=2)

        # start button
        self.start_button = tk.Button(self.root, text='Progress', command=self.progress)
        self.start_button.grid(column=0, row=2, padx=10, pady=10, sticky=tk.E)

        self.stop_button = tk.Button(self.root, text='Stop', command=self.stop)
        self.stop_button.grid(column=1, row=2, padx=10, pady=10, sticky=tk.W)

        self.root.mainloop()

    def update_progress_label(self):
        return f"Current Progress: {self.pb['value']}%"

    def progress(self):
        if self.pb['value'] < 100:
            self.pb['value'] += 20
            self.value_label['text'] = self.update_progress_label()
        else:
            showinfo(message='The progress completed!')

    def stop(self):
        self.pb.stop()
        self.value_label['text'] = self.update_progress_label()


Vorpy()



################################################# Load GUI #############################################################

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

# tk.Label(self.choose_index_subframe, text="Create Index", font=("Times New Roman bold", 12))\
#     .grid(row=3, column=0, columnspan=3)
# self.current_atom_selection = tk.StringVar(self.main, "")
# tk.Label(self.choose_index_subframe, text="Atom").grid(row=4, column=0)
# self.atom_options = tk.OptionMenu(self.choose_index_subframe, self.current_atom_selection, "", *self.sys.atoms)
# self.atom_options.grid(row=5, column=0)
#
# # Molecules
# self.current_mol_selection = tk.StringVar(self.main)
# tk.Label(self.choose_index_subframe, text="Molecule").grid(row=4, column=1)
# self.mol_options = tk.OptionMenu(self.choose_index_subframe, self.current_mol_selection, "", *self.sys.mols)
# self.mol_options.grid(row=5, column=1)
#
# # Residues
# self.current_res_selection = tk.StringVar(self.main, "")
# tk.Label(self.choose_index_subframe, text="Residue").grid(row=4, column=2)
# self.res_options = tk.OptionMenu(self.choose_index_subframe, self.current_res_selection, "", *self.sys.residues)
# self.res_options.grid(row=5, column=2)
#
#
#
#
# # Separator
# ttk.Separator(self.analysis_frame).grid(row=1, columnspan=3, sticky='ew')
#
# # Network name
# tk.Label(self.analysis_frame, textvariable=self.net_name, font=("Times New Roman bold", 20))\
#     .grid(row=2, columnspan=3, padx=10, pady=10)
#
# # Separator
# ttk.Separator(self.analysis_frame).grid(row=3, columnspan=3, sticky='ew')
#
# # Network information List
# anal_net_frame = tk.Frame(self.analysis_frame)
# anal_net_frame.grid(row=4, column=0, padx=10, pady=10)
#
# tk.Label(anal_net_frame, text="Network Information:", font=("Times New Roman bold", 20))\
#     .grid(row=0, column=0, columnspan=2, sticky='nw')
#
# ttk.Separator(anal_net_frame, orient=tk.VERTICAL).grid(row=1, rowspan=3, column=1, sticky='ns')
# ttk.Separator(anal_net_frame, orient=tk.HORIZONTAL).grid(row=2, columnspan=3, column=0, sticky='ew')
#
# # Create the system information sub frame
# anal_sys_info_subframe = tk.Frame(anal_net_frame)
# anal_sys_info_subframe.grid(row=1, column=0, padx=5, pady=5)
#
# tk.Label(anal_sys_info_subframe, text="System Information").grid(row=3, columnspan=3, sticky='n')
# tk.Label(anal_sys_info_subframe, text="Atoms:\nMolecules:\nResidues:\n    ~    \n    ~    \n    ~    \n") \
#     .grid(row=4, column=0)
# tk.Label(anal_sys_info_subframe, text=":\n:\n:\n:\n:\n:").grid(row=4, column=1)
# tk.Label(anal_sys_info_subframe, textvariable=self.sys_pros).grid(row=4, column=2)
#
# # Create the system information sub frame
# anal_cpu_info_subframe = tk.Frame(anal_net_frame)
# anal_cpu_info_subframe.grid(row=1, column=2, padx=5, pady=5)
#
# tk.Label(anal_cpu_info_subframe, text="CPU Information").grid(row=3, columnspan=3, sticky='n')
# tk.Label(anal_cpu_info_subframe, text="My Time:\nCPU Time:\nComputations:\n    ~    \n    ~    \n    ~    \n") \
#     .grid(row=4, column=0)
# tk.Label(anal_cpu_info_subframe, text=":\n:\n:\n:\n:\n:").grid(row=4, column=1)
# self.cpu_info = tk.StringVar(self.main, "   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n")
# tk.Label(anal_cpu_info_subframe, textvariable=self.cpu_info).grid(row=4, column=2)
#
# # Create the system information sub frame
# anal_settings_subframe = tk.Frame(anal_net_frame)
# anal_settings_subframe.grid(row=3, column=0, padx=5, pady=5)
#
# tk.Label(anal_settings_subframe, text="Settings Information").grid(row=3, columnspan=3, sticky='n')
# tk.Label(anal_settings_subframe,
#          text="Box Size:\nResolution:\nMax Vertex:\n    ~    \n    ~    \n    ~    \n").grid(row=4, column=0)
# tk.Label(anal_settings_subframe, text=":\n:\n:\n:\n:\n:").grid(row=4, column=1)
# tk.Label(anal_settings_subframe, textvariable=self.net_sets).grid(row=4, column=2)
#
# # Create the system information sub frame
# anal_outputs_info_subframe = tk.Frame(anal_net_frame)
# anal_outputs_info_subframe.grid(row=3, column=2, padx=5, pady=5)
#
# tk.Label(anal_outputs_info_subframe, text="Output Information").grid(row=3, columnspan=3, sticky='n')
# tk.Label(anal_outputs_info_subframe, text="Vertices:\nSurfaces:\nEdges:\nFull Cells\n   ~   \n   ~   ") \
#     .grid(row=4, column=0)
# tk.Label(anal_outputs_info_subframe, text=":\n:\n:\n:\n:\n:").grid(row=4, column=1)
# self.output_info = tk.StringVar(self.main, "   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n")
# tk.Label(anal_outputs_info_subframe, textvariable=self.output_info).grid(row=4, column=2)
#
# # Seperator
# ttk.Separator(self.analysis_frame, orient=tk.VERTICAL).grid(row=3, column=1, rowspan=3, sticky='ns')
#
# # Export information List
# group_info_frame = tk.Frame(self.analysis_frame)
# group_info_frame.grid(row=4, column=2, sticky='ns', padx=10, pady=10)
#
# tk.Label(group_info_frame, text="Group Information:", font=("Times New Roman bold", 18))\
#     .grid(row=0, column=0, sticky='nw', columnspan=5)
# tk.Label(group_info_frame, text="Group 1").grid(row=1, column=2)
# ttk.Separator(group_info_frame, orient=tk.VERTICAL).grid(row=1, rowspan=3, column=1, sticky='ns')
# ttk.Separator(group_info_frame, orient=tk.VERTICAL).grid(row=1, rowspan=3, column=3, sticky='ns')
#
# tk.Label(group_info_frame, text="Group 2").grid(row=1, column=4)
# ttk.Separator(group_info_frame).grid(row=2, columnspan=5, sticky='ew')
# tk.Label(group_info_frame, text="Surface Area:\n\nVolume:\n\nBoundary atoms:\n\nOuter Atoms:\n")\
#     .grid(row=3, column=0, sticky='w')
# self.g1_info = tk.StringVar(self.main, "   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n   ~   ")
# tk.Label(group_info_frame, textvariable=self.g1_info).grid(row=3, column=2)
# self.g2_info = tk.StringVar(self.main, "   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n   ~   \n   ~   ")
# tk.Label(group_info_frame, textvariable=self.g2_info).grid(row=3, column=4)
# ttk.Separator(group_info_frame).grid(row=4, column=0, columnspan=5, sticky='ew')
# tk.Label(group_info_frame, text="Interface:").grid(row=6, column=1)
# tk.StringVar()
# tk.Label(group_info_frame, text="   ~   \n   ~   \n   ~   ")
#
# self.interface_info = tk.StringVar(self.main, )
#
# ttk.Separator(self.analysis_frame).grid(row=5, columnspan=5, sticky='ew')
#
# # Main frame Network object gathering
# export_obj_frame = tk.Frame(self.analysis_frame)
# export_obj_frame.grid(row=6, column=0, columnspan=3, padx=10, pady=10)
#
# # Group sub frame
# self.using_group1 = tk.BooleanVar(self.main, True)
# self.using_group2 = tk.BooleanVar(self.main, not self.using_group1.get())
# anal_check_subframe = tk.Frame(export_obj_frame)
# anal_check_subframe.grid(column=0, row=0, padx=10, pady=10)
#
# tk.Label(anal_check_subframe, text="Groups:", font=("Times New Roman bold", 20)).grid(row=0)
# self.group1_selections = tk.StringVar(self.main, "")
# self.group2_selections = tk.StringVar(self.main, "")
# tk.Checkbutton(anal_check_subframe, text="Group 1", font=("Times New Roman", 15), variable=self.using_group1,
#                command=self.flip_g2).grid(row=1, column=0)
# tk.Label(anal_check_subframe, text="Selections").grid(row=2, column=0)
# tk.Label(anal_check_subframe, textvariable=self.group1_selections).grid(row=3, column=0)
# tk.Checkbutton(anal_check_subframe, text="Group 2", font=("Times New Roman", 15), variable=self.using_group2,
#                command=self.flip_g1).grid(row=4, column=0)
# tk.Label(anal_check_subframe, text="Selections").grid(row=5, column=0)
# tk.Label(anal_check_subframe, textvariable=self.group2_selections).grid(row=6, column=0)
#
# # Choose index
#
# # Seperator
# ttk.Separator(export_obj_frame, orient=tk.VERTICAL).grid(column=1, row=0, sticky='ns', rowspan=3)
# # Choose index sub frame
# self.choose_index_subframe = tk.Frame(export_obj_frame)
# self.choose_index_subframe.grid(column=2, row=0, sticky='ns', padx=10, pady=10)
#
# # Set the label header for the "choose index" sub frame
# tk.Label(self.choose_index_subframe, text="Get Atoms:", font=("Times New Roman bold", 20))\
#     .grid(row=0, columnspan=3, sticky='nw')
# tk.Label(self.choose_index_subframe, text="Choose Index:", font=("Times New Roman bold", 15))\
#     .grid(row=1, column=0)
# self.current_ndx = tk.StringVar(self.main, "[None]")
# tk.OptionMenu(self.choose_index_subframe, variable=self.current_ndx, value=self.index_list)\
#     .grid(row=1, column=1)
# tk.Button(self.choose_index_subframe, text="Browse", command=self.load_index).grid(row=1, column=2)
#
# # Create the index
# # Atoms
# tk.Label(self.choose_index_subframe, text="Create Index", font=("Times New Roman bold", 12))\
#     .grid(row=3, column=0, columnspan=3)
# self.current_atom_selection = tk.StringVar(self.main, "")
# tk.Label(self.choose_index_subframe, text="Atom").grid(row=4, column=0)
# self.atom_options = tk.OptionMenu(self.choose_index_subframe, self.current_atom_selection, "", *self.sys.atoms)
# self.atom_options.grid(row=5, column=0)
#
# # Molecules
# self.current_mol_selection = tk.StringVar(self.main)
# tk.Label(self.choose_index_subframe, text="Molecule").grid(row=4, column=1)
# self.mol_options = tk.OptionMenu(self.choose_index_subframe, self.current_mol_selection, "", *self.sys.mols)
# self.mol_options.grid(row=5, column=1)
#
# # Residues
# self.current_res_selection = tk.StringVar(self.main, "")
# tk.Label(self.choose_index_subframe, text="Residue").grid(row=4, column=2)
# self.res_options = tk.OptionMenu(self.choose_index_subframe, self.current_res_selection, "", *self.sys.residues)
# self.res_options.grid(row=5, column=2)
#
# # Selection Button
# tk.Button(self.choose_index_subframe, text="Reset Group", command=self.reset_group).grid(row=6, column=0)
# tk.Button(self.choose_index_subframe, text="Undo Selection", command=self.undo_selection).grid(row=6, column=1)
# tk.Button(self.choose_index_subframe, text="Add Selection", command=self.add_selection).grid(row=6, column=2)
#
# # Seperator
# ttk.Separator(export_obj_frame, orient=tk.VERTICAL).grid(column=3, row=0, sticky='ns', rowspan=3)
#
# # Export object sub frame
# export_obj_subframe = tk.Frame(export_obj_frame)
# export_obj_subframe.grid(column=4, row=0, padx=10, pady=10, sticky='ns')
#
# tk.Label(export_obj_subframe, text="Export:", font=("Times New Roman bold", 20))\
#     .grid(column=0, row=0, columnspan=2, sticky='nw')
# tk.Label(export_obj_subframe, text="Change Output Directory:").grid(row=1, column=0)
# tk.Button(export_obj_subframe, text="Browse", command=self.change_output_directory).grid(row=1, column=1)
# self.output_dir_str = tk.StringVar(self.main, os.getcwd()[:12] + "..." + "/Data/User_data")
# tk.Label(export_obj_subframe, textvariable=self.output_dir_str).grid(row=2, column=0, columnspan=2)
# self.export_info = tk.BooleanVar(self.main, False)
# tk.Checkbutton(export_obj_subframe, text="Export Info", variable=self.export_info)\
#     .grid(row=3, column=0, columnspan=2)
#
# tk.Button(export_obj_subframe, text="Export", font=("Times New Roman bold", 20),
#           command=self.export_selections).grid(row=4, column=0, columnspan=2, rowspan=2)


# Update progressbar function. Updates the color of the progressbar depending on the stage that vorpy is in
def update_progress_canvas(self, task_ndx=None):
    # Set up the canvas variable
    canvas = self.progress_canvas
    # Set up the dimensions
    w, h, = int(self.canvas_width / 4), int(self.canvas_height / 2)

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
    canvas.create_polygon(0, 0, w - h / 2, 0, w + h / 2, h, w - h / 2, 2 * h, 0, 2 * h, fill=colors[0], outline='black')
    canvas.create_polygon(w - h / 2, 0, 2 * w - h / 2, 0, 2 * w + h / 2, h, 2 * w - h / 2, 2 * h, w - h / 2, 2 * h,
                          w + h / 2, h, fill=colors[1],
                          outline='black')
    canvas.create_polygon(2 * w - h / 2, 0, 3 * w - h / 2, 0, 3 * w + h / 2, h, 3 * w - h / 2, 2 * h, 2 * w - h / 2,
                          2 * h, 2 * w + h / 2, h,
                          fill=colors[2], outline='black')
    canvas.create_polygon(3 * w - h / 2, 0, 4 * w, 0, 4 * w, 2 * h, 3 * w - h / 2, 2 * h, 3 * w + h / 2, h,
                          fill=colors[3], outline='black')
    # Add the text on the polygons
    canvas.create_text(0.5 * w, h, text="Finding\nVertices", font=("Times New Roman bold", 10))
    canvas.create_text(1.5 * w, h, text="Connecting\nNetwork", font=("Timescv New Roman bold", 10))
    canvas.create_text(2.5 * w, h, text="Building\nSurfaces", font=("Times New Roman bold", 10))
    canvas.create_text(3.5 * w, h, text="Analyzing\nNetwork", font=("Times New Roman bold", 10))




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