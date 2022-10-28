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
        # self.load_frame.pack()

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
        tk.Button(load_system_frame, text="Add Atoms ", command=self.build_atoms_button)\
            .grid(row=2, column=0, sticky='w')
        tk.Button(load_system_frame, text="Load System ", command=self.load_sys_button).grid(row=2, column=1)
        tk.Button(load_system_frame, text="Load Frames ", command=self.load_frames_button)\
            .grid(row=2, column=2, sticky='e')
        # Add the system information header and labels
        tk.Label(load_system_frame, text="System Information").grid(row=3, column=1, sticky='n')
        tk.Label(load_system_frame, text="Atoms:\nMolecules:\nSolute:\n   ~   \n   ~   \n   ~   \n")\
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

        """######################################### Analysis Frame #################################################"""

        # Set up the analysis frame
        self.analysis_frame = tk.Frame(self.main, name="analysis")
        self.analysis_frame.pack()

        # Header
        tk.Label(self.analysis_frame, text="VorPy", font=("Times New Roman bold", 20), justify=tk.CENTER)\
            .grid(row=0, columnspan=2, padx=10, pady=10)

        # System information subframe
        self.sys_info_frm = tk.Frame(self.analysis_frame, highlightthickness=1, highlightcolor='blue')
        self.sys_info_frm.grid(row=1, column=0)
        # Header
        tk.Label(self.sys_info_frm, text="System", font=("bold underlined", 15), justify=tk.CENTER)\
            .grid(column=0, columnspan=3)

        # Network information subframe
        self.net_info_frm = tk.Frame(self.analysis_frame, highlightthickness=1)
        self.net_info_frm.grid(row=1, column=1)
        # Header
        tk.Label(self.net_info_frm, text="Network", font=("Times New Roman bold underlined", 15), justify=tk.CENTER)\
            .grid(column=0, columnspan=3)

        # Select atoms subframe
        self.select_atoms_frm = tk.Frame(self.analysis_frame)
        self.select_atoms_frm.grid(row=2, columnspan=2)

        # Header
        tk.Label(self.select_atoms_frm, text="Select Atoms", font=("underlined bold", 20)).grid(columnspan=3)

        # Current selection textbox
        self.select_atoms_current_selection_frm = tk.Frame(self.select_atoms_frm)
        self.select_atoms_current_selection_frm.grid(row=1, column=0)
        tk.Label(self.select_atoms_current_selection_frm, text="Current Selection").grid(row=0, column=0)
        self.selected_atoms_str = tk.StringVar(self.main)
        tk.Label(self.select_atoms_current_selection_frm, textvariable=self.selected_atoms_str, highlightcolor='white',
                 highlightthickness=1, justify=tk.CENTER).grid(row=1, column=0, sticky='ns')

        # Select atoms subsub-frame
        self.select_atoms_subfrm = tk.Frame(self.select_atoms_frm)
        self.select_atoms_subfrm.grid(row=1, column=1, sticky='nsew')

        # Molecule Dropdown
        self.sys.mol_names, self.sys.res_names, self.sys.atom_names = [["o", "sdds", "s", "G"] for _ in range(3)]
        tk.Label(self.select_atoms_subfrm, text="Molecules").grid()
        self.current_mol_selection = tk.StringVar(self.main)
        self.mol_options = tk.OptionMenu(self.select_atoms_subfrm, self.current_mol_selection, "", *self.sys.mol_names)
        self.mol_options.grid(row=1)
        tk.Button(self.select_atoms_subfrm, text="Add", command=self.add_mol_button).grid(row=1, column=1)

        tk.Label(self.select_atoms_subfrm, text="Residues").grid(row=2)
        self.current_res_selection = tk.StringVar(self.main)
        self.res_options = tk.OptionMenu(self.select_atoms_subfrm, self.current_res_selection, "", *self.sys.res_names)
        self.res_options.grid(row=3)
        tk.Button(self.select_atoms_subfrm, text="Add", command=self.add_res_button).grid(row=3, column=1)

        tk.Label(self.select_atoms_subfrm, text="Atoms").grid()
        self.current_atom_selection = tk.StringVar(self.main)
        self.atom_options = tk.OptionMenu(self.select_atoms_subfrm, self.current_atom_selection, "", *self.sys.atom_names)
        self.atom_options.grid(row=5)
        tk.Button(self.select_atoms_subfrm, text="Add", command=self.add_atom_button).grid(row=5, column=1)

        # Select atoms buttons subsub-frame
        self.select_atoms_buttons_subfrm = tk.Frame(self.select_atoms_frm)
        self.select_atoms_buttons_subfrm.grid(column=2, row=1)

        tk.Button(self.select_atoms_buttons_subfrm, text="Load Index", command=self.load_index)\
            .grid(row=0, padx=10, pady=10)
        tk.Button(self.select_atoms_buttons_subfrm, text="Undo Last", command=self.undo_last)\
            .grid(row=1, padx=10, pady=10)
        tk.Button(self.select_atoms_buttons_subfrm, text="Reset All", command=self.reset_all)\
            .grid(row=2, padx=10, pady=10)

        # Export selections
        self.export_selections_subfrm = tk.Frame(self.analysis_frame)
        self.export_selections_subfrm.grid(row=3, columnspan=2)

        # Header
        tk.Label(self.export_selections_subfrm, text="Export Selections", font=("underlined bold", 20))\
            .grid(columnspan=3)

        # Export cell subframe
        self.export_cell_subfrm = tk.Frame(self.export_selections_subfrm)
        self.export_cell_subfrm.grid(row=1)
        self.cell_atoms = []
        self.cell_atoms_names = tk.StringVar(self.main)
        tk.Label(self.export_cell_subfrm, textvariable=self.cell_atoms_names).grid()

        # Export interface subframe
        self.export_interface_subfrm = tk.Label(self.export_selections_subfrm)
        self.export_interface_subfrm.grid(row=1, column=1)


        # Export information subframe
        self.export_info_subfrm = tk.Label(self.export_selections_subfrm)
        self.export_info_subfrm.grid(row=1, column=2)

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
        self.mol_options = tk.OptionMenu(self.select_atoms_subfrm, self.current_mol_selection, *self.sys.mol_names)
        self.mol_options.grid(row=5, column=1)

        self.res_options.destroy()
        self.res_options = tk.OptionMenu(self.select_atoms_subfrm, self.current_res_selection, *self.sys.res_names)
        self.res_options.grid(row=5, column=2)

        self.atom_options.destroy()
        self.atom_options = tk.OptionMenu(self.select_atoms_subfrm, self.current_atom_selection, *self.sys.atom_names)
        self.atom_options.grid(row=5, column=0)

    # Load frames function.
    def load_frames_button(self):
        # File grabber pop up
        file_path = filedialog.askopenfilename()
        self.frame_files.append(file_path)

    def build_atoms_button(self):
        pass

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
        self.sys.set_output_directory()

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



    def reset_all(self):
        self.selected_atoms = []
        self.selected_atoms_names = []
        self.selected_atoms_last = 0
        self.selected_atoms_str.set("")

    def undo_last(self):
        self.selected_atoms.pop()
        self.selected_atoms_names.pop()
        self.selected_atoms_str.set("\n".join(self.selected_atoms_names))

    def add_mol_button(self):
        self.selected_atoms.append(self.sys.mols[self.sys.mol_names.index(self.current_mol_selection.get()) - 1])

    def add_res_button(self):
        self.selected_atoms.append(self.sys.residues[self.sys.res_names.index(self.current_res_selection.get()) - 1])

    def add_atom_button(self):
        self.selected_atoms.append(self.sys.atoms[self.sys.atom_names.index(self.current_atom_selection.get()) - 1])


Vorpy()