import customtkinter as ctk
import tkinter as tk
import tkinter.ttk as ttk
from tkinter import filedialog
from System.system import System, Network
from System.output import *
from System.group import *


class Vorpy(ctk.CTk):
    def __init__(self):
        # Get all the characteristics from the CTk class
        super().__init__()

        self.sys = None
        self.sys_file = None
        self.net_file = None
        self.g1 = None
        self.g2 = None

        self.elms = ['h', 'he', 'li', 'be', 'b', 'c', 'n', 'o', 'f', 'ne', 'na', 'mg', 'al', 'si', 'p', 's', 'cl', 'ar',
             'k' , 'ca', 'sc', 'ti', 'v' , 'cr', 'mn', 'fe', 'co', 'ni', 'cu', 'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr',
             'rb', 'sr', 'y' , 'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', 'cd', 'in', 'sn', 'sb', 'te', 'i' , 'xe',
             'cs', 'ba', 'la', 'hf', 'ta', 'w' , 're', 'os', 'ir', 'pt', 'au', 'hg', 'tl', 'pb', 'bi', 'po', 'at', 'rn',
             'fr', 'ra', 'ac', 'rf', 'db', 'sg', 'bh', 'hs', 'mt', 'ds', 'rg', 'cn', 'nh', 'fl', 'mc', 'lv', 'ts', 'og',
             'ce', 'pr', 'nd', 'pm', 'sm', 'eu', 'gd', 'tb', 'dy', 'ho', 'er', 'tm', 'yb', 'lu',
             'th', 'pa', 'u' , 'np', 'pu', 'am', 'cm', 'bk', 'cf', 'es', 'fm', 'md', 'no', 'lr']

        # Set the geometry
        self.geometry("700x700")
        self.title("vorPy")
        self.minsize(700, 700)

        # Set up the font variables

        # Main Header frame
        self.header_frame = ctk.CTkFrame(self)
        self.header_frame.grid(columnspan=4, sticky='nsew')
        ctk.CTkLabel(self.header_frame, text="Vorpy", anchor="center").grid(sticky='e')
        # System name header
        self.sys_name = ctk.StringVar(self, "System")
        ctk.CTkLabel(self.header_frame, textvariable=self.sys_name, text_font=("bold", 20)).grid(row=0, column=1, sticky='e')



        # Seperator
        ttk.Separator(self, orient=tk.HORIZONTAL).grid(row=2, columnspan=4, sticky='ew')

        # Data frame
        self.data_frame = ctk.CTkFrame(self)
        self.data_frame.grid(row=3, sticky='nsew')

        # Data Header
        ctk.CTkLabel(self.data_frame, text="Data").grid()

        # Seperator
        ttk.Separator(self.data_frame, orient=tk.HORIZONTAL).grid(row=1, sticky='ew')

        # System information
        self.sys_data_subfrm = ctk.CTkFrame(self.data_frame)
        self.sys_data_subfrm.grid(row=2, padx=5, pady=5)

        # System information header
        ctk.CTkLabel(self.sys_data_subfrm, text="System:").grid(columnspan=3)
        # Add the system information header and labels
        tk.Label(self.sys_data_subfrm, text="Atoms\nMolecules\nSolute\n   ~   \n   ~   \n   ~  ") \
            .grid(row=1, column=0, sticky='w')
        tk.Label(self.sys_data_subfrm, text=":\n:\n:\n:\n:\n:").grid(row=1, column=1)
        # Add the system information to match the labels
        self.sys_data = ctk.StringVar(self, "~\n~\n~\n~\n~\n~")
        tk.Label(self.sys_data_subfrm, textvariable=self.sys_data).grid(row=1, column=2)

        # Seperator
        ttk.Separator(self.data_frame, orient=tk.HORIZONTAL).grid(row=3, sticky='ew')

        # Network information
        self.net_data_subfrm = ctk.CTkFrame(self.data_frame)
        self.net_data_subfrm.grid(row=4, padx=5, pady=5)

        # System information header
        ctk.CTkLabel(self.net_data_subfrm, text="Network:").grid(columnspan=3)
        # Add the system information header and labels
        tk.Label(self.net_data_subfrm, text="Surfaces\nEdges\nVertices\nAtom/Cell\n   ~   \n   ~  ") \
            .grid(row=1, column=0, sticky='w')
        tk.Label(self.net_data_subfrm, text=":\n:\n:\n:\n:\n:").grid(row=1, column=1)
        # Add the system information to match the labels
        self.net_data = ctk.StringVar(self, "~\n~\n~\n~\n~\n~")
        tk.Label(self.net_data_subfrm, textvariable=self.net_data).grid(row=1, column=2)

        # Seperator
        ttk.Separator(self.data_frame, orient=tk.HORIZONTAL).grid(row=5, sticky='ew')

        # Network information
        self.build_data_subfrm = ctk.CTkFrame(self.data_frame)
        self.build_data_subfrm.grid(row=6, padx=5, pady=5)

        # System information header
        ctk.CTkLabel(self.build_data_subfrm, text="Build:").grid(columnspan=3)
        # Add the system information header and labels
        tk.Label(self.build_data_subfrm, text="Time\nCPU Time\nResolution\nBox Size\nMax Vert Rad\n   ~  ") \
            .grid(row=1, column=0, sticky='w')
        tk.Label(self.build_data_subfrm, text=":\n:\n:\n:\n:\n:").grid(row=1, column=1)
        # Add the system information to match the labels
        self.build_data = ctk.StringVar(self, "~\n~\n~\n~\n~\n~")
        tk.Label(self.build_data_subfrm, textvariable=self.build_data).grid(row=1, column=2)

        # Seperator
        ttk.Separator(self, orient=tk.VERTICAL).grid(row=3, column=1, sticky='ns')

        # View frame
        self.view_frame = ctk.CTkFrame(self)
        self.view_frame.grid(row=3, column=3, sticky='nsew')


        #

        ############################################### Load Frame #####################################################

        self.load_frame = ctk.CTkFrame(self.view_frame)
        self.load_frame.grid(row=1, columnspan=3, sticky='nesw', padx=20, pady=20)

        # System sub frame
        self.load_sys_frame = ctk.CTkFrame(self.load_frame)
        self.load_sys_frame.grid(row=0, column=0, padx=10, pady=10)

        # System Header
        ctk.CTkLabel(self.load_sys_frame, text="System").grid(row=0, column=0, columnspan=3)

        # Load system stuff
        ctk.CTkLabel(self.load_sys_frame, text="System File").grid(row=1, column=0)
        self.sys_file_str = ctk.StringVar(self)
        ctk.CTkLabel(self.load_sys_frame, textvariable=self.sys_file_str, width=150).grid(row=1, column=1)
        ctk.CTkButton(self.load_sys_frame, text="Browse", command=self.load_sys_button).grid(row=1, column=2, padx=4, pady=4)

        # Load index file stuff
        ctk.CTkLabel(self.load_sys_frame, text="Index File").grid(row=2, column=0)
        self.ndx_file_str = ctk.StringVar(self)
        ctk.CTkLabel(self.load_sys_frame, textvariable=self.ndx_file_str, width=150).grid(row=2, column=1)
        ctk.CTkButton(self.load_sys_frame, text="Browse", command=self.load_ndx_button).grid(row=2, column=2, padx=4, pady=4)

        # System sub frame
        self.load_net_frame = ctk.CTkFrame(self.load_frame)
        self.load_net_frame.grid(row=1, column=0, padx=10, pady=10)

        # System Header
        ctk.CTkLabel(self.load_net_frame, text="Network").grid(row=0, column=0, columnspan=3)

        # Load system stuff
        ctk.CTkLabel(self.load_net_frame, text="Network File").grid(row=1, column=0)
        self.net_file_str = ctk.StringVar(self)
        ctk.CTkLabel(self.load_net_frame, textvariable=self.net_file_str, width=150).grid(row=1, column=1)
        ctk.CTkButton(self.load_net_frame, text="Browse", command=self.load_net_button).grid(row=1, column=2, padx=4,
                                                                                             pady=4)

        # Load index file stuff
        ctk.CTkLabel(self.load_net_frame, text="Vertex File").grid(row=2, column=0)
        self.vert_file_str = ctk.StringVar(self)
        ctk.CTkLabel(self.load_net_frame, textvariable=self.vert_file_str, width=150).grid(row=2, column=1)
        ctk.CTkButton(self.load_net_frame, text="Browse", command=self.load_verts_button).grid(row=2, column=2, padx=4,
                                                                                             pady=4)


        ############################################### Build frame ####################################################

        self.build_frame = ctk.CTkFrame(self.view_frame)
        self.build_frame.grid(row=1, columnspan=3, sticky='nsew', padx=20, pady=20)

        # Loaded information frame
        self.loaded_subfrm = ctk.CTkFrame(self.build_frame)
        self.loaded_subfrm.grid(row=0, column=0, padx=10, pady=10, columnspan=2)

        # Loaded subframe header
        ctk.CTkLabel(self.loaded_subfrm, text="Loaded").grid(row=0, column=0, columnspan=3)

        # Files information subframe
        self.loaded_files_subfrm = ctk.CTkFrame(self.loaded_subfrm)
        self.loaded_files_subfrm.grid(row=1, column=0, padx=5, pady=5)

        # Files Header
        ctk.CTkLabel(self.loaded_files_subfrm, text="Files").grid()

        # Files prompts
        ctk.CTkLabel(self.loaded_files_subfrm, text="Network\nVertices\nOther\n.").grid(row=1)

        # loaded vertices subframe
        self.loaded_verts_subfrm = ctk.CTkFrame(self.loaded_subfrm)
        self.loaded_verts_subfrm.grid(row=1, column=1, padx=5, pady=5)

        # Vertices header
        ctk.CTkLabel(self.loaded_verts_subfrm, text="Vertices").grid()

        ctk.CTkLabel(self.loaded_verts_subfrm, text="Max Vertex\nBox Size")

        # Loaded surfaces subframe
        self.loaded_surfs_subfrm = ctk.CTkFrame(self.loaded_subfrm)
        self.loaded_surfs_subfrm.grid(row=1, column=2, padx=5, pady=5)

        # Surfaces header
        ctk.CTkLabel(self.loaded_surfs_subfrm, text="Surfaces").grid(row=0, column=0)
        ctk.CTkLabel(self.loaded_surfs_subfrm, text="Resolution").grid(row=1, column=0)

        # Build sub frame
        self.build_subfrm = ctk.CTkFrame(self.build_frame)
        self.build_subfrm.grid(row=1, columnspan=2)

        # Build sub frame header
        ctk.CTkLabel(self.build_subfrm, text="Build").grid(row=0, column=0)

        # Change Atom radius frame
        self.change_atom_rad_frame = ctk.CTkFrame(self.build_subfrm)
        self.change_atom_rad_frame.grid(row=1, column=0, padx=10, pady=10)

        # Change the atom's radius
        ctk.CTkLabel(self.change_atom_rad_frame, text="Change Radius:").grid(row=0, column=0, padx=2)
        self.cur_elem_rad = ctk.DoubleVar(self, 1.0)
        ctk.CTkEntry(self.change_atom_rad_frame, textvariable=self.cur_elem_rad, width=60).grid(row=0, column=1, pady=2)
        ctk.CTkLabel(self.change_atom_rad_frame, text=u'\u212B', width=25).grid(row=0, column=2, sticky='w')
        self.cur_elem = ctk.StringVar(self)
        ctk.CTkOptionMenu(self.change_atom_rad_frame, variable=self.cur_elem, values=self.elms, width=60).grid(row=0, column=3, padx=2, pady=2)

        ctk.CTkButton(self.change_atom_rad_frame, text="Change", command=self.change_atom_radius).grid(row=0, column=4, padx=2)

        # Estimated time sub frame
        self.est_subfrm = ctk.CTkFrame(self.build_frame)
        self.est_subfrm.grid(row=2, padx=5, pady=5, columnspan=2)

        # Estimated time header
        ctk.CTkLabel(self.est_subfrm, text="Estimated Time").grid()

        # Set surfaces subframe
        self.set_surf_subfrm = ctk.CTkFrame(self.build_frame)
        self.set_surf_subfrm.grid(row=3, column=0, padx=5, pady=5)

        # Set surfaces header
        ctk.CTkLabel(self.set_surf_subfrm, text="Surface Settings").grid(row=0, column=0)

        # Surface Resolution
        ctk.CTkLabel(self.set_surf_subfrm, text="Surface Resolution:").grid(row=2, column=0)
        self.surf_res = ctk.DoubleVar(self, 0.1)
        ctk.CTkEntry(self.set_surf_subfrm, textvariable=self.surf_res, width=60).grid(row=2, column=1)
        ctk.CTkLabel(self.set_surf_subfrm, text=u'\u212B', width=25).grid(row=2, column=2)

        # Flat faces
        self.flat_faces = ctk.BooleanVar(self, False)
        ctk.CTkCheckBox(self.set_surf_subfrm, text="Flat Faces", variable=self.flat_faces).grid(row=3, columnspan=3)

        # Reuse surfaces
        self.reuse_surfs = ctk.BooleanVar(self, False)
        ctk.CTkCheckBox(self.set_surf_subfrm, text="Reuse Surfaces", variable=self.reuse_surfs).grid(row=4, columnspan=3)

        # Set vertices subframe
        self.set_vert_subfrm = ctk.CTkFrame(self.build_frame)
        self.set_vert_subfrm.grid(row=3, column=1, padx=5, pady=5)

        # Set vertices header
        ctk.CTkLabel(self.set_vert_subfrm, text="Vertex Settings").grid(row=0, column=0, columnspan=3)

        # Maximum vertex
        ctk.CTkLabel(self.set_vert_subfrm, text="Max Vert:").grid(row=1, column=0)
        self.max_vert = ctk.DoubleVar(self, 5)
        ctk.CTkEntry(self.set_vert_subfrm, textvariable=self.max_vert, width=60).grid(row=1, column=1)
        ctk.CTkLabel(self.set_vert_subfrm, text=u'\u212B', width=25).grid(row=1, column=2)

        # Maximum box
        ctk.CTkLabel(self.set_vert_subfrm, text="Box Size:").grid(row=2, column=0)
        self.box_size = ctk.DoubleVar(self, 1.5)
        ctk.CTkEntry(self.set_vert_subfrm, textvariable=self.box_size, width=60).grid(row=2, column=1)
        ctk.CTkLabel(self.set_vert_subfrm, text="X", width=25).grid(row=2, column=2)

        # Number of solute layers
        ctk.CTkLabel(self.set_vert_subfrm, text="# SOL Layers:").grid(row=3, column=0)
        self.num_sol_lyrs = ctk.DoubleVar(self, 5)
        ctk.CTkEntry(self.set_vert_subfrm, textvariable=self.num_sol_lyrs, width=60, ).grid(row=3, column=1)

        # Use Loaded verts checkbox
        self.use_loaded_verts = ctk.BooleanVar(self, False)
        ctk.CTkCheckBox(self.set_vert_subfrm, text="Use Loaded Vertices", variable=self.use_loaded_verts).grid(row=4, columnspan=3)

        # Build network button
        ctk.CTkButton(self.build_frame, text="Build", command=self.build_net_button).grid(row=4, columnspan=2)


        ############################################### Save frame #####################################################

        self.save_frame = ctk.CTkFrame(self.view_frame)
        self.save_frame.grid(row=1, columnspan=3, sticky='nsew', padx=20, pady=20)

        # Output directory information frame
        self.out_dir_subfrm = ctk.CTkFrame(self.save_frame)
        self.out_dir_subfrm.grid(row=0, column=0, padx=10, pady=10)

        ctk.CTkLabel(self.out_dir_subfrm, text="Output directory: ").grid(row=0, column=0)
        self.out_dir_str = ctk.StringVar(self.out_dir_subfrm)
        self.out_dir_label = ctk.CTkLabel(self.out_dir_subfrm, textvariable=self.out_dir_subfrm)
        self.out_dir_label.grid(row=0, column=1)
        self.out_dir_button = ctk.CTkButton(self.out_dir_label, text="Browse", command=self.set_out_dir_button)
        self.out_dir_button.grid(row=0, column=2, sticky='w', padx=5, pady=5)

        # Choose specific list subframe
        self.choose_list_subfrm = ctk.CTkFrame(self.save_frame)
        self.choose_list_subfrm.grid(row=1, sticky='nsew', padx=10, pady=10)

        # Choose Header
        ctk.CTkLabel(self.choose_list_subfrm, text="Choose").grid(row=0, columnspan=4)

        # Show radio button
        self.choose_mol_bool = ctk.BooleanVar(self, True)
        self.mol_radio_button = ctk.CTkRadioButton(self.choose_list_subfrm, text="Molecules", command=self.set_show_list_mol, state=self.choose_mol_bool)
        self.mol_radio_button.grid(row=1, column=0, padx=10)
        self.choose_res_bool = ctk.BooleanVar(self, False)
        self.res_radio_button = ctk.CTkRadioButton(self.choose_list_subfrm, text="Residues", command=self.set_show_list_res, state=self.choose_res_bool)
        self.res_radio_button.grid(row=1, column=1, padx=10)
        self.choose_ndx_bool = ctk.BooleanVar(self, False)
        self.ndx_radio_button = ctk.CTkRadioButton(self.choose_list_subfrm, text="Indices", command=self.set_show_list_ndx, state=self.choose_ndx_bool)
        self.ndx_radio_button.grid(row=1, column=2, padx=10)
        self.choose_atom_bool = ctk.BooleanVar(self, False)
        self.atom_radio_button = ctk.CTkRadioButton(self.choose_list_subfrm, text="Atoms", command=self.set_show_list_atom, state=self.choose_atom_bool)
        self.atom_radio_button.grid(row=1, column=3, padx=10)

        # Choose lists
        self.current_selection = ctk.StringVar(self)
        self.current_selection_atoms = []
        self.mol_list, self.res_list, self.atom_list, self.ndx_list = [], [], [], []
        self.choose_mol_list = ctk.CTkOptionMenu(master=self.choose_list_subfrm, variable=self.current_selection, values=self.mol_list, width=400)
        self.choose_mol_list.grid(row=2, columnspan=4, sticky='ew', padx=10, pady=6)
        self.choose_res_list = ctk.CTkOptionMenu(master=self.choose_list_subfrm, variable=self.current_selection, values=self.res_list, width=400)
        self.choose_res_list.grid(row=2, columnspan=4, sticky='ew', padx=10, pady=6)
        self.choose_ndx_list = ctk.CTkOptionMenu(master=self.choose_list_subfrm, variable=self.current_selection, values=self.ndx_list, width=400)
        self.choose_ndx_list.grid(row=2, columnspan=4, sticky='ew', padx=10, pady=6)
        self.choose_atom_list = ctk.CTkOptionMenu(master=self.choose_list_subfrm, variable=self.current_selection, values=self.atom_list, width=400)
        self.choose_atom_list.grid(row=2, columnspan=4, sticky='ew', padx=10, pady=6)

        # Group frame
        self.group_frame = ctk.CTkFrame(self.save_frame)
        self.group_frame.grid(row=2, padx=10, pady=10)

        # Group 1 sub frame
        self.g1_subfrm = ctk.CTkFrame(self.group_frame)
        self.g1_subfrm.grid(row=0, column=0, padx=10, pady=10)

        # Group 1 Header
        ctk.CTkLabel(self.g1_subfrm, text="Group 1").grid(row=0, columnspan=3)

        # Group 1 Buttons
        ctk.CTkButton(self.g1_subfrm, text="Add", command=self.add_g1_button, width=40).grid(row=1, column=0, sticky='nsew', padx=2)
        ctk.CTkButton(self.g1_subfrm, text="Undo", command=self.undo_g1_button, width=30).grid(row=1, column=1, sticky='nsew', padx=2)
        ctk.CTkButton(self.g1_subfrm, text="Reset", command=self.reset_g1_button, width=30).grid(row=1, column=2, sticky='nsew', padx=2)

        # Group 1 selection list
        self.g1_sele_str = ctk.StringVar(self, "\n\n\n\n")
        ctk.CTkLabel(self.g1_subfrm, textvariable=self.g1_sele_str, width=200).grid(row=2, columnspan=3, sticky='nsew')

        # Group 1 export button
        ctk.CTkButton(self.group_frame, text="Export\nGroup", command=self.export_g1_button).grid(row=1, column=0)

        # Seperator
        ttk.Separator(self.group_frame, orient=tk.VERTICAL).grid(row=0, column=1, rowspan=1, sticky='ns')

        # Group 2 sub frame
        self.g2_subfrm = ctk.CTkFrame(self.group_frame)
        self.g2_subfrm.grid(row=0, column=2, padx=10, pady=10)

        # Group 2 Header
        ctk.CTkLabel(self.g2_subfrm, text="Group 2").grid(row=0, columnspan=3)

        # Group 2 buttons
        ctk.CTkButton(self.g2_subfrm, text="Add", command=self.add_g2_button, width=40).grid(row=1, column=0, sticky='nsew', padx=2)
        ctk.CTkButton(self.g2_subfrm, text="Undo", command=self.undo_g2_button, width=30).grid(row=1, column=1, sticky='nsew', padx=2)
        ctk.CTkButton(self.g2_subfrm, text="Reset", command=self.reset_g2_button, width=30).grid(row=1, column=2, sticky='nsew', padx=2)

        # Group 2 selection list
        self.g2_sele_str = ctk.StringVar(self, "\n\n\n\n")
        ctk.CTkLabel(self.g2_subfrm, textvariable=self.g2_sele_str, width=200).grid(row=2, columnspan=3, sticky='nsew')

        # Group 2 export button
        ctk.CTkButton(self.group_frame, text="Export\nGroup", command=self.export_g2_button).grid(row=1, column=2)

        # Export interface button
        ctk.CTkButton(self.group_frame, text="Export\nInterface", command=self.export_iface_button).grid(row=2, columnspan=3, padx=5, pady=5)

        ######################################### Top buttons ##########################################################

        self.load_frame_button = ctk.CTkButton(self.view_frame, text="Load", command=self.load_frame.tkraise)
        self.load_frame_button.grid(row=0, sticky='nse', pady=10)
        self.build_frame_button = ctk.CTkButton(self.view_frame, text="Build", command=self.build_frame.tkraise)
        self.build_frame_button.grid(row=0, column=1, sticky='news', pady=10)
        self.save_frame_button = ctk.CTkButton(self.view_frame, text="Save", command=self.save_frame.tkraise)
        self.save_frame_button.grid(row=0, column=2, sticky='nsw', pady=10)
        self.load_frame.tkraise()

    ############################################## Functions ###########################################################

    ############################################## Load Frame ##########################################################

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
        self.sys.load_sys_file(file_path)
        # Set the system information
        # We want to get the number of atoms, the number of molecules, etc
        myStr = str(len(self.sys.atoms)) + '\n' + str(len(self.sys.mols)) + '\n' + str(len(self.sys.residues)) + \
                "\n   ~   \n   ~   \n   ~   "
        # Set the variables
        self.sys_data.set(myStr)
        # Set the molecule names
        self.mol_list = self.sys.mol_names
        self.res_list = self.sys.res_names
        self.ndx_list = self.sys.ndx_names
        self.atom_list = self.sys.atom_names


        self.choose_mol_list = ctk.CTkOptionMenu(self.choose_list_subfrm, values=self.sys.mol_names, variable=self.current_selection, width=400)
        self.choose_mol_list.grid(row=2, columnspan=4, sticky='ew', padx=10, pady=6)
        self.choose_res_list = ctk.CTkOptionMenu(self.choose_list_subfrm, values=self.sys.res_names, variable=self.current_selection, width=400)
        self.choose_res_list.grid(row=2, columnspan=4, sticky='ew', padx=10, pady=6)
        self.choose_ndx_list = ctk.CTkOptionMenu(self.choose_list_subfrm, values=self.sys.ndx_names, variable=self.current_selection, width=400)
        self.choose_ndx_list.grid(row=2, columnspan=4, sticky='ew', padx=10, pady=6)
        self.choose_atom_list = ctk.CTkOptionMenu(self.choose_list_subfrm, values=self.sys.atom_names, variable=self.current_selection, width=400)
        self.choose_atom_list.grid(row=2, columnspan=4, sticky='ew', padx=10, pady=6)
        # Set the output directory
        self.out_dir_str.set(set_output_dir(self.sys))
        # Raise the buttons
        self.load_frame_button.tkraise()
        self.build_frame_button.tkraise()
        self.save_frame_button.tkraise()

    def load_ndx_button(self):
        # File grabber pop up
        file_path = filedialog.askopenfilename()
        # Set the system's index file
        self.sys.ndx_file = file_path
        # Load the indices into the system
        self.sys.load_ndx()

    # Load network button function. Pulls up the file browser and lets the user select their vorpy saved system
    def load_net_button(self):
        # Check to see if a network exists already
        if self.sys.net is None:
            self.sys.net = Network(self, self.sys.atoms)
        # File grabber pop up
        file_path = filedialog.askopenfilename()
        self.net_file = file_path
        # Check to see if there is a system file
        if len(self.sys.atoms) < 1:
            return
        else:
            self.sys.load_net(self.net_file)
        self.sys.set_output_directory()
        self.net_data.set(str(len(self.sys.net.surfs)) + "\n" + str(len(self.sys.net.edges)) + "\n" +
                          str(len(self.sys.net.verts)) + "\n~\n~\n~")

    def load_verts_button(self):
        # Get the filepath for the vertices
        # File grabber pop up
        file_path = filedialog.askopenfilename()
        self.sys.load_verts(vert_file=file_path)
        self.net_data.set(str(len(self.sys.net.surfs)) + "\n" + str(len(self.sys.net.edges)) + "\n" +
                          str(len(self.sys.net.verts)) + "\n~\n~\n~")

    ################################################# Build Frame ######################################################

    def change_atom_radius(self):
        # When pressed, the current atom selection's radius changes
        old_rad = self.sys.radii[1][self.sys.radii[0].index(self.cur_elem.get())]
        self.sys.radii[1][self.sys.radii[0].index(self.cur_elem.get())] = self.cur_elem_rad.get()
        # print an update
        my_elem = self.cur_elem.get() + " "
        my_elem = my_elem[0].upper() + my_elem[1:]
        print("\r{} radius changed from {} to {}".format(my_elem, old_rad,
                                                         self.cur_elem_rad.get()), end="")

    def build_net_button(self):
        if self.sys is None:
            return
        print(self.flat_faces.get())

        self.sys.build_network(max_vert=self.max_vert.get(), surf_res=self.surf_res.get(),
                               flat_faces=self.flat_faces.get(), box_size=self.box_size.get(),
                               find_verts=self.use_loaded_verts.get())
        self.net_data.set(str(len(self.sys.net.surfs)) + "\n" + str(len(self.sys.net.edges)) + "\n" +
                          str(len(self.sys.net.verts)) + "\n    ~    \n    ~    \n    ~   ")

        self.build_data.set(str(self.sys.net.my_time) + "\n" + str(self.sys.net.cpu_time) + "\n" + str(self.sys.net.min_dist) +
                            "\n" + str(self.sys.net.box_size) + "\n" + str(self.sys.net.max_vert))


    ############################################### Export Frame #######################################################

    def set_out_dir_button(self):
        """
        Sets the output directory for the export files and folders
        :return:
        """
        # File grabber pop up
        file_path = filedialog.askdirectory()
        # Create the System
        if file_path:
            # Set the system's output directory and the string version of the file location
            self.sys.output_directory = file_path
            self.out_dir_str.set(file_path[:12] + ' ... ' + file_path[-12:])

    def set_show_list_mol(self):
        """
        Shows the molecule list in the dropdown menu of the choose section of the save frame
        :return:
        """
        self.choose_mol_list.tkraise()
        self.choose_mol_bool.set(True)
        self.choose_ndx_bool.set(False)
        self.choose_res_bool.set(False)
        self.choose_atom_bool.set(False)

    def set_show_list_res(self):
        """
        Shows the reside list in the dropdown menu of the choose section of the save frame
        :return:
        """
        self.choose_res_list.tkraise()
        self.choose_mol_bool.set(False)
        self.choose_ndx_bool.set(False)
        self.choose_res_bool.set(True)
        self.choose_atom_bool.set(False)

    def set_show_list_ndx(self):
        """
        Shows the indices list in the dropdown menu of the choose section of the save frame
        :return:
        """
        self.choose_ndx_list.tkraise()
        self.choose_mol_bool.set(False)
        self.choose_ndx_bool.set(True)
        self.choose_res_bool.set(False)
        self.choose_atom_bool.set(False)

    def set_show_list_atom(self):
        """
        Shows the atoms list in the dropdown menu of the choose section of the save frame
        :return:
        """
        self.choose_atom_list.tkraise()
        self.choose_mol_bool.set(False)
        self.choose_ndx_bool.set(False)
        self.choose_res_bool.set(False)
        self.choose_atom_bool.set(True)

    def get_current_selection_atoms(self):
        """
        Finds and sets the current selection atoms
        :return:
        """
        if self.choose_mol_bool.get():
            self.current_selection_atoms = self.sys.mols[self.sys.mol_names.index(self.current_selection.get())]
        elif self.choose_res_bool.get():
            self.current_selection_atoms = self.sys.residues[self.sys.res_names.index(self.current_selection.get())]
        elif self.choose_ndx_bool.get():
            self.current_selection_atoms = self.sys.ndxs[self.sys.ndx_names.index(self.current_selection.get())]
        elif self.choose_atom_bool.get():
            self.current_selection_atoms = self.sys.atoms[self.sys.atom_names.index(self.current_selection.get())]


    def add_g1_button(self):
        """
        Adds information to the group 1 group and analyzes the new information
        :return:
        """
        # Get the current selection atoms for reference
        self.get_current_selection_atoms()
        # If group 1 has not been created, make a class and set the string with the current selection
        if self.g1 is None:
            self.g1 = Group(self.sys.net, self.current_selection_atoms)
            self.g1.select_strs = [self.current_selection.get()]
        # If the group exists, use the add selection method with the current selection
        else:
            self.g1.add_sele(self.current_selection_atoms, self.current_selection.get())
        # Make an interface if both groups are populated
        if self.g2 is not None:
            self.g1.bff, self.g2.bff = self.g2, self.g1
            self.g2.get_info()
        # Get the information about both the cell and the interface
        self.g1.get_info()
        # Set the selections string
        self.g1_sele_str.set("\n".join(self.g1.select_strs))

    def add_g2_button(self):
        """
        Adds information to the group 2 group and analyzes the new information
        :return:
        """
        # Get the current selection atoms for reference
        self.get_current_selection_atoms()
        # If group 2 has not been created, make a class and set the string with the current selection
        if self.g2 is None:
            self.g2 = Group(self.sys.net, self.current_selection_atoms)
            self.g2.select_strs = [self.current_selection.get()]
        # If the group exists, use the add selection method with the current selection
        else:
            self.g2.add_sele(self.current_selection_atoms, self.current_selection.get())
        # Make an interface if both groups are populated
        if self.g1 is not None:
            self.g2.bff, self.g1.bff = self.g1, self.g2
            self.g1.get_info()
        # Get the information about both the cell and the interface
        self.g2.get_info()
        # Set the selections string
        self.g2_sele_str.set("\n".join(self.g2.select_strs))

    def undo_g1_button(self):
        self.g1.undo_sele()
        self.g1_sele_str.set("\n".join(self.g1.select_strs))

    def undo_g2_button(self):
        self.g2.undo_sele()
        self.g2_sele_str.set("\n".join(self.g2.select_strs))

    def reset_g1_button(self):
        self.g1 = Group(self.sys.net, [])
        self.g1_sele_str.set("")

    def reset_g2_button(self):
        self.g2 = Group(self.sys.net, [])
        self.g2_sele_str.set("")

    def export_g1_button(self):
        export_body(self.g1, info_file=True, outer_atoms=False)

    def export_g2_button(self):
        export_body(self.g2, info_file=True, outer_atoms=False)

    def export_iface_button(self):
        export_iface([self.g1, self.g2], info_file=True, interface_atoms=False)
