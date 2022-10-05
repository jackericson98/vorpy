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
