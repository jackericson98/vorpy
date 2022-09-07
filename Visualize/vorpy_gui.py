# Outside Imports
import os
import tkinter as tk
from tkinter import filedialog
# Internal Imports
from System.system import System


class Vorpy:
    """Vorpy GUI class. When instantiated the Gui will launch"""
    def __init__(self, width=450, height=600):
        # Set up the window
        self.vp_main = tk.Tk()
        self.width = width
        self.height = height
        self.vp_main.geometry(str(width) + "x" + str(height))
        self.vp_main.title('vorpy')
        self.file = ""
        self.sys_box_size = 2
        self.sys_resolution = 0.1
        self.vorpy_directory = os.getcwd()
        self.output_directory = os.getcwd()
        self.verts_loaded = False

        # Instantiate the System
        self.sys = None

        # Set up the strings
        self.file_name = tk.StringVar(self.vp_main, "No System Selected")
        self.dd_name = tk.StringVar(self.vp_main)
        self.dd_var = tk.StringVar(self.vp_main)
        self.box_size = tk.StringVar(self.vp_main, str(self.sys_box_size))
        self.resolution = tk.StringVar(self.vp_main, str(self.sys_resolution))
        self.vta = tk.BooleanVar(self.vp_main)
        self.parallelize = tk.BooleanVar(self.vp_main)
        self.output_dir_str = tk.StringVar(self.vp_main, os.getcwd()[:12] + ' ... ' + os.getcwd()[-12:])

        self.output_all = tk.BooleanVar(self.vp_main, True)
        self.output_sys = tk.BooleanVar(self.vp_main, False)
        self.output_atoms = tk.BooleanVar(self.vp_main, False)
        self.output_mols = tk.BooleanVar(self.vp_main, False)
        self.output_surfs = tk.BooleanVar(self.vp_main, False)
        self.output_analysis = tk.BooleanVar(self.vp_main, False)
        self.output_verts = tk.BooleanVar(self.vp_main, False)
        self.output_pdb = tk.BooleanVar(self.vp_main, False)

        # Set up the main frame
        f = tk.Frame(self.vp_main)
        f.place(relx=0.5, rely=0.5, anchor='c')
        # Header:
        tk.Label(f, text="VorPy", font=('Times New Roman bold', 40)).grid(row=0, column=0, columnspan=3)
        tk.Label(f, text="").grid(row=1)
        # Input info
        tk.Label(f, text="Load", font=('Times New Roman bold', 20)).grid(row=2, column=0, columnspan=2, sticky='w')
        tk.Label(f, text="System: ", font=('Times New Roman', 15)).grid(row=3, column=0, sticky='w')
        tk.Label(f, textvariable=self.file_name, font=('Times New Roman', 15))\
            .grid(row=3, column=1, columnspan=2, sticky='w')
        tk.Button(f, text="Load Vertices", command=self.add_vertices).grid(row=4, column=0, sticky='e')
        tk.Button(f, text="Build System", command=self.build_sys_button).grid(row=4, column=1)
        tk.Button(f, text="Load System", command=self.load_sys_button).grid(row=4, column=2, sticky='w')
        tk.Label(f, text="").grid(row=5)

        # Settings
        tk.Label(f, text='Settings', font=('Times New Roman bold', 20)).grid(row=6, column=0, columnspan=2, sticky='w')
        tk.Label(f, text="Container Size: ").grid(row=7, column=0, sticky='w')
        tk.Entry(f, textvariable=self.box_size).grid(row=7, column=1, sticky='w')
        tk.Label(f, text="Resolution: ").grid(row=8, column=0, sticky='w')
        tk.Entry(f, textvariable=self.resolution).grid(row=8, column=1, sticky='w')
        tk.Checkbutton(f, text="Voronota System ", variable=self.vta, onvalue=True, offvalue=False)\
            .grid(row=9, column=0, sticky='w')
        tk.Checkbutton(f, text="Parallelize ", variable=self.parallelize, onvalue=True, offvalue=False)\
            .grid(row=10, column=0, sticky='w')
        tk.Label(f, text="").grid(row=11)

        # Outputs
        tk.Label(f, text="Outputs", font=("Times New Roman bold", 20)).grid(row=12, column=0, columnspan=2, sticky='w')
        tk.Label(f, text="Output Directory: ").grid(row=13, column=0, sticky='w')
        tk.Label(f, textvariable=self.output_dir_str).grid(row=13, column=1, sticky='w')
        tk.Button(f, text="Change", command=self.change_output_directory).grid(row=13, column=2, sticky='e')

        tk.Label(f, text="Output Files").grid(row=14, column=0, sticky='w')
        tk.Checkbutton(f, text="All", variable=self.output_all, onvalue=True, offvalue=False)\
            .grid(row=15, column=0, sticky='w')
        tk.Checkbutton(f, text="System", variable=self.output_sys, onvalue=True, offvalue=False)\
            .grid(row=15, column=1, sticky='w')
        tk.Checkbutton(f, text="System pdb ", variable=self.output_pdb, onvalue=True, offvalue=False) \
            .grid(row=15, column=2, sticky='w')
        tk.Checkbutton(f, text="Atom Cells", variable=self.output_atoms, onvalue=True, offvalue=False)\
            .grid(row=16, column=0, sticky='w')
        tk.Checkbutton(f, text="Surfaces", variable=self.output_surfs, onvalue=True, offvalue=False) \
            .grid(row=16, column=1, sticky='w')
        tk.Checkbutton(f, text="Molecule Interfaces", variable=self.output_mols, onvalue=True, offvalue=False)\
            .grid(row=16, column=2, sticky='w')
        tk.Checkbutton(f, text="Analysis ", variable=self.output_analysis, onvalue=True, offvalue=False)\
            .grid(row=17, column=0, sticky='w')
        tk.Checkbutton(f, text="Vertices ", variable=self.output_verts, onvalue=True, offvalue=False)\
            .grid(row=17, column=1, sticky='w')
        tk.Label(f, text="").grid(row=18)

        # Run program
        tk.Button(f, text="Build System", font=("Times New Roman bold", 20), command=self.build_network_button) \
            .grid(row=19, column=0, columnspan=3)

        # End the loop
        self.vp_main.mainloop()

    def load_sys_button(self):
        # File grabber pop up
        file_path = filedialog.askopenfilename()
        # Create the System
        if file_path:
            self.file = file_path
            self.sys = self.sys = System(file_path)
        filename = ""
        i = -1
        while self.file[i] != "/":
            filename = filename + self.file[i]
            i -= 1

        self.file_name.set(filename[::-1][:-4])

    def build_sys_button(self):
        pass

    def add_vertices(self):
        if self.sys is None:
            ErrorBox("Please select a system to add the vertices to")
        # File grabber pop up
        file_path = filedialog.askopenfilename()
        self.sys.add_verts(file_path)
        self.verts_loaded = True

    def change_output_directory(self):
        # File grabber pop up
        file_path = filedialog.askdirectory()
        # Create the System
        if file_path:
            self.output_directory = file_path
            self.output_dir_str.set(os.getcwd()[:12] + ' ... ' + os.getcwd()[-12:])
        else:
            ErrorBox("Directory not changed")

    # Build network button function.
    def build_network_button(self):
        # Create the System
        if self.sys is None:
            ErrorBox("Please select a system")
        # Create a new directory of one has not been indicated
        if self.output_directory[-len(self.sys.name):] != self.sys.name:
            # If the system doesn't have a name
            if self.sys.name == '':
                self.sys.name = "User_Data"
            i = 0
            while True:
                try:
                    i_str = str(i)
                    if i == 0:
                        i_str = ""
                    os.mkdir(os.getcwd() + "/" + self.sys.name + i_str)
                    break
                except FileExistsError:
                    i += 1
            self.output_directory = os.getcwd() + "/" + self.sys.name + i_str
        # Build the network
        self.sys.build_network(get_verts=not self.verts_loaded, export_verts=(self.output_verts.get() or self.output_all.get()),
                               directory=self.output_directory, box_size=float(self.box_size.get()),
                               min_dist=float(self.resolution.get()))
        # Analyze the network
        if self.output_analysis.get() or self.output_all.get():
            self.sys.analyze()

        # Export the system
        self.sys.export(directory=self.output_directory, export_all=self.output_all.get(),
                        export_sys=self.output_sys.get(), export_mols=self.output_mols.get(),
                        export_atoms=self.output_atoms.get(), export_analysis=self.output_analysis.get(),
                        export_surfs=self.output_surfs.get())
        os.chdir(self.vorpy_directory)


class ErrorBox:
    """Error box class. Used to indicate if an error has occurred and takes the message as an input"""
    def __init__(self, error_message):
        self.eroot = tk.Tk()
        self.eroot.geometry("300x75")

        self.error_message = tk.StringVar(self.eroot, error_message)
        self.eroot.title("vorpy")
        error_frame = tk.Frame(self.eroot)
        error_frame.place(relx=0.5, rely=0.45, anchor='c')

        tk.Label(error_frame, text="Error!", font=("Times New Roman bold", 15)).grid(row=0, column=0, sticky='w')
        tk.Label(error_frame, textvariable=self.error_message).grid(row=1, column=0, sticky='e')
        tk.Button(error_frame, text="OK", command=self.eroot.quit()).grid(row=2, column=0, sticky='e')

        # End the loop
        self.eroot.mainloop()


if __name__ == '__main__':
    os.chdir("..")
    Vorpy()
