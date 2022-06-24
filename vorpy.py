# Outside Imports
import tkinter as tk
from tkinter import filedialog, Button, CENTER
# Internal Imports
from build_network import build_network
from build_mesh import build_meshes
from visualize import *
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)


class Vorpy:
    """Vorpy GUI class. When instantiated the Gui will launch"""
    def __init__(self, width=1000, height=800):
        # Set up the window
        self.vp_main = tk.Tk()
        self.width = width
        self.height = height
        self.vp_main.geometry(str(width) + "x" + str(height))
        self.vp_main.title('vorpy')

        # Instantiate the system
        self.sys = None

        # Set up the strings
        self.name = tk.StringVar()
        self.name.set("No File Selected")

        # Text boxes:
        self.head = tk.Label(text="VorPy", font=('Helvetica bold', 40))
        self.filename = tk.Label(self.vp_main, textvariable=self.name, font=('Times New Roman', 20))

        # Plot
        self.fig = Figure(figsize=(self.width / 100, self.height / 100), dpi=50)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.vp_main)

        # Buttons:
        self.get_file = tk.Button(text="Load Molecule", command=self.load_molecule_button)
        self.make_network = tk.Button(text="Build Network", command=self.build_network_button)
        self.make_meshes = tk.Button(text="Build Meshes", command=self.build_meshes_button)
        self.exit = tk.Button(text="Exit", command=self.vp_main.destroy)
        self.atoms_butt = None
        self.verts_butt = None

        # Place the items
        self.head.place(x=width/2 - width/12, y=height/16)
        self.filename.place(x=width*19/32, y=7/32*height)
        self.get_file.place(x=3/16*width, y=3/8*height)
        self.make_network.place(x=3/16*width, y=1/2*height)
        self.make_meshes.place(x=3/16*width, y=5/8*height)
        self.atoms_butt = tk.Button(text="Atoms", command=self.show_atoms)
        self.verts_butt = tk.Button(text="Vertices", command=self.show_verts)
        self.atoms_butt.place(x=3 / 4 * self.width - 3 / 16 * self.width, y=3 / 4 * self.height + 1 / 16 * self.height)
        self.verts_butt.place(x=3 / 4 * self.width + 1 / 16 * self.width, y=3 / 4 * self.height + 1 / 16 * self.height)

        # End the loop
        self.vp_main.mainloop()

    def load_molecule_button(self):
        # File grabber pop up
        file_path = filedialog.askopenfilename()
        # Create the system
        self.sys = System(file_path)

        # Set the name in reverse order since the letters were added backwards
        self.name.set(self.sys.file_name)

    # Build network button function.
    def build_network_button(self):
        # Build the network
        build_network(self.sys)  # Try statment with voronota_verts

        # Print out all the vertices
        for i in range(len(self.sys.net.verts)):
            print(self.sys.net.verts[i])

    def build_meshes_button(self):
        # Build the meshes
        build_meshes(self.sys)
        # Plot the surfaces
        plot_surfs(surfs=self.sys.net.surfs)

    def show_atoms(self):
        if self.sys:
            plot_atoms(self.sys.atoms, fig=self.fig)
        # Create Canvas
        self.canvas.draw()

        self.canvas.get_tk_widget().pack()

        # creating the Matplotlib toolbar
        # toolbar = NavigationToolbar2Tk(self.canvas, self.vp_main)
        #
        # toolbar.update()

        # placing the toolbar on the Tkinter window
        self.canvas.get_tk_widget().place(x=self.width * 7 / 16, y=self.height * 5 / 16)

    def show_verts(self):

        plot_verts(self.sys.net.verts, fig=self.fig)
        # Create Canvas
        self.canvas.draw()

        self.canvas.get_tk_widget().pack()


class ErrorBox:
    """Error box class. Used to indicate if an error has occurred and takes the message as an input"""
    def __init__(self, error_message):
        self.eroot = tk.Tk()

        self.eroot.title("vorpy")
        self.error = tk.Label(text=error_message, font="none 14 bold")

        self.error.config(anchor=CENTER)
        self.error.pack(padx=20, pady=20)
        # End the loop
        self.eroot.mainloop()


Vorpy()
