# Outside Imports
import tkinter as tk
from tkinter import filedialog, Button, CENTER
# Internal Imports
from file import get_data
from build_network import build_network
from build_mesh import build_meshes
from visualize import *
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)


class Vorpy:
    """Vorpy GUI class. When instantiated the Gui will launch"""
    def __init__(self, width=1000, height=800):
        # Set up the window
        self.window = tk.Tk()
        self.width = width
        self.height = height
        self.window.geometry(str(width) + "x" + str(height))
        self.window.title('vorpy')

        # Instantiate the system
        self.sys = None

        # Set up the strings
        self.name = tk.StringVar()
        self.name.set("No File Selected")

        # Text boxes:
        self.head = tk.Label(text="VorPy", font=('Helvetica bold', 40))
        self.filename = tk.Label(self.window, textvariable=self.name, font=('Times New Roman', 20))

        # Buttons:
        self.get_file = tk.Button(text="Load Molecule", command=self.load_molecule_button)
        self.make_network = tk.Button(text="Build Network", command=self.build_network_button)
        self.make_meshes = tk.Button(text="Build Meshes", command=self.build_meshes_button)
        self.exit = tk.Button(text="Exit", command=self.window.destroy)
        self.plot_button = tk.Button(text="Plot", command=self.plot)
        self.atoms_butt = None
        self.verts_butt = None

        # Place the items
        self.head.place(x=width/2 - width/12, y=height/16)
        self.filename.place(x=width*19/32, y=7/32*height)
        self.get_file.place(x=3/16*width, y=3/8*height)
        self.make_network.place(x=3/16*width, y=1/2*height)
        self.make_meshes.place(x=3/16*width, y=5/8*height)
        self.plot_button.place(x=3/4*width - 1/8*width, y=3/4*height + 1/8*height)

        # End the loop
        self.window.mainloop()

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
        build_network(self.sys)
        # Print out all the vertices
        for i in range(len(self.sys.net.verts)):
            print(self.sys.net.verts[i])

    def build_meshes_button(self):
        # Build the meshes
        build_meshes(self.sys)
        # Plot the surfaces
        plot_surfs(surfs=self.sys.net.surfs)

    def show_atoms(self, fig, canvas):

        plot_atoms(self.sys.atoms, fig=fig)
        # Create Canvas
        canvas.draw()

        # canvas.get_tk_widget().pack()

        # creating the Matplotlib toolbar
        toolbar = NavigationToolbar2Tk(canvas, self.window)

        toolbar.update()

        # placing the toolbar on the Tkinter window
        canvas.get_tk_widget().place(x=self.width * 7 / 16, y=self.height * 5 / 16)

    def show_verts(self, fig, canvas):

        canvas.flush_events()
        plot_atoms(self.sys.atoms, fig=fig)
        plot_verts(self.sys.net.verts, fig=fig)
        # Create Canvas
        canvas.draw()

        # canvas.get_tk_widget().pack()

    def plot(self):
        # Create figure
        fig = Figure(figsize=(self.width / 100, self.height / 100), dpi=50)
        canvas = FigureCanvasTkAgg(fig, master=self.window)
        # Use the plot atoms method to create plot contents
        self.atoms_butt = tk.Button(text="Atoms", command=self.show_atoms(fig, canvas))
        self.verts_butt = tk.Button(text="Vertices", command=self.show_verts(fig, canvas))

        self.atoms_butt.place(x=3 / 4 * self.width - 3 / 16 * self.width, y=3 / 4 * self.height + 1 / 16 * self.height)
        self.verts_butt.place(x=3 / 4 * self.width + 1 / 16 * self.width, y=3 / 4 * self.height + 1 / 16 * self.height)

        # creating the Matplotlib toolbar
        toolbar = NavigationToolbar2Tk(canvas, self.window)

        toolbar.update()

        # placing the toolbar on the Tkinter window
        canvas.get_tk_widget().place(x=self.width * 7 / 16, y=self.height * 5 / 16)
        canvas.flush_events()


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
