# Outside Imports
import tkinter as tk
from tkinter import filedialog, Button
# Internal Imports
from load_system import read_pdb
from build_network import build_network
from build_mesh import build_meshes
from visualize import *


class Vorpy:
    """Vorpy GUI class. When instantiated the Gui will launch"""
    def __init__(self, mySys=None):
        # Set up the window
        self.root = tk.Tk()
        self.root.title('vorpy')
        # Instantiate the system
        self.sys = mySys
        # Set up the strings
        self.name = tk.StringVar()
        self.name.set("No File Selected")
        # Set up the labels
        self.head = tk.Label(text="VorPy", font=('Helvetica bold', 40))
        self.filename = tk.Label(self.root, textvariable=self.name)
        # Set up the buttons
        self.get_file = tk.Button(text="Load Molecule", command=self.load_molecule_button)
        self.make_network = tk.Button(text="Build Network?", command=self.build_network_button)
        self.make_meshes = tk.Button(text="Build Meshes?", command=self.build_meshes_button)
        self.exit = tk.Button(text="Exit", command=self.root.destroy)
        # Place the items
        self.head.grid(row=0)
        self.filename.grid(row=2)
        self.get_file.grid(row=3, column=0)
        self.make_network.grid(row=3, column=2)
        self.make_meshes.grid(row=3, column=4)
        self.exit.grid(row=5)
        # End the loop
        self.root.mainloop()

    def load_molecule_button(self):
        # File grabber pop up
        file_path = filedialog.askopenfilename()
        # Create the system
        self.sys = read_pdb(file_path)
        # Grab the name of the file from the input file
        if self.sys.info['header'] == '':
            self.name.set(self.sys.info["header"][0][0].capitalize())
        # If the input file did not have a name grab the name of the file itself
        else:
            filename = ""
            i = -1
            # Go through each char in the path from the back and stop at the first slash
            while file_path[i] != "/":
                # When the first period is encountered reset the file name
                if file_path[i] == '.':
                    filename = ''
                else:
                    filename = filename + file_path[i]
                i -= 1
            # Set the name in reverse order since the letters were added backwards
            self.name.set(''.join(reversed(filename)))

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


# Driver code
Vorpy()
