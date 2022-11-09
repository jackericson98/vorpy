import os

import customtkinter as ctk
import tkinter as tk
from tkinter import filedialog
from System.system import System


class Vorpy(ctk.CTk):
    def __init__(self):
        # Get all the characteristics from the CTk class
        super().__init__()

        self.sys = None
        self.sys_file = None
        self.net_file = None

        # Set the geometry
        self.geometry("800x500")
        self.title("Vorpy")
        self.minsize(400, 300)
        self.main_frame = ctk.CTkFrame(self)
        self.main_frame.grid(sticky='nsew')
        self.main_frame.columnconfigure(0, weight=1)
        self.main_frame.columnconfigure(1, weight=1)

        self.sys_frame = ctk.CTkFrame(self.main_frame)
        self.sys_frame.grid(row=0, column=0, padx=20, pady=20, sticky='nsew', columnspan=2)
        self.sys_frame.rowconfigure(0, weight=1)
        self.sys_frame.columnconfigure(0, weight=1)

        # System
        self.sys_name = tk.StringVar(self, "No system selected")
        self.sys_label = ctk.CTkLabel(self.sys_frame, textvariable=self.sys_name)
        self.sys_label.grid(row=0, column=0, padx=10, pady=10)
        self.load_sys_button = ctk.CTkButton(self.sys_frame, text="Load System", command=self.load_sys_button)
        self.load_sys_button.grid(row=1, column=0, padx=10, pady=10)

        # Settings frame
        self.settings_frame = ctk.CTkFrame(self.main_frame)
        self.settings_frame.grid(row=1, column=0, sticky='nsew', padx=20, pady=20)

        # Settings
        self.res_val = tk.DoubleVar(self, 0.1)
        self.res_val_entry = ctk.CTkEntry(self.settings_frame, textvariable=self.res_val, width=20)
        self.res_val_entry.grid(row=0, column=0, padx=10, pady=10, sticky='ew')
        ctk.CTkLabel(self.settings_frame, text=u'\u212B').grid(row=0, column=1, sticky='w')

        self.max_vert = tk.DoubleVar(self, 5.0)
        self.max_vert_entry = ctk.CTkEntry(self.settings_frame, textvariable=self.max_vert)
        self.max_vert_entry.grid(row=1, column=0, padx=10, pady=10, sticky='ew')
        ctk.CTkLabel(self.settings_frame, text=u'\u212B').grid(row=1, column=1, sticky='w')

        self.box_size = tk.DoubleVar(self, 1.3)
        self.box_size_entry = ctk.CTkEntry(self.settings_frame, textvariable=self.box_size)
        self.box_size_entry.grid(row=2, column=0, padx=10, pady=10, sticky='ew')
        ctk.CTkLabel(self.settings_frame, text='x').grid(row=2, column=1, sticky='w')

        self.find_sol = tk.BooleanVar(self, True)
        ctk.CTkCheckBox(self.settings_frame, text="Find Solute Vertices", variable=self.find_sol).grid(row=0, column=2)


        self.net_frame = ctk.CTkFrame(self.main_frame)
        self.net_frame.grid(row=1, column=1, padx=20, pady=20, sticky='nsew')

        # Network
        self.net_name = tk.StringVar(self, "No network selected")
        self.net_label = ctk.CTkLabel(self.net_frame, textvariable=self.net_name)
        self.net_label.grid(row=0, column=0, padx=10, pady=10)
        self.load_net_button = ctk.CTkButton(self.net_frame, text="Load Network", command=self.load_net_button)
        self.load_net_button.grid(row=1, column=0, padx=10, pady=10)


        # Run frame
        self.run_frame = ctk.CTkFrame(self.main_frame)
        self.run_frame.grid(row=2, column=0, columnspan=2, sticky='nsew', padx=20, pady=20)
        self.run_button = ctk.CTkButton(self.run_frame, text="Run", command=self.build_net_button)
        self.run_button.pack()

    # Load system button function. Calls the file browser and sets the system
    def load_sys_button(self):
        # Reset the system
        self.sys = System()
        # Re-connect the gui and the system
        self.sys.gui = self
        # File grabber pop up
        file_path = filedialog.askopenfilename()
        # Set the file path
        if file_path is not None:
            self.sys_file = file_path
        else:
            return
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

    # Load network button function. Pulls up the file browser and lets the user select their vorpy saved system
    def load_net_button(self):
        # File grabber pop up
        file_path = filedialog.askopenfilename()
        self.net_file = file_path
        # Check to see if there is a system file
        if len(self.sys.atoms) < 1:
            return
        else:
            self.sys.load_net(self.net_file)
        self.sys.set_output_directory()
        self.net_name.set(self.sys.net.name)

    def build_net_button(self):
        if self.sys is None:
            return

        self.sys.build_network()


if __name__ == "__main__":
    os.chdir("..")
    app = Vorpy()
    app.mainloop()

