import os
import sys
from pathlib import Path

# Add the project root directory to the Python path
project_root = Path(__file__).resolve().parents[2]
sys.path.append(str(project_root))

import tkinter as tk
from tkinter import filedialog
from System.system import System
from System.Group.group import Group
from Visualize.GUI.info.info_frame import SystemFrame
from Visualize.GUI.groups.groups_frame import GroupsFrame
from Visualize.GUI.settings.settings_frame import create_settings_section
from Visualize.GUI.settings.surface.surface_settings_window import SurfaceSettingsWindow
from Visualize.GUI.help.help_window import HelpWindow

"""
This GUI operates the whole VorPy interface. Once running, the command line that it was run out of will inform you on 
your progress. The GUI is only for launching the program.

GUI Options:

GUI returns
"""


class VorPyGUI(tk.Tk):
    def __init__(self):
        # Initialize the parent class first
        super().__init__()
        
        # Create a default system
        self.sys = System(simple=True)
        
        # Set window title
        self.title("VorPy")
        
        # Font classes
        self.fonts = {
            'title': ("Arial", 24, "bold"),
            'subtitle': ("Arial", 12),
            'class 1': ("Arial", 16),
            'class 2': ("Arial", 10),
            'class 3': ("Arial", 12, "bold"),
            'class 4': ("Arial", 14)
        }
        
        # Title Section
        title_frame = tk.Frame(self, pady=10)
        title_frame.pack(fill="x")
        
        title_label = tk.Label(title_frame, text="VorPy", font=self.fonts['title'])
        title_label.pack()
        
        subtitle_label = tk.Label(title_frame, text="Comprehensive Voronoi Diagram Calculation Tool", 
                                font=self.fonts['subtitle'])
        subtitle_label.pack(pady=(0, 10))
        
        # System Information Section (Full Width)
        self.info_frame = tk.Frame(self, height=200)
        self.info_frame.pack(fill="x", padx=10, pady=(0, 10))
        self.create_information_section(self.info_frame)
        
        # Selection Frame (Holds Groups and Settings)
        selection_frame = tk.Frame(self)
        selection_frame.pack(expand=True, fill="both", padx=10, pady=(0, 10))
        
        # Configure grid weights for selection frame
        selection_frame.grid_columnconfigure(0, weight=1)  # Groups section
        selection_frame.grid_columnconfigure(1, weight=1)  # Settings section
        selection_frame.grid_rowconfigure(0, weight=1)
        
        # Groups Section (Left Column)
        groups_frame = tk.Frame(selection_frame)
        groups_frame.grid(row=0, column=0, padx=(0, 10), sticky="nsew")
        self.groups_frame = GroupsFrame(self, groups_frame)
        
        # Settings Section (Right Column)
        settings_frame = tk.Frame(selection_frame)
        settings_frame.grid(row=0, column=1, sticky="nsew")
        
        # Create settings section
        create_settings_section(self, settings_frame)
        
        # Run and Cancel Buttons
        button_frame = tk.Frame(self, pady=10)
        button_frame.pack()
        
        run_button = tk.Button(button_frame, text="Run", command=self.run_program, font=self.fonts['class 2'])
        run_button.pack(side="left", padx=5)
        
        cancel_button = tk.Button(button_frame, text="Cancel", command=self.quit, font=self.fonts['class 2'])
        cancel_button.pack(side="left", padx=5)
        
        help_button = tk.Button(button_frame, text="Help", command=self.open_help, font=self.fonts['class 2'])
        help_button.pack(side="left", padx=5)

    def create_information_section(self, frame):
        SystemFrame(self, frame)

    def create_settings_section(self, frame):
        create_settings_section(self, frame)

    def open_surface_settings_gui(self):
        """Open the surface settings window."""
        SurfaceSettingsWindow(self)

    def choose_ball_file(self):
        self.sys.files['base_file'] = filedialog.askopenfilename(title='Choose Base File')
        self.info_frame.update()
        print(f"Ball file selected: {self.sys.files['base_file']}")

    def choose_output_directory(self):
        self.sys.files['dir'] = filedialog.askdirectory(title='Choose Output Directory')
        print(f"Output directory selected: {self.sys.files['dir']}")

    def add_group(self):
        print("Adding a new group...")

    def run_program(self):
        """
        This sends a system to start running networks on all groups
        """
        return self.sys

    def open_help(self):
        """Open the help window."""
        HelpWindow(self)


if __name__ == "__main__":
    os.chdir('../..')
    app = VorPyGUI()
    app.mainloop()

