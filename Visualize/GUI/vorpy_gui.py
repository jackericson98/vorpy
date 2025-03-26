import os
import sys
from pathlib import Path

# Add the project root directory to the Python path
project_root = Path(__file__).resolve().parents[2]
sys.path.append(str(project_root))

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from System.system import System
from System.Group.group import Group
from Visualize.GUI.info.info_frame import SystemFrame
from Visualize.GUI.settings.group_settings_frame import GroupSettingsFrame
from Visualize.GUI.settings.surface.color_settings_window import ColorSettingsWindow
from Visualize.GUI.help.help_window import HelpWindow

"""
This GUI operates the whole VorPy interface. Once running, the command line that it was run out of will inform you on 
your progress. The GUI is only for launching the program.

GUI Options:

GUI returns
"""


class VorPyGUI(tk.Tk):
    def __init__(self, system=None):
        # Initialize the parent class first
        super().__init__()
        
        # Create a default system
        self.sys = System(simple=True, name="No System Chosen")
        
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

        # Set up the files dictionary
        self.files = {'sys_name': 'No File Loaded', 'base_file': '', 'other_files': [], 'dir': '' }
        
        # Set up the dictionaries for build and export settings
        self.build_settings = {
            'max_vert': 40,
            'box_size': 1.25,
            'net_type': 'aw',
            'color_settings': {
                'surf_col': 'plasma',
                'surf_scheme': 'mean_curv',
                'surf_fact': 'log',
                'vert_col': 'red',
                'edge_col': 'grey'
            }
        }
        self.export_settings = {'dir': None, 'type': 'med', 'other':{}}

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
        
        # Settings Frame (Full Width)
        settings_frame = tk.Frame(self)
        settings_frame.pack(expand=True, fill="both", padx=10, pady=(0, 10))
        
        # Create group settings section
        self.group_settings_frame = GroupSettingsFrame(settings_frame, self)
        self.group_settings_frame.pack(fill="both", expand=True)
        
        # Run and Cancel Buttons
        button_frame = tk.Frame(self, pady=10)
        button_frame.pack()
        
        run_button = ttk.Button(button_frame, text="Run", command=self.run_program)
        run_button.pack(side="left", padx=5)
        
        cancel_button = ttk.Button(button_frame, text="Cancel", command=self.quit)
        cancel_button.pack(side="left", padx=5)
        
        help_button = ttk.Button(button_frame, text="Help", command=self.open_help)
        help_button.pack(side="left", padx=5)

        print_button = ttk.Button(button_frame, text="Print", command=self.print_system)
        print_button.pack(side="left", padx=5)

    def create_information_section(self, frame):
        SystemFrame(self, frame)

    def open_surface_settings_gui(self):
        """Open the surface settings window."""
        ColorSettingsWindow(self)

    def choose_ball_file(self):
        """Open file dialog to select a ball file."""
        filename = filedialog.askopenfilename(
            title="Select Ball File",
            filetypes=[("Ball files", "*.pdb"), ("All files", "*.*")]
        )
        if filename:
            self.ball_file = filename
            self.sys.ball_file = filename
            self.sys.name = os.path.basename(filename)  # Update system name to filename
            self.files['sys_name'].set(self.sys.name)  # Update the display

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

    def print_system(self):
        """Print the system."""
        print(self.files)
        print(self.build_settings)
        print(self.export_settings)
        print(self.group_settings_frame.get_current_group_settings())

    def update_surface_settings_display(self):
        """Update the display of surface settings in the main GUI."""
        # Update the surface settings display in the build frame
        if hasattr(self, 'build_frame'):
            self.build_frame.update_surface_settings_display()


if __name__ == "__main__":
    os.chdir('../..')
    # create the system
    sys = System(file="./Data/test_data/cambrin.pdb", name="Test System")
    app = VorPyGUI(sys)
    app.mainloop()

