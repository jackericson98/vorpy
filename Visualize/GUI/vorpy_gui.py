import os
import sys
from pathlib import Path
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog

# Add the project root directory to the Python path
project_root = Path(__file__).resolve().parents[2]
sys.path.append(str(project_root))

from System.system import System
from Visualize.GUI.system.system_frame import SystemFrame
from Visualize.GUI.group.groups_frame import GroupsFrame
from Visualize.GUI.help.help_window import HelpWindow
from System.Group.group import Group


class VorPyGUI(tk.Tk):
    def __init__(self):
        # Initialize the parent class first
        super().__init__()
        
        # Create a default system
        self.sys = System(simple=True, name="No System Chosen")
        self.ball_file = None
        
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
        self.exports = {'set_atoms': True, 'info': True, 'pdb': True, 'mol': False, 'cif': False, 'xyz': False,
                        'txt': False}
        self.radii_changes = []
        
        self.group_settings = {}

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
        self.info_frame.pack(fill="x", padx=5)
        self.create_information_section(self.info_frame)
        
        # Settings Frame (Full Width)
        settings_frame = tk.Frame(self)
        settings_frame.pack(expand=True, fill="both", padx=10)
        
        # Create group settings section
        self.group_settings_frame = GroupsFrame(settings_frame, self, self.group_settings)
        self.group_settings_frame.pack(fill="both", expand=True)
        
        # Run and Cancel Buttons
        button_frame = tk.Frame(self, pady=10)
        button_frame.pack()

        help_button = ttk.Button(button_frame, text="Help", command=self.open_help)
        help_button.pack(side="left", padx=5)
        
        print_button = ttk.Button(button_frame, text="Print", command=self.print_system)
        print_button.pack(side="left", padx=5)

        run_button = ttk.Button(button_frame, text="Run All", command=self.run_program)
        run_button.pack(side="right", padx=5)

        cancel_button = ttk.Button(button_frame, text="Cancel", command=self.quit)
        cancel_button.pack(side="right", padx=5)

    def create_information_section(self, frame):
        SystemFrame(self, frame)

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
            self.files['sys_name'].set(self.sys.name.upper())  # Update the display

    def choose_output_directory(self):
        self.sys.files['dir'] = filedialog.askdirectory(title='Choose Output Directory')
        print(f"Output directory selected: {self.sys.files['dir']}")

    def run_group(self, group_name):
        """
        This runs a group from the group settings dictionary
        """
        settings = self.group_settings[group_name]
        build_settings = settings['build_settings'].get_settings()

        # Create a dictionary to convert the net type to something that can be interpreted
        net_type_dict = {'Additively Weighted': 'aw', 'Power': 'pow', 'Primitive': 'prm'}

        # Create the group
        group = Group(
            self.sys,
            name=group_name,
            atoms=settings['selections'].selections['balls'],
            residues=settings['selections'].selections['residues'],
            chains=settings['selections'].selections['chains'],
            molecules=settings['selections'].selections['molecules'],
            build_net=True,
            surf_res=float(build_settings['color_settings']['surf_res']),
            box_size=float(build_settings['box_size']),
            max_vert=float(build_settings['max_vert']),
            net_type=net_type_dict[build_settings['net_type']],
            surf_col=build_settings['color_settings']['surf_col'],
            surf_scheme=build_settings['color_settings']['surf_scheme'],
            scheme_factor=build_settings['color_settings']['surf_fact'],
            vert_col=build_settings['color_settings']['vert_col'],
            edge_col=build_settings['color_settings']['edge_col']
        )

        # Export the group
        exports = settings['export_settings'].get_settings()
        # Check if the export directory is chosen
        if exports['directory'] == 'Default Output Directory' or not os.path.exists(exports['directory']):
            exports['directory'] = None
        # Set the group's directory
        group.dir = exports['directory']

        # If the size is not custom export the given size information
        if exports['size'] == 'Small':
            group.exports(info=True, shell_surfs=True, logs=True)
        elif exports['size'] == 'Medium':
            group.exports(shell_surfs=True, surfs=True, shell_edges=True, edges=True, shell_verts=True, verts=True,
                            logs=True, atoms=True, surr_atoms=True)
        elif exports['size'] == 'Large':
            # Export the group exports
            group.exports(shell_verts=True, shell_edges=True, shell_surfs=True, info=True, edges=True, verts=True,
                            atoms=True, surr_atoms=True, logs=True, atom_surfs=True, atom_edges=True, atom_verts=True)
        else:
            cust = exports['custom_settings']
            group.exports(info=cust['info'], logs=cust['logs'], atoms=cust['group_vars']['pdb'],
                          sep_surfs=cust['surfs_separate'], sep_edges=cust['edges_separate'],
                          sep_verts=cust['verts_separate'], atom_surfs=cust['surfs_cell'],
                          atom_edges=cust['edges_cell'], atom_verts=cust['verts_cell'], surfs=cust['surfs_all'],
                          edges=cust['edges_all'], verts=['verts_all'], shell_surfs=cust['surfs_shell'],
                          shell_edges=cust['edges_shell'], shell_verts=cust['verts_shell'],
                          surr_atoms=cust['surrounding_vars']['pdb'])

    def run_program(self):
        """
        This sends a system to start running networks on all groups
        """

        # Create a group if None exists
        if len(self.group_settings) == 0:
            self.sys.create_group()

        # Create the groups with the correct settings
        for group_name in self.group_settings:
            self.run_group(group_name)
        # Export the system exports
        self.sys.exports(pdb=self.exports['pdb'], mol=self.exports['mol'], cif=self.exports['cif'],
                         xyz=self.exports['xyz'], txt=self.exports['txt'], info=self.exports['info'],
                         set_atoms=self.exports['set_atoms'])

        return self.sys

    def open_help(self):
        """Open the help window."""
        HelpWindow(self)

    def print_system(self):
        """Print the system."""
        print(self.files)
        print(self.exports)
        for group in self.group_settings:
            print(group)
            print(self.group_settings[group]['build_settings'].get_settings())
            print(self.group_settings[group]['export_settings'].get_settings())
            print(self.group_settings[group]['selections'].selections)
        print(self.radii_changes)

    def update_surface_settings_display(self):
        """Update the display of surface settings in the main GUI."""
        # Update the surface settings display in the build frame
        if hasattr(self, 'build_frame'):
            self.build_frame.update_surface_settings_display()


if __name__ == "__main__":
    os.chdir('../..')
    # create the system
    app = VorPyGUI()
    app.mainloop()

