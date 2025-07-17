import os
import sys
from pathlib import Path
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from PIL import Image, ImageTk

# Add the project root directory to the Python path
project_root = Path(__file__).resolve().parents[3]
sys.path.append(str(project_root))

from vorpy.src.system.system import System
from vorpy.src.GUI.system.system_frame import SystemFrame
from vorpy.src.GUI.group.groups_frame import GroupsFrame
from vorpy.src.GUI.help.help_window import HelpWindow
from vorpy.src.group import Group


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

        # Set the output directory
        self.output_dir = None

        # Set up the files dictionary
        self.files = {'sys_name': 'No File Loaded', 'base_file': '', 'other_files': [], 'dir': ''}
        self.exports = {'set_atoms': True, 'info': True, 'pdb': True, 'mol': False, 'cif': False, 'xyz': False,
                        'txt': False}
        self.radii_changes = []
        
        self.group_settings = {}

        # Title Section with Image Template
        title_frame = tk.Frame(self, pady=10)
        title_frame.pack(fill="x")

        # --- IMAGE TEMPLATE START ---
        # To use an image, place your image file (e.g., 'VorpyIcon.png') in the appropriate directory.
        # This template uses Pillow to allow resizing and background removal (transparency).
        # Make sure to install Pillow: pip install pillow
        try:
            # Load the image using a path relative to this file's location

            script_dir = os.path.dirname(os.path.abspath(__file__))
            img_path = os.path.join(script_dir, "Images", "VorpyIcon.png")
            img = Image.open(img_path)

            # Resize the image to be much smaller (e.g., 40x40 pixels)
            img = img.resize((50, 50), Image.LANCZOS)

            # If the image has a background and is PNG, try to remove it by converting white to transparent
            if img.mode != "RGBA":
                img = img.convert("RGBA")
            datas = img.getdata()
            newData = []
            for item in datas:
                # Detect white-ish pixels (tune threshold as needed)
                if item[0] > 240 and item[1] > 240 and item[2] > 240:
                    # Set alpha to 0 (transparent)
                    newData.append((255, 255, 255, 0))
                else:
                    newData.append(item)
            img.putdata(newData)

            self.logo_img = ImageTk.PhotoImage(img)

            # --- Set the window/taskbar icon ---
            # Try to set the icon for the window and taskbar
            # On Windows, .ico is best; on Linux, .png is usually fine
            # We'll try to use the same image for both, but .ico is preferred for best compatibility
            # If you have a .ico file, use it here
            script_dir = os.path.dirname(os.path.abspath(__file__))
            icon_path_ico = os.path.join(script_dir, "Images", "VorpyIcon.ico")
            
            try:
                if os.path.exists(icon_path_ico):
                    self.iconbitmap(icon_path_ico)
                else:
                    # For Linux/Mac, or if .ico not available, use .png
                    # Tkinter's iconphoto works cross-platform
                    self.iconphoto(True, self.logo_img)
                self.tray_icon = Image.open(icon_path_ico)
            except Exception as e:
                print("Could not set window/taskbar icon:", e)
        except Exception as e:
            # If Pillow is not available or image not found, just skip the image
            print("Logo image could not be loaded:", e)
        # --- IMAGE TEMPLATE END ---

        # Title and subtitle in a vertical frame, next to the image
        # Create a horizontal frame for logo and title
        logo_title_frame = tk.Frame(title_frame)
        logo_title_frame.pack(fill="x")
        logo_title_frame.grid_columnconfigure(0, weight=9)
        logo_title_frame.grid_columnconfigure(1, weight=10)

        # Place logo (if loaded) on the left, then title on the right, on the same row
        if hasattr(self, 'logo_img'):
            logo_label = tk.Label(logo_title_frame, image=self.logo_img, bg=title_frame.cget("bg"))
            logo_label.grid(row=0, column=0, padx=(0, 5), sticky="e")
        title_label = tk.Label(logo_title_frame, text="VorPy", font=self.fonts['title'])
        title_label.grid(row=0, column=1, sticky="w")
        
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
        self.system_frame = SystemFrame(self, frame)
        return self.system_frame

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
        self.output_dir = self.sys.files['dir']
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
            group.exports(info=True, shell_surfs=True, logs=True, concave_colors=build_settings['color_settings']['conc_col'])
        elif exports['size'] == 'Medium':
            group.exports(shell_surfs=True, surfs=True, shell_edges=True, edges=True, shell_verts=True, verts=True,
                            logs=True, atoms=True, surr_atoms=True, concave_colors=build_settings['color_settings']['conc_col'])
        elif exports['size'] == 'Large':
            # Export the group exports
            group.exports(shell_verts=True, shell_edges=True, shell_surfs=True, info=True, edges=True, verts=True,
                            atoms=True, surr_atoms=True, logs=True, atom_surfs=True, atom_edges=True, atom_verts=True,
                            concave_colors=build_settings['color_settings']['conc_col'])
        else:
            cust = exports['custom_settings']
            group.exports(info=cust['info'], logs=cust['logs'], atoms=cust['group_vars']['pdb'],
                          sep_surfs=cust['surfs_separate'], sep_edges=cust['edges_separate'],
                          sep_verts=cust['verts_separate'], atom_surfs=cust['surfs_cell'],
                          atom_edges=cust['edges_cell'], atom_verts=cust['verts_cell'], surfs=cust['surfs_all'],
                          edges=cust['edges_all'], verts=['verts_all'], shell_surfs=cust['surfs_shell'],
                          shell_edges=cust['edges_shell'], shell_verts=cust['verts_shell'],
                          surr_atoms=cust['surrounding_vars']['pdb'], concave_colors=build_settings['color_settings']['conc_col'])

    def run_program(self):
        """
        This sends a system to start running networks on all groups
        """

        # Create a group if None exists
        if len(self.group_settings) == 0:
            self.sys.create_group()

        # Set the output directory 
        self.sys.files['dir'] = self.output_dir

        # Update the radii changes in the system
        for change in self.radii_changes:
            self.sys.set_radii(change)

        # Create the groups with the correct settings
        for group_name in self.group_settings:
            self.run_group(group_name)

        # Export the system exports
        self.sys.exports(pdb=self.exports['pdb'], mol=self.exports['mol'], cif=self.exports['cif'],
                         xyz=self.exports['xyz'], txt=self.exports['txt'], info=self.exports['info'],
                         set_atoms=self.exports['set_atoms'])

        # Print where the files were exported to
        print(f"Files were exported to: {self.sys.files['dir']}")
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

