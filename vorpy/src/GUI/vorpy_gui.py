import os
import sys
import queue
import threading
import traceback
from itertools import combinations
from pathlib import Path

try:
    import tkinter as tk
    from tkinter import ttk
    from tkinter import filedialog, messagebox
    from PIL import Image, ImageTk

    HAS_TKINTER = True
except ImportError:
    HAS_TKINTER = False

# Add the project root directory to the Python path
project_root = Path(__file__).resolve().parents[3]
sys.path.append(str(project_root))

from vorpy.src.system.system import System
from vorpy.src.GUI.system.system_frame import SystemFrame
from vorpy.src.GUI.group.groups_frame import GroupsFrame
from vorpy.src.GUI.help.help_window import HelpWindow
from vorpy.src.GUI.support_functions import ToolTip
from vorpy.src.group import Group
from vorpy.src.inputs import read_verts


class VorPyGUI(tk.Tk):
    def __init__(self):
        # Initialize the parent class first
        super().__init__()

        # Create a default system
        self.sys = System(simple=True, name="No System Chosen")
        self.ball_file = None

        # Background solve thread
        self.solve_thread = None
        self.export_thread = None

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
        self.group_solves = {}
        self.interface_requests = []
        self.interface_mode = None

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
            datas = img.get_flattened_data()
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

        # Run status
        self.run_status = tk.StringVar(value="Ready")
        status_label = ttk.Label(self, textvariable=self.run_status)
        status_label.pack(pady=(5, 0))

        self.after(100, self.update_run_status)

        # Run progress
        progress_frame = ttk.Frame(self)
        progress_frame.pack(fill="x", padx=20, pady=(8, 2))

        self.run_status = tk.StringVar(value="Ready")
        self.run_percent = tk.StringVar(value="0.0 %")
        self.run_process = tk.StringVar(value="")

        status_frame = ttk.Frame(progress_frame)
        status_frame.pack(fill="x")

        ttk.Label(status_frame, textvariable=self.run_status, font=self.fonts['class 3']).pack(side="left")
        ttk.Label(status_frame, textvariable=self.run_percent, font=self.fonts['class 3']).pack(side="right")

        self.progress_bar = ttk.Progressbar(progress_frame, orient="horizontal", mode="determinate", maximum=100)
        self.progress_bar.pack(fill="x", pady=(4, 3))

        ttk.Label(progress_frame, textvariable=self.run_process, font=self.fonts['class 2']).pack(anchor="w")

        self.after(100, self.update_run_status)

        # Run and Cancel Buttons
        button_frame = tk.Frame(self, pady=10)
        button_frame.pack()

        help_button = ttk.Button(button_frame, text="Help", command=self.open_help)
        help_button.pack(side="left", padx=5)

        ToolTip(help_button, "Help Button\nLaunches a window to explain each function of the GUI")

        print_button = ttk.Button(button_frame, text="Print", command=self.print_system)
        print_button.pack(side="left", padx=5)

        ToolTip(print_button,
                "Print Button:\nPrints the system information and group\ninformation/settings in the command line.")

        interface_button = ttk.Button(button_frame, text="Interface(s)", command=self.interfaces)
        interface_button.pack(side="left", padx=5)

        ToolTip(
            interface_button,
            "Interface(s) Button\nSelect groups and create one or more interface requests."
        )

        export_button = ttk.Button(button_frame, text="Export", command=self.start_export)
        export_button.pack(side="right", padx=5)

        ToolTip(export_button,
                "Export Button\nExports all solved groups to the system directory.\nGroups will be given individual sub directories\nwithin the directory")

        run_button = ttk.Button(button_frame, text="Solve", command=self.run_program)
        run_button.pack(side="right", padx=5)

        ToolTip(run_button, "Solve Button\nSolves all groups in series with a \nprogress bar in the command line")

    def start_export(self):
        if self.export_thread is not None and self.export_thread.is_alive():
            print("VorPy is already exporting files.")
            return

        self.sys.start_run()

        self.export_thread = threading.Thread(
            target=self._export_worker,
            daemon=True,
        )
        self.export_thread.start()

    def _export_worker(self):
        try:
            self.export_files()
            self.sys.finish_run()
        except Exception as e:
            self.sys.run_active = False
            self.sys.run_process = "Error"
            print(f"\nVorPy export failed: {e}")
            raise

    def queue_progress_update(self, state):
        self.progress_queue.put(("progress", state))

    def process_progress_queue(self):
        try:
            while True:
                event, data = self.progress_queue.get_nowait()

                if event == "progress":
                    self.handle_progress_update(data)

                elif event == "complete":
                    self.handle_solve_complete()

                elif event == "error":
                    self.handle_solve_error(data)

        except queue.Empty:
            pass

        self.after(50, self.process_progress_queue)

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

    def get_group_by_name(self, group_name):
        """Return the Group object stored in self.sys.groups by name."""
        if self.sys.groups is None:
            return None

        for group in self.sys.groups:
            if group.name == group_name:
                return group

        return None

    def choose_output_directory(self):
        self.sys.files['dir'] = filedialog.askdirectory(title='Choose Output Directory')
        self.output_dir = self.sys.files['dir']
        print(f"Output directory selected: {self.sys.files['dir']}")

    def get_interface_group_names(self):
        """Return unique names from the currently available system groups."""
        if not self.sys.groups:
            return []

        names = []
        for group in self.sys.groups:
            if isinstance(group, str):
                name = group
            else:
                name = getattr(group, "name", str(group))

            if name and name not in names:
                names.append(name)

        return names

    def build_all_interface_pairs(self, group_names):
        """Build interface pairs according to the GUI's All-selection rules."""
        if len(group_names) == 0:
            return [("(default)", "surrounding")]

        if len(group_names) == 1:
            return [(group_names[0], "surrounding")]

        if len(group_names) == 2:
            return [(group_names[0], group_names[1])]

        return list(combinations(group_names, 2))

    def interfaces(self):

        dialog = tk.Toplevel(self)
        dialog.title("Interface")
        dialog.resizable(False, False)
        dialog.transient(self)
        dialog.grab_set()

        frame = ttk.Frame(dialog, padding=18)
        frame.grid(row=0, column=0, sticky="nsew")

        ttk.Label(
            frame,
            text="Interface",
            font=self.fonts['class 1']
        ).grid(row=0, column=0, columnspan=3, pady=(0, 3))

        ttk.Label(
            frame,
            text="Select groups for interface",
            font=self.fonts['class 2']
        ).grid(row=1, column=0, columnspan=3, pady=(0, 14))

        system_group_names = self.get_interface_group_names()
        pending_group_names = list(self.group_settings.keys())

        group_names = list(
            dict.fromkeys([
                *system_group_names,
                *pending_group_names,
            ])
        )

        first_options = ["(default)", *group_names]
        second_options = ["surrounding", "solvent", *group_names]

        # Preserve order while avoiding repeated entries.
        first_options = list(dict.fromkeys(first_options))
        second_options = list(dict.fromkeys(second_options))

        first_var = tk.StringVar(value=first_options[0])
        second_var = tk.StringVar(value=second_options[0])

        first_dropdown = ttk.Combobox(
            frame,
            textvariable=first_var,
            values=first_options,
            state="readonly",
            width=24
        )
        first_dropdown.grid(row=2, column=0, padx=(0, 6), sticky="ew")

        second_dropdown = ttk.Combobox(
            frame,
            textvariable=second_var,
            values=second_options,
            state="readonly",
            width=24
        )
        second_dropdown.grid(row=2, column=1, columnspan=2, padx=(6, 0), sticky="ew")

        ttk.Label(
            frame,
            text="Create Interface(s)",
            font=self.fonts['class 3']
        ).grid(row=3, column=0, columnspan=3, pady=(18, 8))

        def finish(mode, interfaces):
            self.interface_mode = mode
            self.interface_requests = interfaces

            dialog.destroy()

        def create_selected():
            group_1 = first_var.get()
            group_2 = second_var.get()

            if group_1 == group_2:
                messagebox.showwarning(
                    "Invalid Interface",
                    "The two interface selections must be different.",
                    parent=dialog
                )
                return

            finish("selected", [(group_1, group_2)])

        def create_all():
            finish("all", [])

        ttk.Button(
            frame,
            text="Cancel",
            command=dialog.destroy
        ).grid(row=4, column=0, sticky="w")

        ttk.Button(
            frame,
            text="All",
            command=create_all
        ).grid(row=4, column=1, padx=6)

        ttk.Button(
            frame,
            text="Selected",
            command=create_selected
        ).grid(row=4, column=2, sticky="e")

        for column in range(3):
            frame.grid_columnconfigure(column, weight=1)

        dialog.protocol("WM_DELETE_WINDOW", dialog.destroy)
        dialog.update_idletasks()

        x = self.winfo_rootx() + (self.winfo_width() - dialog.winfo_width()) // 2
        y = self.winfo_rooty() + (self.winfo_height() - dialog.winfo_height()) // 2
        dialog.geometry(f"+{x}+{y}")

    def run_group(self, group_name, verts=None, make_net=True, build_net=True):

        settings = self.group_settings[group_name]
        build_settings = settings['build_settings'].get_settings()

        net_type_dict = {'Additively Weighted': 'aw', 'Power': 'pow', 'Primitive': 'prm'}

        group = Group(
            self.sys,
            name=group_name,
            atoms=settings['selections'].selections['balls'],
            residues=settings['selections'].selections['residues'],
            chains=settings['selections'].selections['chains'],
            molecules=settings['selections'].selections['molecules'],
            make_net=make_net,
            build_net=build_net,
            surf_res=float(build_settings['color_settings']['surf_res']),
            box_size=float(build_settings['box_size']),
            max_vert=float(build_settings['max_vert']),
            net_type=net_type_dict[build_settings['net_type']],
            surf_col=build_settings['color_settings']['surf_col'],
            surf_scheme=build_settings['color_settings']['surf_scheme'],
            scheme_factor=build_settings['color_settings']['surf_fact'],
            vert_col=build_settings['color_settings']['vert_col'],
            edge_col=build_settings['color_settings']['edge_col'],
            verts=verts
        )

        if self.sys.groups is None:
            self.sys.groups = []

        if group not in self.sys.groups:
            self.sys.groups.append(group)

        self.group_solves[group_name] = bool(build_net)

        return group

    def create_interface_group_definitions(self):
        """
        Create the GUI-defined Group objects without constructing or building
        their complete networks.
        """
        # print("\n" + "=" * 70)
        # print("CREATING INTERFACE GROUP DEFINITIONS")
        # print("=" * 70)

        if self.sys.groups is None:
            self.sys.groups = []

        groups_by_name = {
            group.name: group
            for group in self.sys.groups
            if not isinstance(group, str)
        }

        for group_name in self.group_settings:

            verts = None

            if self.files["other_files"]:
                if "verts" in self.files["other_files"][0]:
                    verts = self.files["other_files"][0]

            group = self.run_group(
                group_name=group_name,
                verts=verts,
                make_net=False,
                build_net=False,
            )

            groups_by_name[group_name] = group

        return groups_by_name

    def resolve_gui_interface_pairs(self, groups_by_name):
        """
        Convert the GUI interface selection into Group object pairs.

        All mode:
            Build every unique pairwise combination of available GUI groups.

        Selected mode:
            Build only the pair selected in the interface dialog.

        System.make_interfaces() expects:
            [(group1_object, group2_object_or_None), ...]
        """
        interface_pairs = []

        # --------------------------------------------------------------
        # All groups
        # --------------------------------------------------------------
        if self.interface_mode == "all":
            groups = list(groups_by_name.values())

            if len(groups) == 0:
                raise ValueError(
                    "Cannot create interfaces because no GUI groups were defined."
                )

            if len(groups) == 1:
                interface_pairs = [
                    (groups[0], None)
                ]
            else:
                interface_pairs = list(
                    combinations(groups, 2)
                )

        # --------------------------------------------------------------
        # Selected groups
        # --------------------------------------------------------------
        elif self.interface_mode == "selected":
            if len(self.interface_requests) != 1:
                raise ValueError(
                    "Selected interface mode requires exactly one interface request."
                )

            first_name, second_name = self.interface_requests[0]

            # Resolve the first side.
            if first_name == "(default)":
                if not groups_by_name:
                    raise ValueError(
                        "The default interface requires at least one group definition."
                    )

                group1 = next(iter(groups_by_name.values()))
            else:
                group1 = groups_by_name.get(first_name)

            if group1 is None:
                raise ValueError(
                    f"Could not resolve interface group '{first_name}'. "
                    f"Available groups: {list(groups_by_name.keys())}"
                )

            # Resolve the second side.
            if second_name in ("surrounding", "solvent"):
                group2 = None
            else:
                group2 = groups_by_name.get(second_name)

            if second_name not in ("surrounding", "solvent") and group2 is None:
                raise ValueError(
                    f"Could not resolve interface group '{second_name}'. "
                    f"Available groups: {list(groups_by_name.keys())}"
                )

            if group1 is group2:
                raise ValueError(
                    f"Cannot create an interface between group "
                    f"'{group1.name}' and itself."
                )

            interface_pairs = [
                (group1, group2)
            ]

        elif self.interface_mode is None:
            pass
        else:
            raise ValueError(
                f"Unknown GUI interface mode: {self.interface_mode}"
            )

        return interface_pairs

    def update_run_status(self):
        elapsed = self.sys.get_run_time()
        h = int(elapsed // 3600)
        m = int((elapsed % 3600) // 60)
        s = elapsed % 60

        if self.sys.run_start is None:
            self.run_status.set("Ready")
            self.run_percent.set("0.0 %")
            self.run_process.set("")
            self.progress_bar["value"] = 0
        else:
            self.run_status.set(f"Run Time: {h}:{m:02d}:{s:05.2f}")
            self.run_percent.set(f"{self.sys.run_progress:.1f} %")
            self.progress_bar["value"] = self.sys.run_progress

            if self.sys.run_network:
                self.run_process.set(f"Network: {self.sys.run_network}  |  {self.sys.run_process}")
            else:
                self.run_process.set(self.sys.run_process)

        self.after(100, self.update_run_status)

    def run_program(self):
        """Start the solve in a background thread."""

        if self.solve_thread is not None and self.solve_thread.is_alive():
            print("VorPy is already solving a system.")
            return

        self.sys.start_run()

        self.solve_thread = threading.Thread(target=self._run_program_worker, daemon=True)
        self.solve_thread.start()

    def _run_program_worker(self):
        try:
            self._run_program()
            self.sys.finish_run()
        except Exception as e:
            self.sys.run_active = False
            self.sys.run_process = "Error"
            print(f"\nVorPy solve failed: {e}")
            raise

    def _run_program(self):
        """
        Solve all groups without exporting.
        """
        self.sys.files['dir'] = self.output_dir

        for change in self.radii_changes:
            self.sys.set_radii(change)

        # No GUI groups specified -> solve the entire system
        if len(self.group_settings) == 0:
            self.sys.create_group()
            group = self.sys.groups[-1]
            group.build()
            return self.sys

        # --------------------------------------------------------------
        # Interface-only build
        # --------------------------------------------------------------
        if self.interface_mode is not None:
            groups_by_name = self.create_interface_group_definitions()
            interface_pairs = self.resolve_gui_interface_pairs(groups_by_name=groups_by_name)
            self.sys.make_interfaces(interface_pairs)
            return self.sys

        # --------------------------------------------------------------
        # Normal GUI group builds
        # --------------------------------------------------------------
        for group_name in self.group_settings:
            verts = None

            if self.files["other_files"]:
                if "verts" in self.files["other_files"][0]:
                    verts = self.files["other_files"][0]

            self.run_group(group_name, verts)

        return self.sys

    def export_single_group(self, group, settings):
        """
        Export an already solved or loaded group using that group's export settings.
        """
        build_settings = settings['build_settings'].get_settings()
        exports = settings['export_settings'].get_settings()

        # Resolve parent export directory
        export_parent = exports.get("directory")

        if export_parent in {None, "", "Default Output Directory"}:
            export_parent = self.sys.files.get("dir")

        if export_parent is None:
            raise ValueError("No valid export directory is set.")

        os.makedirs(export_parent, exist_ok=True)

        # System files go to the parent export directory
        self.sys.files["dir"] = export_parent

        # Group files go to parent/group_name
        group_dir = os.path.join(export_parent, group.name)
        os.makedirs(group_dir, exist_ok=True)

        group.dir = group_dir

        concave_colors = build_settings['color_settings']['conc_col']
        round_to = exports.get("round_to", 3)
        file_type = exports.get("file_type", "off")

        if exports['size'] == 'Small':
            group.exports(
                info=True,
                shell_surfs=True,
                logs=True,
                concave_colors=concave_colors,
                round_to=round_to,
                file_type=file_type,
            )

        elif exports['size'] == 'Medium':
            group.exports(
                shell_surfs=True,
                surfs=True,
                shell_edges=True,
                edges=True,
                shell_verts=True,
                verts=True,
                logs=True,
                atoms=True,
                surr_atoms=True,
                concave_colors=concave_colors,
                round_to=round_to,
                file_type=file_type,
            )

        elif exports['size'] == 'Large':
            group.exports(
                shell_verts=True,
                shell_edges=True,
                shell_surfs=True,
                info=True,
                edges=True,
                verts=True,
                atoms=True,
                surr_atoms=True,
                logs=True,
                atom_surfs=True,
                atom_edges=True,
                atom_verts=True,
                concave_colors=concave_colors,
                round_to=round_to,
                file_type=file_type,
            )

        else:
            cust = exports['custom_settings']

            group.exports(
                info=cust['info'],
                logs=cust['logs'],
                atoms=cust['group_vars']['pdb'],
                sep_surfs=cust['surfs_separate'],
                sep_edges=cust['edges_separate'],
                sep_verts=cust['verts_separate'],
                atom_surfs=cust['surfs_cell'],
                atom_edges=cust['edges_cell'],
                atom_verts=cust['verts_cell'],
                surfs=cust['surfs_all'],
                edges=cust['edges_all'],
                verts=cust['verts_all'],
                shell_surfs=cust['surfs_shell'],
                shell_edges=cust['edges_shell'],
                shell_verts=cust['verts_shell'],
                surr_atoms=cust['surrounding_vars']['pdb'],
                concave_colors=concave_colors,
                round_to=round_to,
                file_type=file_type,
            )

        print(f"Exported group: {group.name}")

    def export_interface_files(self):
        """
        Export all completed interface networks and definition-level information
        for the groups participating in those interfaces.
        """
        interfaces = getattr(self.sys, "ifaces", None) or []

        built_interfaces = [
            iface
            for iface in interfaces
            if iface is not None and iface.net is not None
        ]

        if not built_interfaces:
            messagebox.showwarning(
                "No Interfaces Solved",
                "No solved interface networks are available to export.\n\n"
                "Please solve the requested interfaces before exporting.",
            )
            return False

        export_parent = self.sys.files.get("dir")

        # Use the first interface's Group 1 GUI export directory when available.
        first_group_name = built_interfaces[0].group1.name
        first_group_settings = self.group_settings.get(first_group_name)

        if first_group_settings is not None:
            export_settings = (
                first_group_settings["export_settings"].get_settings()
            )

            requested_directory = export_settings.get("directory")

            if requested_directory not in {
                None,
                "",
                "Default Output Directory",
            }:
                export_parent = requested_directory

        if export_parent is None:
            raise ValueError("No valid interface export directory is set.")

        os.makedirs(export_parent, exist_ok=True)
        self.sys.files["dir"] = export_parent

        for interface_index, iface in enumerate(built_interfaces):
            group_settings = self.group_settings.get(iface.group1.name)

            round_to = 3

            if group_settings is not None:
                export_settings = (
                    group_settings["export_settings"].get_settings()
                )
                round_to = export_settings.get("round_to", 3)

            # Set a deterministic interface directory rather than allowing
            # set_dir() to create a new numbered directory on every export.
            iface.dir = os.path.join(
                export_parent,
                iface.name,
            )

            os.makedirs(iface.dir, exist_ok=True)

            iface.export(
                all_=True,
                group_info=True,
                round_to=round_to,
            )

        return True

    def export_group_files(self, group_name=None, export_system=False):
        """
        Export one group if group_name is provided, otherwise export all solved/loaded groups.
        """
        if group_name is not None:
            group_names = [group_name]
        else:
            group_names = list(self.group_settings.keys())

        for name in group_names:
            group = self.get_group_by_name(name)

            solved = self.group_solves.get(name, False)
            loaded_from_logs = group is not None and group.net is not None and getattr(group.net, "loaded_from_logs",
                                                                                       False)

            if group is None or group.net is None or (not solved and not loaded_from_logs):
                messagebox.showwarning(
                    "Group Not Solved",
                    f"Group '{name}' has not been solved or loaded from logs yet.\n\n"
                    f"Please solve the group before exporting."
                )
                continue

            settings = self.group_settings[name]
            self.export_single_group(group, settings)

            group_dir = group.dir if group.dir is not None else self.sys.files['dir']
            print(f"Group '{group.name}' files were exported to: {group_dir}")

        if export_system:
            self.sys.exports(
                pdb=self.exports['pdb'],
                mol=self.exports['mol'],
                cif=self.exports['cif'],
                xyz=self.exports['xyz'],
                txt=self.exports['txt'],
                info=self.exports['info'],
                set_atoms=self.exports['set_atoms'],
            )

    def export_files(self):
        """
        Export either completed interfaces or completed full groups, followed by
        system-level files.
        """
        interfaces = getattr(self.sys, "ifaces", None) or []

        built_interfaces = [
            iface
            for iface in interfaces
            if iface is not None and iface.net is not None
        ]

        if built_interfaces:
            interface_exported = self.export_interface_files()

            if not interface_exported:
                return

        else:
            self.export_group_files(
                group_name=None,
                export_system=False,
            )

        self.sys.exports(
            pdb=self.exports['pdb'],
            mol=self.exports['mol'],
            cif=self.exports['cif'],
            xyz=self.exports['xyz'],
            txt=self.exports['txt'],
            info=self.exports['info'],
            set_atoms=self.exports['set_atoms'],
        )

        print(
            f"System files were exported to: "
            f"{self.sys.files['dir']}"
        )

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
