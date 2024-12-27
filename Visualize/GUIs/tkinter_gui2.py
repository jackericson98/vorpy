import tkinter as tk
from tkinter import ttk, filedialog
from System.system import System
import os

class VorPyGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.sys = System(simple=True)

        self.title("VorPy")

        # Font classes
        self.fonts = {'class 1': ("Arial", 16), 'class 2': ("Arial", 10), 'class 3': ("Arial", 12, "bold"), 'class 4': ("Arial", 14)}

        # Title Section
        title_label = tk.Label(self, text="VorPy", font=self.fonts['class 1'], pady=10)
        title_label.pack()

        # Notebook for Main Sections
        self.main_notebook = ttk.Notebook(self)
        self.main_notebook.pack(expand=True, fill="both")

        # Information Section
        info_frame = ttk.Frame(self.main_notebook)
        self.main_notebook.add(info_frame, text="Information")

        self.create_information_section(info_frame)

        # Settings Section
        settings_frame = ttk.Frame(self.main_notebook)
        self.main_notebook.add(settings_frame, text="Settings")

        self.create_settings_section(settings_frame)

        # Run and Cancel Buttons
        button_frame = tk.Frame(self, pady=10)
        button_frame.pack()

        run_button = tk.Button(button_frame, text="Run", command=self.run_program, font=self.fonts['class 2'])
        run_button.pack(side="left", padx=5)

        cancel_button = tk.Button(button_frame, text="Cancel", command=self.quit, font=self.fonts['class 2'])
        cancel_button.pack(side="left", padx=5)

    def create_information_section(self, parent):
        # System Information
        sys_info_frame = ttk.LabelFrame(parent, text=" Information ")
        sys_info_frame.pack(fill="both", padx=10, pady=5)

        tk.Label(sys_info_frame, text="System Information", font=self.fonts['class 3'], anchor="w").grid(row=0, column=0, columnspan=3, sticky="w")

        # System Name
        tk.Label(sys_info_frame, text=self.sys.name, font=self.fonts['class 4']).grid(row=1, column=0, columnspan=2, sticky="w")

        # Ball File
        tk.Label(sys_info_frame, text="...", font=self.fonts['class 2']).grid(row=2, column=0, sticky="w")
        tk.Label(sys_info_frame, text="None", font=self.fonts['class 2']).grid(row=2, column=1, sticky="w")
        tk.Button(sys_info_frame, text="Ball File", command=self.choose_ball_file).grid(row=2, column=2, sticky="e")

        # Other Files
        tk.Label(sys_info_frame, text="Other Files:", font=self.fonts['class 2']).grid(row=3, column=0, sticky="w")
        tk.Label(sys_info_frame, text="None", font=self.fonts['class 2']).grid(row=3, column=1, sticky="w")

        # System Type
        tk.Label(sys_info_frame, text="Type:", font=self.fonts['class 2']).grid(row=4, column=0, sticky="w")
        tk.Label(sys_info_frame, text=self.sys.type.capitalize(), font=self.fonts['class 2']).grid(row=4, column=1, sticky="w")

        # Number of Balls
        tk.Label(sys_info_frame, text="# of Balls:", font=self.fonts['class 2']).grid(row=5, column=0, sticky="w")
        tk.Label(sys_info_frame, text=len(self.sys.balls), font=self.fonts['class 2']).grid(row=5, column=1, sticky="w")

        # System Data Placeholder
        tk.Label(sys_info_frame, text="System Data:", font=self.fonts['class 2']).grid(row=6, column=0, sticky="w")
        tk.Label(sys_info_frame, text="", font=self.fonts['class 2']).grid(row=6, column=1, sticky="w")

        # Output Directory
        tk.Label(sys_info_frame, text="...", font=self.fonts['class 2']).grid(row=7, column=0, sticky="w")
        tk.Label(sys_info_frame, text="None", font=self.fonts['class 2']).grid(row=7, column=1, sticky="w")
        tk.Button(sys_info_frame, text="Choose Output Directory", command=self.choose_output_directory).grid(row=7, column=2, sticky="e")

        # Group Information Placeholder
        group_info_frame = ttk.LabelFrame(parent, text=" Group Information ")
        group_info_frame.pack(fill="both", padx=10, pady=5)

        tk.Button(group_info_frame, text="Add Group", command=self.add_group).pack(anchor="ne")

    def create_settings_section(self, parent):
        # Build Settings
        build_settings_frame = ttk.LabelFrame(parent, text=" Build Settings ")
        build_settings_frame.pack(fill="both", padx=10, pady=5)

        # Add tabs for groups
        build_notebook = ttk.Notebook(build_settings_frame)
        build_notebook.pack(fill="both", expand=True)

        for group in self.sys.groups:
            group_frame = ttk.Frame(build_notebook)
            build_notebook.add(group_frame, text=group.name)

            # Network Type Dropdown
            tk.Label(group_frame, text="Network Type", font=self.fonts['class 2']).grid(row=0, column=0, sticky="w")
            network_type = ttk.Combobox(group_frame, values=["Additively Weighted", "Power", "Primitive"])
            network_type.set("Additively Weighted")
            network_type.grid(row=0, column=1, sticky="w")

            # Maximum Vertex
            tk.Label(group_frame, text="Maximum Vertex", font=self.fonts['class 2']).grid(row=1, column=0, sticky="w")
            max_vertex = tk.Entry(group_frame)
            max_vertex.insert(0, "15")
            max_vertex.grid(row=1, column=1, sticky="w")

            # Box Multiplier
            tk.Label(group_frame, text="Box Multiplier", font=self.fonts['class 2']).grid(row=2, column=0, sticky="w")
            box_multiplier = tk.Entry(group_frame)
            box_multiplier.insert(0, "1.25")
            box_multiplier.grid(row=2, column=1, sticky="w")

            # Atomic Radii/Masses Button
            tk.Label(group_frame, text="Atomic Radii/Masses", font=self.fonts['class 2']).grid(row=3, column=0, sticky="w")
            tk.Button(group_frame, text="Open", command=self.open_atomic_radii_gui).grid(row=3, column=1, sticky="w")

            # Surface Settings Button
            tk.Label(group_frame, text="Surface Settings", font=self.fonts['class 2']).grid(row=4, column=0, sticky="w")
            tk.Button(group_frame, text="Change", command=self.open_surface_settings_gui).grid(row=4, column=1, sticky="w")

        # Export Settings
        export_settings_frame = ttk.LabelFrame(parent, text=" Export Settings ")
        export_settings_frame.pack(fill="both", padx=10, pady=5)

        # Radio Buttons for Presets
        tk.Label(export_settings_frame, text="Presets", font=self.fonts['class 3']).pack(anchor="w")
        presets_frame = tk.Frame(export_settings_frame)
        presets_frame.pack(anchor="w")

        preset_var = tk.StringVar(value="Medium")

        def update_checkboxes():
            for checkbox, state in checkbox_states.items():
                checkbox.set(state[preset_var.get()])

        tk.Radiobutton(presets_frame, text="Small", variable=preset_var, value="Small", command=update_checkboxes).pack(side="left")
        tk.Radiobutton(presets_frame, text="Medium", variable=preset_var, value="Medium", command=update_checkboxes).pack(side="left")
        tk.Radiobutton(presets_frame, text="Large", variable=preset_var, value="Large", command=update_checkboxes).pack(side="left")
        tk.Radiobutton(presets_frame, text="All", variable=preset_var, value="All", command=update_checkboxes).pack(side="left")

        # Checkboxes for Output Options
        tk.Label(export_settings_frame, text="System Outputs", font=self.fonts['class 2']).pack(anchor="w")
        system_outputs_frame = tk.Frame(export_settings_frame)
        system_outputs_frame.pack(anchor="w")

        checkbox_states = {
            "PDB": {"Small": False, "Medium": True, "Large": True, "All": True},
            "Set Ball Radius": {"Small": False, "Medium": True, "Large": True, "All": True},
            "Info": {"Small": False, "Medium": True, "Large": True, "All": True},
            "Logs": {"Small": True, "Medium": True, "Large": True, "All": True},
            "Balls": {"Small": False, "Medium": True, "Large": True, "All": True},
            "Surrounding Balls": {"Small": False, "Medium": False, "Large": True, "All": True},
            "Shell: Surfaces": {"Small": False, "Medium": True, "Large": True, "All": True},
            "Shell: Edges": {"Small": False, "Medium": False, "Large": True, "All": True},
            "Shell: Vertices": {"Small": False, "Medium": False, "Large": True, "All": True},
            "Cells: Surfaces": {"Small": False, "Medium": False, "Large": True, "All": True},
            "Cells: Edges": {"Small": False, "Medium": False, "Large": True, "All": True},
            "Cells: Vertices": {"Small": False, "Medium": False, "Large": True, "All": True},
            "All: Surfaces": {"Small": False, "Medium": False, "Large": False, "All": True},
            "All: Edges": {"Small": False, "Medium": False, "Large": False, "All": True},
            "All: Vertices": {"Small": False, "Medium": False, "Large": False, "All": True}
        }

        checkboxes = {}
        for label_text in checkbox_states.keys():
            var = tk.BooleanVar()
            checkbox = tk.Checkbutton(system_outputs_frame, text=label_text, variable=var, anchor="w")
            checkbox.pack(anchor="w")
            checkboxes[label_text] = var

        update_checkboxes()

    def open_atomic_radii_gui(self):
        print("Opening Atomic Radii GUI...")

    def open_surface_settings_gui(self):
        print("Opening Surface Settings GUI...")

    def run_program(self):
        print("Running the program...")

    def choose_ball_file(self):
        self.sys.files['base_file'] = filedialog.askopenfilename(title='Choose Base File')
        print(f"Ball file selected: {self.sys.files['base_file']}")

    def choose_output_directory(self):
        self.sys.files['dir'] = filedialog.askdirectory(title='Choose Output Directory')
        print(f"Output directory selected: {self.sys.files['dir']}")

    def add_group(self):
        print("Adding a new group...")

if __name__ == "__main__":

    os.chdir('../..')
    app = VorPyGUI()
    app.mainloop()
