import tkinter as tk
from tkinter import ttk


def export_frame(self, parent):
    # Export Settings
    export_settings_frame = ttk.LabelFrame(parent, text=" Export Settings ")
    export_settings_frame.pack(fill="both", padx=10, pady=5)

    presets_frame = tk.Frame(export_settings_frame)
    presets_frame.pack(anchor="w")

    preset_var = tk.StringVar(value="Medium")

    checkbox_states = {
        "Small": {"Shell: Surfaces": False, "Shell: Edges": False, "Shell: Vertices": False,
                  "Cells: Surfaces": False, "Cells: Edges": False, "Cells: Vertices": False,
                  "All: Surfaces": False, "All: Edges": False, "All: Vertices": False},
        "Medium": {"Shell: Surfaces": True, "Shell: Edges": False, "Shell: Vertices": False,
                   "Cells: Surfaces": False, "Cells: Edges": False, "Cells: Vertices": False,
                   "All: Surfaces": True, "All: Edges": False, "All: Vertices": False},
        "Large": {"Shell: Surfaces": True, "Shell: Edges": True, "Shell: Vertices": True,
                  "Cells: Surfaces": True, "Cells: Edges": True, "Cells: Vertices": True,
                  "All: Surfaces": True, "All: Edges": True, "All: Vertices": True},
        "All": {"Shell: Surfaces": True, "Shell: Edges": True, "Shell: Vertices": True,
                "Cells: Surfaces": True, "Cells: Edges": True, "Cells: Vertices": True,
                "All: Surfaces": True, "All: Edges": True, "All: Vertices": True}
    }

    def update_checkboxes():
        for checkbox, state in checkbox_states[preset_var.get()].items():
            checkboxes[checkbox].set(state)

    tk.Radiobutton(presets_frame, text="Small", variable=preset_var, value="Small", command=update_checkboxes).pack(
        side="left")
    tk.Radiobutton(presets_frame, text="Medium", variable=preset_var, value="Medium", command=update_checkboxes).pack(
        side="left")
    tk.Radiobutton(presets_frame, text="Large", variable=preset_var, value="Large", command=update_checkboxes).pack(
        side="left")
    tk.Radiobutton(presets_frame, text="All", variable=preset_var, value="All", command=update_checkboxes).pack(
        side="left")

    # System Outputs Section
    tk.Label(export_settings_frame, text="System Outputs", font=self.fonts['class 2']).pack(anchor="w", pady=(10, 0))
    system_outputs_frame = tk.Frame(export_settings_frame)
    system_outputs_frame.pack(anchor="w", pady=(0, 10))

    # Grid for Shell, Cells, All options
    shell_cells_all_frame = ttk.LabelFrame(export_settings_frame, text=" Network Components ")
    shell_cells_all_frame.pack(fill="both", padx=10, pady=10)

    tk.Label(shell_cells_all_frame, text="", font=self.fonts['class 2'], width=5).grid(row=0, column=0)
    tk.Label(shell_cells_all_frame, text="Surfaces", font=self.fonts['class 2'], width=5).grid(row=0, column=1)
    tk.Label(shell_cells_all_frame, text="Edges", font=self.fonts['class 2'], width=5).grid(row=0, column=2)
    tk.Label(shell_cells_all_frame, text="Vertices", font=self.fonts['class 2'], width=5).grid(row=0, column=3)

    components = ["Shell", "Cells", "All"]
    checkboxes = {}

    for i, component in enumerate(components):
        tk.Label(shell_cells_all_frame, text=component, font=self.fonts['class 2'], width=15).grid(row=i + 1, column=0)
        for j, sub_component in enumerate(["Surfaces", "Edges", "Vertices"]):
            var = tk.BooleanVar()
            checkbox = tk.Checkbutton(shell_cells_all_frame, variable=var)
            checkbox.grid(row=i + 1, column=j + 1)
            checkboxes[f"{component}: {sub_component}"] = var

    update_checkboxes()
