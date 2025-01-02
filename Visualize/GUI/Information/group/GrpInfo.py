import tkinter as tk
from tkinter import ttk


"""
This file updates the group info. If you add more than one group, the gui will update it.

Options include:
    Defaults:
        Centered Text saying no molecule is loaded, but spaced out for options
    Defaults Once a file is loaded:
        1. Molecule: No Sol - this means as soon as a molecule file is loaded anything marked SOL,
           HOH, WAT, etc... is excluded
        2. Foam:


"""


def group_frame(root, parent):
    """
    Builds the frame for group configuration within the parent frame.

    Args:
        root: The main GUI application object, containing system and font settings.
        parent: The parent frame where this frame will be added.
    """
    # Build Settings
    group_settings_frame = ttk.LabelFrame(parent, text=" Groups ")
    group_settings_frame.pack(fill="both", padx=10, pady=5)

    # Add tabs for groups
    group_notebook = ttk.Notebook(group_settings_frame)
    group_notebook.pack(fill="both", expand=True)

    # Default: If no file is loaded
    if not hasattr(root.sys, 'groups') or root.sys.groups is None:
        no_groups_frame = ttk.Frame(group_notebook)
        group_notebook.add(no_groups_frame, text="No Groups")
        tk.Label(no_groups_frame, text="No file is loaded or no groups available.", font=root.fonts['class 2']).pack(
            padx=20, pady=20)
        return

    # Default: There are no groups
    if len(root.sys.groups) == 0:
        no_groups_frame = ttk.Frame(group_notebook)
        group_notebook.add(no_groups_frame, text="No Groups")
        tk.Label(no_groups_frame, text="System Must Be Loaded", font=root.fonts['class 2']).pack(
            padx=20, pady=20)
        return

    # Create tabs for each group
    for group in root.sys.groups:
        group_frame = ttk.Frame(group_notebook)
        group_notebook.add(group_frame, text=group.name)

        # Network Type Dropdown
        tk.Label(group_frame, text="Network Type", font=root.fonts['class 2']).grid(row=0, column=0, sticky="w")
        network_type = ttk.Combobox(group_frame, values=["Additively Weighted", "Power", "Primitive"])
        network_type.set("Additively Weighted")
        network_type.grid(row=0, column=1, sticky="w")

        # Maximum Vertex
        tk.Label(group_frame, text="Maximum Vertex", font=root.fonts['class 2']).grid(row=1, column=0, sticky="w")
        max_vertex = tk.Entry(group_frame)
        max_vertex.insert(0, "15")
        max_vertex.grid(row=1, column=1, sticky="w")

        # Box Multiplier
        tk.Label(group_frame, text="Box Multiplier", font=root.fonts['class 2']).grid(row=2, column=0, sticky="w")
        box_multiplier = tk.Entry(group_frame)
        box_multiplier.insert(0, "1.25")
        box_multiplier.grid(row=2, column=1, sticky="w")

        # Atomic Radii/Masses Button
        tk.Label(group_frame, text="Atomic Radii/Masses", font=root.fonts['class 2']).grid(row=3, column=0, sticky="w")
        tk.Button(group_frame, text="Open", command=root.open_atomic_radii_gui).grid(row=3, column=1, sticky="w")

        # Surface Settings Button
        tk.Label(group_frame, text="Surface Settings", font=root.fonts['class 2']).grid(row=4, column=0, sticky="w")
        tk.Button(group_frame, text="Change", command=root.open_surface_settings_gui).grid(row=4, column=1, sticky="w")
