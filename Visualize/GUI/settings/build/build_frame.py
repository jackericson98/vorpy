import tkinter as tk
from tkinter import ttk


def build_frame(root, parent):
    """
    Builds the frame for group configuration within the parent frame.

    Args:
        root: The main GUI application object, containing system and font settings.
        parent: The parent frame where this frame will be added.
    """
    # Build Settings
    build_settings_frame = ttk.LabelFrame(parent, text=" Build Settings ")
    build_settings_frame.pack(fill="both", padx=10, pady=5)

    # Top settings row
    max_vert_label = tk.Label(build_settings_frame, text="Max Vert Rad (0.01-500)")
    max_vert_label.grid(row=0, column=0, padx=5, pady=5, sticky="w")

    max_vert_entry = tk.Entry(build_settings_frame, width=10)
    max_vert_entry.grid(row=0, column=1, padx=5, pady=5, sticky="w")

    max_box_label = tk.Label(build_settings_frame, text="Max Box Multi (1-200)")
    max_box_label.grid(row=1, column=0, padx=5, pady=5, sticky="w")

    max_box_entry = tk.Entry(build_settings_frame, width=10)
    max_box_entry.grid(row=1, column=1, padx=5, pady=5, sticky="w")

    # Network Type row
    network_type_label = tk.Label(build_settings_frame, text="Network Type")
    network_type_label.grid(row=2, column=0, padx=5, pady=5, sticky="w")

    network_types = [
        "Additively Weighted",
        "Power Diagram",
        "Primitive (Delaunay)",
        "Compare (1 & 2)",
        "Compare (1 & 3)",
        "Compare (2 & 3)",
        "Compare (all)"
    ]
    
    network_type_var = tk.StringVar()
    network_type_var.set(network_types[0])  # Set default value
    
    network_type_combo = ttk.Combobox(build_settings_frame, 
                                    textvariable=network_type_var,
                                    values=network_types,
                                    state="readonly",
                                    width=20)
    network_type_combo.grid(row=2, column=1, columnspan=3, padx=5, pady=5, sticky="w")

    # Surface Settings row
    surface_label = tk.Label(build_settings_frame, text="Surface Settings")
    surface_label.grid(row=3, column=0, padx=5, pady=5, sticky="w")

    # Surface settings values (to be filled in later)
    surface_values = ["Setting 1", "Setting 2", "Setting 3", "Setting 4"]
    surface_text = ", ".join(surface_values)
    surface_values_label = tk.Label(build_settings_frame, text=surface_text)
    surface_values_label.grid(row=3, column=1, columnspan=3, padx=5, pady=5, sticky="w")

    # Button to open surface settings window
    surface_button = ttk.Button(build_settings_frame, 
                              text="Change",
                              command=root.open_surface_settings_gui)
    surface_button.grid(row=4, column=0, columnspan=3, padx=5, pady=5)

    # Subframes for future development
    subframe1 = ttk.LabelFrame(build_settings_frame, text=" Subframe 1 ")
    subframe1.grid(row=5, column=0, columnspan=4, padx=10, pady=10, sticky="nsew")

    subframe2 = ttk.LabelFrame(build_settings_frame, text=" Subframe 2 ")
    subframe2.grid(row=6, column=0, columnspan=4, padx=10, pady=10, sticky="nsew")

    # Allow subframes to stretch
    build_settings_frame.columnconfigure(0, weight=1)
    build_settings_frame.columnconfigure(1, weight=1)
    build_settings_frame.columnconfigure(2, weight=1)
    build_settings_frame.columnconfigure(3, weight=1)
    build_settings_frame.rowconfigure(2, weight=1)
    build_settings_frame.rowconfigure(3, weight=1)
