import os
import sys
from pathlib import Path

# Add the project root directory to the Python path
project_root = Path(__file__).resolve().parents[4]
sys.path.append(str(project_root))

import tkinter as tk
from tkinter import ttk
from Visualize.GUI.settings.surface.surface_settings_window import SurfaceSettingsWindow


class BuildFrame(ttk.LabelFrame):
    """
    A frame for build settings configuration.
    """
    def __init__(self, parent, gui):
        super().__init__(parent, text="Build Settings", padding="10")
        self.gui = gui
        
        # Create and pack widgets
        self._create_widgets()
        
    def _create_widgets(self):
        """Create and pack all widgets in the frame."""
        # Max Vert Rad
        ttk.Label(self, text="Max Vert Rad:").grid(row=0, column=0, sticky="e", padx=5, pady=5)
        self.max_vert_rad = ttk.Entry(self, width=10)
        self.max_vert_rad.grid(row=0, column=1, sticky="w", padx=5, pady=5)
        self.max_vert_rad.insert(0, "0.01")
        
        # Max Box Multi
        ttk.Label(self, text="Max Box Multi:").grid(row=1, column=0, sticky="e", padx=5, pady=5)
        self.max_box_multi = ttk.Entry(self, width=10)
        self.max_box_multi.grid(row=1, column=1, sticky="w", padx=5, pady=5)
        self.max_box_multi.insert(0, "1")
        
        # Network Type
        ttk.Label(self, text="Network Type:").grid(row=2, column=0, sticky="e", padx=5, pady=5)
        self.network_type = ttk.Combobox(self, values=["Additively Weighted", "Power Diagram", "Primitive (Delaunay)"], 
                                       state="readonly", width=20)
        self.network_type.grid(row=2, column=1, sticky="w", padx=5, pady=5)
        self.network_type.set("Additively Weighted")
        
        # Surface Settings
        ttk.Label(self, text="Surface Settings").grid(row=3, column=0, columnspan=2, sticky="w", padx=5, pady=(10, 5))
        
        # Surface Values
        self.surface_values_label = ttk.Label(self, text="Setting 1: 0.0\nSetting 2: 0.0\nSetting 3: 0.0\nSetting 4: 0.0")
        self.surface_values_label.grid(row=4, column=0, columnspan=3, sticky="w", padx=5, pady=5)
        
        # Change Button
        change_button = ttk.Button(self, text="Change", command=self.open_surface_settings_gui)
        change_button.grid(row=5, column=0, columnspan=2, pady=10)
        
    def open_surface_settings_gui(self):
        """Open the surface settings window."""
        SurfaceSettingsWindow(self.gui)


if __name__ == "__main__":
    root = tk.Tk()
    build_frame = BuildFrame(root, None)
    build_frame.pack(fill="both", expand=True)
    root.mainloop()

    
