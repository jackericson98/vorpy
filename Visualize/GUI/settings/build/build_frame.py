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
        self.settings = {'surf_res': 0.2, 'surf_col': 'plasma', 'surf_scheme': 'mean curvature', 'max_vert': 40, 'box_size': 1.25, 'net_type': 'additively weighted', 
                         'build_type': 'all', 'num_splits': 0, 'print_metrics': True, 'ball_type': 'all', 'foam_box': False,'atom_rad': None, 'scheme_factor': 'log',
                         'edge_col': 'grey', 'vert_col': 'red'}
        
        # Create and pack widgets
        self._create_widgets()
        
    def _create_widgets(self):
        """Create and pack all widgets in the frame."""
        # Max Vert Rad
        ttk.Label(self, text="Probe Radius").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        self.max_vert_rad = ttk.Entry(self, width=15)
        self.max_vert_rad.grid(row=0, column=1, sticky="e", padx=5, pady=5)
        self.max_vert_rad.insert(0, str(self.settings['max_vert']))
        self.max_vert_rad.bind('<KeyRelease>', self._update_max_vert)
        
        # Max Box Multi
        ttk.Label(self, text="Outer Reach").grid(row=1, column=0, sticky="w", padx=5, pady=5)
        self.max_box_multi = ttk.Entry(self, width=15)
        self.max_box_multi.grid(row=1, column=1, sticky="e", padx=5, pady=5)
        self.max_box_multi.insert(0, str(self.settings['box_size']))
        self.max_box_multi.bind('<KeyRelease>', self._update_max_box)
        
        # Network Type
        ttk.Label(self, text="Network Type:").grid(row=2, column=0, sticky="w", padx=5, pady=5)
        self.network_type = ttk.Combobox(self, values=["1. Additively Weighted", "2. Power Diagram", "3. Primitive (Delaunay)", "4. Combo (1 & 2)", "5. Combo (1 & 3)", "6. Combo (2 & 3)", "7. All 3"], 
                                       state="readonly", width=15)
        self.network_type.grid(row=2, column=1, sticky="e", padx=5, pady=5)
        self.network_type.set(self.settings['net_type'])
        self.network_type.bind('<<ComboboxSelected>>', self._update_net_type)
        
        # Surface Settings
        ttk.Label(self, text="Color Settings").grid(row=3, column=0, columnspan=1, sticky="w", padx=5, pady=5)
        
        # # Surface Values
        # # Surface Color
        # ttk.Label(self, text="Surf. Color:").grid(row=4, column=0, sticky="e", padx=5, pady=2)
        # ttk.Label(self, text=self.settings['surf_col']).grid(row=4, column=1, sticky="w", padx=5, pady=2)
        
        # # Surface Scheme
        # ttk.Label(self, text="Surf. Scheme:").grid(row=5, column=0, sticky="e", padx=5, pady=2)
        # ttk.Label(self, text=self.settings['surf_scheme']).grid(row=5, column=1, sticky="w", padx=5, pady=2)
        
        # # Scheme Factor
        # ttk.Label(self, text="Surf. Scheme Factor:").grid(row=6, column=0, sticky="e", padx=5, pady=2)
        # ttk.Label(self, text=self.settings['scheme_factor']).grid(row=6, column=1, sticky="w", padx=5, pady=2)
        
        # # Edge Color
        # ttk.Label(self, text="Edge Color:").grid(row=7, column=0, sticky="e", padx=5, pady=2)
        # ttk.Label(self, text=self.settings['edge_col']).grid(row=7, column=1, sticky="w", padx=5, pady=2)
        
        # # Vertex Color
        # ttk.Label(self, text="Vertex Color:").grid(row=8, column=0, sticky="e", padx=5, pady=2)
        # ttk.Label(self, text=self.settings['vert_col']).grid(row=8, column=1, sticky="w", padx=5, pady=2)
        
        # Change Button
        change_button = ttk.Button(self, text="Change", command=self.open_surface_settings_gui)
        change_button.grid(row=3, column=1, columnspan=1, pady=5, sticky="e")
        
    def _update_max_vert(self, event=None):
        """Update max_vert setting when entry changes."""
        try:
            value = float(self.max_vert_rad.get())
            if 0.01 <= value <= 500:
                self.settings['max_vert'] = value
        except ValueError:
            pass  # Ignore invalid input
            
    def _update_max_box(self, event=None):
        """Update max_box setting when entry changes."""
        try:
            value = float(self.max_box_multi.get())
            if 1 <= value <= 200:
                self.settings['box_size'] = value
        except ValueError:
            pass  # Ignore invalid input
            
    def _update_net_type(self, event=None):
        """Update net_type setting when combobox selection changes."""
        self.settings['net_type'] = self.network_type.get()
        
    def open_surface_settings_gui(self):
        """Open the surface settings window."""
        SurfaceSettingsWindow(self.gui)


if __name__ == "__main__":
    root = tk.Tk()
    build_frame = BuildFrame(root, None)
    build_frame.pack(fill="both", expand=True)
    root.mainloop()

    
