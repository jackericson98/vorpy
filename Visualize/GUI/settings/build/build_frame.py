import os
import sys
from pathlib import Path

# Add the project root directory to the Python path
project_root = Path(__file__).resolve().parents[4]
sys.path.append(str(project_root))

import tkinter as tk
from tkinter import ttk
from Visualize.GUI.settings.surface.color_settings_window import ColorSettingsWindow


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
        # Configure grid weights for the frame
        self.grid_columnconfigure(1, weight=1)  # Make the middle column expand
        
        # Network Type
        ttk.Label(self, text="Network Type:").grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.network_type = ttk.Combobox(self, values=["Delaunay", "Gabriel", "Relative Neighborhood", "Beta Skeleton"], 
                                        state="readonly", width=14, justify="right")
        self.network_type.set(self.gui.build_settings['net_type'])
        self.network_type.grid(row=0, column=1, columnspan=1, sticky="e", padx=5, pady=2)
        self.network_type.bind('<<ComboboxSelected>>', self._update_net_type)
        
        # Probe Radius
        ttk.Label(self, text="Probe Radius:").grid(row=1, column=0, sticky="w", padx=5, pady=2)
        self.max_vert_rad = ttk.Entry(self, width=17, justify="right")
        self.max_vert_rad.insert(0, self.gui.build_settings['max_vert'])
        self.max_vert_rad.grid(row=1, column=1, sticky="e", padx=5, pady=2)
        self.max_vert_rad.bind('<KeyRelease>', self._update_max_vert)
        ttk.Label(self, text=u"\u212b").grid(row=1, column=2, sticky="e", padx=5, pady=2)
        
        # Outer Reach
        ttk.Label(self, text="Outer Reach:").grid(row=2, column=0, sticky="w", padx=5, pady=2)
        self.max_box_multi = ttk.Entry(self, width=17, justify="right")
        self.max_box_multi.insert(0, self.gui.build_settings['box_size'])
        self.max_box_multi.grid(row=2, column=1, sticky="e", padx=5, pady=2)
        self.max_box_multi.bind('<KeyRelease>', self._update_max_box)
        ttk.Label(self, text="x ").grid(row=2, column=2, sticky="e", padx=5, pady=2)
        
        # Surface Settings
        ttk.Label(self, text="Color Settings").grid(row=3, column=0, sticky="w", padx=5, pady=5)
        
        # Change Button
        change_button = ttk.Button(self, text="Change", command=self.open_surface_settings_gui)
        change_button.grid(row=3, column=1, columnspan=1, sticky="e", padx=5, pady=5)
        
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
        ColorSettingsWindow(self.gui)

    def update_surface_settings_display(self):
        """Update the display of surface settings."""
        # Update the surface settings values
        self.surface_values_label.config(
            text=f"Surface Color: {self.gui.build_settings['color_settings']['surf_col']}\n"
                 f"Surface Scheme: {self.gui.build_settings['color_settings']['surf_scheme']}\n"
                 f"Scheme Factor: {self.gui.build_settings['color_settings']['surf_fact']}\n"
                 f"Edge Color: {self.gui.build_settings['color_settings']['edge_col']}\n"
                 f"Vertex Color: {self.gui.build_settings['color_settings']['vert_col']}"
        )


if __name__ == "__main__":
    root = tk.Tk()
    build_frame = BuildFrame(root, None)
    build_frame.pack(fill="both", expand=True)
    root.mainloop()

    
