import sys
from pathlib import Path

# Add the project root directory to the Python path
project_root = Path(__file__).resolve().parents[4]
sys.path.append(str(project_root))

import tkinter as tk
from tkinter import ttk
from vorpy.src.GUI.group.build.color_settings_window import ColorSettingsWindow


class BuildFrame(ttk.LabelFrame):
    """
    A frame for build settings configuration.
    """
    def __init__(self, parent, gui):
        super().__init__(parent, text="Build Settings")
        self.gui = gui
        
        # Initialize default settings for this group
        self.settings = {
            'max_vert': 40,
            'box_size': 1.25,
            'net_type': 'Additively Weighted',
            'color_settings': {
                'surf_res': '0.2',
                'surf_col': 'plasma',
                'surf_scheme': 'Mean Curvature',
                'surf_fact': 'Linear',
                'vert_col': 'red',
                'edge_col': 'grey',
                'conc_col': True
            }
        }
        self.color_settings = None
        
        # Create widgets after settings are initialized
        self._create_widgets()
        
    def _create_widgets(self):
        """Create and pack all widgets in the frame."""
        # Configure grid weights for the frame
        self.grid_columnconfigure(1, weight=1)  # Make the middle column expand
        
        # Network Type
        network_frame = ttk.Frame(self)
        network_frame.pack(fill="x", padx=5, pady=5)
        
        ttk.Label(network_frame, text="Network Type:").pack(side="left")
        self.network_type = ttk.Combobox(network_frame, values=['Additively Weighted', 'Power', 'Primitive'],
                                         state="readonly")
        self.network_type.pack(side="right", padx=5)
        self.network_type.set(self.settings['net_type'])
        self.network_type.bind('<<ComboboxSelected>>', self._on_network_type_change)
        
        # Max Vertices Entry
        vertices_frame = ttk.Frame(self)
        vertices_frame.pack(fill="x", padx=5, pady=5)
        
        ttk.Label(vertices_frame, text="Probe Distance:").pack(side="left")
        self.max_vertices = ttk.Entry(vertices_frame, width=10)
        self.max_vertices.pack(side="right", padx=5)
        self.max_vertices.insert(0, str(self.settings['max_vert']))
        self.max_vertices.bind('<KeyRelease>', self._on_max_vertices_change)
        
        # Box Size Entry
        box_frame = ttk.Frame(self)
        box_frame.pack(fill="x", padx=5, pady=5)
        
        ttk.Label(box_frame, text="Box Size:").pack(side="left")
        self.box_size = ttk.Entry(box_frame, width=10)
        self.box_size.pack(side="right", padx=5)
        self.box_size.insert(0, str(self.settings['box_size']))
        self.box_size.bind('<KeyRelease>', self._on_box_size_change)
        
        # Surface Settings Button
        surface_button = ttk.Button(self, text="Outputs Color Settings", command=self._open_surface_settings)
        surface_button.pack(pady=5)
        
    def _on_network_type_change(self, event=None):
        self.settings['net_type'] = self.network_type.get()
    
    def _on_max_vertices_change(self, event=None):
        try:
            self.settings['max_vert'] = int(self.max_vertices.get())
        except ValueError:
            pass
    
    def _on_box_size_change(self, event=None):
        try:
            self.settings['box_size'] = float(self.box_size.get())
        except ValueError:
            pass
    
    def _open_surface_settings(self):
        """Open the surface settings window."""
        self.color_settings = ColorSettingsWindow(self, self.settings['color_settings'])
    
    def get_settings(self):
        """Return the current build settings."""
        if self.color_settings is not None:
            color_settings = self.color_settings.get_settings()
        else:
            color_settings = self.settings['color_settings']
        self.settings = {
            'net_type': self.network_type.get(),
            'max_vert': self.max_vertices.get(),
            'box_size': self.box_size.get(),
            'color_settings': color_settings
        }
        return self.settings
    
    def copy_settings_from(self, other_frame):
        """Copy settings from another build frame."""
        # Copy the settings
        self.settings = other_frame.get_settings()
        
        # Update the UI to reflect the new settings
        self.network_type.set(self.settings['net_type'])
        self.max_vertices.delete(0, tk.END)
        self.max_vertices.insert(0, str(self.settings['max_vert']))
        self.box_size.delete(0, tk.END)
        self.box_size.insert(0, str(self.settings['box_size']))


if __name__ == "__main__":
    root = tk.Tk()
    build_frame = BuildFrame(root, None)
    build_frame.pack(fill="both", expand=True)
    root.mainloop()

    
