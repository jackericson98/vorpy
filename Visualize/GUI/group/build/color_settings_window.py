import tkinter as tk
from tkinter import ttk


class ColorSettingsWindow(tk.Toplevel):
    """
    A window for configuring surface settings.
    """
    def __init__(self, build_frame, default_settings):
        super().__init__(build_frame)
        
        # Configure window
        self.title("Color Settings")
        self.geometry("300x250")  # Adjusted size to fit content
        self.resizable(False, False)
        self.build_frame = build_frame
        
        # Make window modal
        self.transient(build_frame)
        self.grab_set()

        # Set up the main settings
        color_settings = self.build_frame.settings['color_settings']
        self.surf_res = tk.StringVar(value=color_settings['surf_res'])
        self.surf_col = tk.StringVar(value=color_settings['surf_col'].capitalize())
        self.surf_scheme = tk.StringVar(value=color_settings['surf_scheme'].capitalize())
        self.surf_fact = tk.StringVar(value=color_settings['surf_fact'].capitalize())
        self.vert_col = tk.StringVar(value=color_settings['vert_col'].capitalize())
        self.edge_col = tk.StringVar(value=color_settings['edge_col'].capitalize())
        
        # Store the current settings
        self.current_settings = color_settings.copy()
        
        # Create main frame with proper padding
        main_frame = ttk.Frame(self, padding="10")
        main_frame.pack(fill="both", expand=True)
        
        # Surface Settings
        settings_frame = ttk.LabelFrame(main_frame, text="Color Settings", padding="5")
        settings_frame.pack(fill="both", expand=True, padx=5, pady=5)
        


        # Configure grid weights for better layout
        settings_frame.grid_columnconfigure(1, weight=1)
        

        ttk.Label(settings_frame, text="Surface Resolution").grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.surf_res = ttk.Entry(settings_frame, width=15)
        self.surf_res.insert(0, default_settings['surf_res'])
        self.surf_res.grid(row=0, column=1, sticky="w", padx=5, pady=2)

        # Create the surface colorway entry
        ttk.Label(settings_frame, text="Surface Colorway").grid(row=1, column=0, sticky="w", padx=5, pady=2)
        self.surf_col = ttk.Entry(settings_frame, width=15)
        self.surf_col.insert(0, default_settings['surf_col'].capitalize())
        self.surf_col.grid(row=1, column=1, sticky="w", padx=5, pady=2)

        # Create settings rows
        ttk.Label(settings_frame, text="Surface Coloring Factor").grid(row=2, column=0, sticky="w", padx=5, pady=2)
        # Create dropdown for surface factor
        self.surf_fact = ttk.Combobox(settings_frame, values=['Log', 'Linear', 'Exponential', 'Squared', 'Cubed'], 
                                    state="readonly", width=15)
        self.surf_fact.set(default_settings['surf_fact'].capitalize())
        self.surf_fact.grid(row=2, column=1, sticky="w", padx=5, pady=2)
        self.surf_fact.bind('<<ComboboxSelected>>')

        # Create dropdown for surface scheme with translations
        self.surf_scheme = ttk.Combobox(settings_frame, values=['Mean Curvature', 'Gaussian Curvature', 'Distance', 'Overlapping', 'No Scheme'], state="readonly", width=15)
        # Set the initial value based on the current setting
        scheme_translations = {
            'mean_curv': 'Mean Curvature',
            'mean curvature': 'Mean Curvature',
            'gaus_curv': 'Gaussian Curvature',
            'gaussian curvature': 'Gaussian Curvature',
            'dist': 'Distance',
            'distance': 'Distance',
            'olap': 'Overlapping',
            'overlapping': 'Overlapping',
            'none': 'No Scheme',
            'no scheme': 'No Scheme'
        }
        current_value = default_settings['surf_scheme']
        display_value = scheme_translations[current_value.lower()]
        self.surf_scheme.set(display_value)
        self.surf_scheme.grid(row=3, column=1, sticky="w", padx=5, pady=2)
        self.surf_scheme.bind('<<ComboboxSelected>>')

        # Create the vertex colorway entry
        ttk.Label(settings_frame, text="Vertex Color").grid(row=4, column=0, sticky="w", padx=5, pady=2)
        self.vert_col = ttk.Entry(settings_frame, width=15)
        self.vert_col.insert(0, default_settings['vert_col'].capitalize())
        self.vert_col.grid(row=4, column=1, sticky="w", padx=5, pady=2)

        # Create the edge colorway entry
        ttk.Label(settings_frame, text="Edge Color").grid(row=5, column=0, sticky="w", padx=5, pady=2)
        self.edge_col = ttk.Entry(settings_frame, width=15)
        self.edge_col.insert(0, default_settings['edge_col'].capitalize())
        self.edge_col.grid(row=5, column=1, sticky="w", padx=5, pady=2)
        
        # Buttons frame with proper spacing
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill="x", pady=(10, 0))
        
        # OK and Cancel buttons
        ttk.Button(button_frame, text="OK", command=self._on_ok).pack(side="right", padx=5)
        ttk.Button(button_frame, text="Cancel", command=self._on_cancel).pack(side="right", padx=5)

    def _on_ok(self):
        """Handle OK button click."""
        # Update the current settings with the widget values
        self.current_settings = {
            'surf_res': float(self.surf_res.get()),
            'surf_col': self.surf_col.get().lower(),
            'surf_scheme': self.surf_scheme.get(),
            'surf_fact': self.surf_fact.get().lower(),
            'vert_col': self.vert_col.get().lower(),
            'edge_col': self.edge_col.get().lower()
        }
        
        # Update the build frame settings
        self.build_frame.settings['color_settings'] = self.current_settings
        
        # Close the window
        self.destroy()

    def get_settings(self):
        """Return the current color settings."""
        return self.current_settings

    def _on_cancel(self):
        """Handle Cancel button click."""
        self.destroy()

if __name__ == "__main__":
    root = tk.Tk()
    ColorSettingsWindow(root)
    root.mainloop()

