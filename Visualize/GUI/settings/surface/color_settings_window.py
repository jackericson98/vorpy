import tkinter as tk
from tkinter import ttk


class ColorSettingsWindow(tk.Toplevel):
    """
    A window for configuring surface settings.
    """
    def __init__(self, parent):
        super().__init__(parent)
        
        # Configure window
        self.title("Color Settings")
        self.geometry("300x220")  # Adjusted size to fit content
        self.resizable(False, False)
        self.gui = parent
        
        # Make window modal
        self.transient(parent)
        self.grab_set()

        # Set up the main settings
        self.surf_col = tk.StringVar(value=self.gui.build_settings['color_settings']['surf_col'].capitalize())
        self.surf_scheme = tk.StringVar(value=self.gui.build_settings['color_settings']['surf_scheme'].capitalize())
        self.surf_fact = tk.StringVar(value=self.gui.build_settings['color_settings']['surf_fact'].capitalize())
        self.vert_col = tk.StringVar(value=self.gui.build_settings['color_settings']['vert_col'].capitalize())
        self.edge_col = tk.StringVar(value=self.gui.build_settings['color_settings']['edge_col'].capitalize())
        
        # Create main frame with proper padding
        main_frame = ttk.Frame(self, padding="10")
        main_frame.pack(fill="both", expand=True)
        
        # Surface Settings
        settings_frame = ttk.LabelFrame(main_frame, text="Color Settings", padding="5")
        settings_frame.pack(fill="both", expand=True, padx=5, pady=5)
        
        # Configure grid weights for better layout
        settings_frame.grid_columnconfigure(1, weight=1)
        
        # Create settings rows
        self.settings = {}
        self._create_setting_row(settings_frame, "Surface Colorway", 0, self.surf_col)
        self._create_setting_row(settings_frame, "Surface Coloring Scheme", 1, self.surf_scheme)
        self._create_setting_row(settings_frame, "Surface Coloring Factor", 2, self.surf_fact)
        self._create_setting_row(settings_frame, "Vertex Color", 3, self.vert_col)
        self._create_setting_row(settings_frame, "Edge Color", 4, self.edge_col)
        
        # Buttons frame with proper spacing
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill="x", pady=(10, 0))
        
        # OK and Cancel buttons
        ttk.Button(button_frame, text="OK", command=self._on_ok).pack(side="right", padx=5)
        ttk.Button(button_frame, text="Cancel", command=self._on_cancel).pack(side="right", padx=5)
        
        # Center the window on the parent
        self.update_idletasks()
        width = self.winfo_width()
        height = self.winfo_height()
        x = (self.winfo_screenwidth() // 2) - (width // 2)
        y = (self.winfo_screenheight() // 2) - (height // 2)
        self.geometry(f"{width}x{height}+{x}+{y}")

    def _create_setting_row(self, parent, label_text, row, value, command=None):
        """Create a row with a label and input widget."""
        ttk.Label(parent, text=label_text).grid(row=row, column=0, sticky="w", padx=5, pady=2)
        
        if label_text == "Surface Coloring Factor":
            # Create dropdown for surface factor
            self.surf_fact = ttk.Combobox(parent, values=['Log', 'Linear', 'Exponential', 'Squared', 'Cubed'], 
                                        state="readonly", width=15)
            self.surf_fact.set(value.get())
            self.surf_fact.grid(row=row, column=1, sticky="w", padx=5, pady=2)
            if command:
                self.surf_fact.bind('<<ComboboxSelected>>', command)
        elif label_text == "Surface Coloring Scheme":
            # Create dropdown for surface scheme with translations
            self.surf_scheme = ttk.Combobox(parent, 
                                          values=['Mean Curvature', 'Gaussian Curvature', 'Distance', 'Overlapping', 'No Scheme'],
                                          state="readonly", width=15)
            # Set the initial value based on the current setting
            scheme_translations = {
                'mean_curv': 'Mean Curvature',
                'gaus_curv': 'Gaussian Curvature',
                'dist': 'Distance',
                'olap': 'Overlapping',
                'none': 'No Scheme'
            }
            current_value = value.get()
            display_value = scheme_translations.get(current_value, current_value)
            self.surf_scheme.set(display_value)
            self.surf_scheme.grid(row=row, column=1, sticky="w", padx=5, pady=2)
            if command:
                self.surf_scheme.bind('<<ComboboxSelected>>', command)
        else:
            # Create entry for other settings
            entry = ttk.Entry(parent, width=15)
            entry.insert(0, value.get())
            entry.grid(row=row, column=1, sticky="w", padx=5, pady=2)
            if command:
                entry.bind('<KeyRelease>', command)
            return entry

    def _on_ok(self):
        """Handle OK button click."""
        # Update the settings in the main GUI
        self.gui.build_settings['color_settings']['surf_col'] = self.surf_col.get().lower()
        
        # Translate the surface scheme value back to the internal format
        scheme_translations = {
            'Mean Curvature': 'mean_curv',
            'Gaussian Curvature': 'gaus_curv',
            'Distance': 'dist',
            'Overlapping': 'olap',
            'No Scheme': 'none'
        }
        self.gui.build_settings['color_settings']['surf_scheme'] = scheme_translations.get(self.surf_scheme.get(), self.surf_scheme.get())
        
        self.gui.build_settings['color_settings']['surf_fact'] = self.surf_fact.get().lower()
        self.gui.build_settings['color_settings']['vert_col'] = self.vert_col.get().lower()
        self.gui.build_settings['color_settings']['edge_col'] = self.edge_col.get().lower()
        
        # Update the display in the main GUI
        self.gui.update_surface_settings_display()
        self.destroy()

    def _on_cancel(self):
        """Handle Cancel button click."""
        self.destroy() 