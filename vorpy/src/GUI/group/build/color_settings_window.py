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
        self.geometry("300x280")  # Adjusted size to fit content (increased for checkbox)
        self.resizable(False, False)
        self.build_frame = build_frame

        # Make window modal
        self.transient(build_frame)
        self.grab_set()

        # Set up the main settings
        color_settings = self.build_frame.settings['color_settings']
        self.surf_res = tk.StringVar(value=color_settings['surf_res'])
        self.surf_col = tk.StringVar(value=color_settings['surf_col'])
        self.surf_scheme = tk.StringVar(value=color_settings['surf_scheme'].capitalize())
        self.surf_fact = tk.StringVar(value=color_settings['surf_fact'].capitalize())
        self.vert_col = tk.StringVar(value=color_settings['vert_col'].capitalize())
        self.edge_col = tk.StringVar(value=color_settings['edge_col'].capitalize())
        # Add concave_colors variable
        self.concave_colors = tk.BooleanVar(value=color_settings.get('conc_col'))

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

        # Surface resolution variable
        ttk.Label(settings_frame, text="Surface Resolution").grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.surf_res = ttk.Entry(settings_frame, width=15)
        self.surf_res.insert(0, default_settings['surf_res'])
        self.surf_res.grid(row=0, column=1, sticky="w", padx=5, pady=2)

        # Create dropdown for surface scheme with translations
        ttk.Label(settings_frame, text="Coloring Scheme").grid(row=1, column=0, sticky="w", padx=5, pady=2)
        self.surf_scheme = ttk.Combobox(
            settings_frame,
            values=[
                'Mean Curvature',
                'Gaussian Curvature',
                'Integrated Mean Curvature',
                'Integrated Mean Curvature Squared',
                'Integrated Gaussian Curvature',
                'Distance',
                'Overlapping',
                'No Scheme',
                'Average Mean Curvature',
                'Average Gaussian Curvature',
                'Maximum Mean Curvature',
                'Maximum Gaussian Curvature'
            ],
            state="readonly",
            width=32
        )
        # Set the initial value based on the current setting
        self.scheme_translations = {
            'mean': 'mean',
            'mean_curv': 'mean',
            'mean curvature': 'mean',
            'gauss': 'gauss',
            'gaus_curv': 'gauss',
            'gaussian curvature': 'gauss',
            'dist': 'dist',
            'distance': 'dist',
            'olap': 'olap',
            'overlapping': 'olap',
            'none': 'none',
            'no scheme': 'none',
            'avg_mean': 'avg_mean',
            'average mean curvature': 'avg_mean',
            'avg_gauss': 'avg_gauss',
            'average gaussian curvature': 'avg_gauss',
            'max_mean': 'max_mean',
            'maximum mean curvature': 'max_mean',
            'max_gauss': 'max_gauss',
            'maximum gaussian curvature': 'max_gauss',

            # Integrated curvature coloring
            'int_mean_curv': 'int_mean_curv',
            'int_mean': 'int_mean_curv',
            'integrated mean curvature': 'int_mean_curv',

            'int_mean_curv_sq': 'int_mean_curv_sq',
            'int_mean_sq': 'int_mean_curv_sq',
            'integrated mean curvature squared': 'int_mean_curv_sq',

            'int_gauss_curv': 'int_gauss_curv',
            'int_gauss': 'int_gauss_curv',
            'integrated gaussian curvature': 'int_gauss_curv',

            'representative surface energy': 'surf_energy',
            'surface energy': 'surf_energy',
            'surf_energy': 'surf_energy',
        }
        current_value = default_settings['surf_scheme']

        scheme_display = {
            'mean': 'Mean Curvature',
            'gauss': 'Gaussian Curvature',
            'dist': 'Distance',
            'olap': 'Overlapping',
            'ins_out': 'Overlapping',
            'none': 'No Scheme',
            'avg_mean': 'Average Mean Curvature',
            'avg_gauss': 'Average Gaussian Curvature',
            'max_mean': 'Maximum Mean Curvature',
            'max_gauss': 'Maximum Gaussian Curvature',
            'int_mean_curv': 'Integrated Mean Curvature',
            'int_mean_curv_sq': 'Integrated Mean Curvature Squared',
            'int_gauss_curv': 'Integrated Gaussian Curvature',
            'surf_energy': 'Representative Surface Energy',
        }

        self.surf_scheme.set(
            scheme_display.get(current_value.lower(), current_value)
        )
        self.surf_scheme.grid(row=1, column=1, sticky="w", padx=5, pady=2)
        self.surf_scheme.bind('<<ComboboxSelected>>')

        # Create the surface colorway entry
        ttk.Label(settings_frame, text="Surface Colorway").grid(row=2, column=0, sticky="w", padx=5, pady=2)
        self.surf_col = ttk.Entry(settings_frame, width=15)
        self.surf_col.insert(0, default_settings['surf_col'])
        self.surf_col.grid(row=2, column=1, sticky="w", padx=5, pady=2)

        # Create settings rows
        ttk.Label(settings_frame, text="Surface Coloring Factor").grid(row=3, column=0, sticky="w", padx=5, pady=2)
        # Create dropdown for surface factor
        self.surf_fact = ttk.Combobox(settings_frame, values=['Log', 'Linear', 'Exponential', 'Squared', 'Cubed'],
                                      state="readonly", width=15)
        self.surf_fact.set(default_settings['surf_fact'].capitalize())
        self.surf_fact.grid(row=3, column=1, sticky="w", padx=5, pady=2)
        self.surf_fact.bind('<<ComboboxSelected>>')

        # Add Concave Colors checkbox
        ttk.Label(settings_frame, text="Concave Colors").grid(row=4, column=0, sticky="w", padx=5, pady=2)
        self.concave_colors_checkbox = ttk.Checkbutton(
            settings_frame, variable=self.concave_colors
        )
        self.concave_colors_checkbox.grid(row=4, column=1, columnspan=2, sticky="w", padx=5, pady=2)

        # Create the vertex colorway entry
        ttk.Label(settings_frame, text="Vertex Color").grid(row=5, column=0, sticky="w", padx=5, pady=2)
        self.vert_col = ttk.Entry(settings_frame, width=15)
        self.vert_col.insert(0, default_settings['vert_col'].capitalize())
        self.vert_col.grid(row=5, column=1, sticky="w", padx=5, pady=2)

        # Create the edge colorway entry
        ttk.Label(settings_frame, text="Edge Color").grid(row=6, column=0, sticky="w", padx=5, pady=2)
        self.edge_col = ttk.Entry(settings_frame, width=15)
        self.edge_col.insert(0, default_settings['edge_col'].capitalize())
        self.edge_col.grid(row=6, column=1, sticky="w", padx=5, pady=2)

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
            'surf_scheme': self.scheme_translations[self.surf_scheme.get().lower()],
            'surf_fact': self.surf_fact.get().lower(),
            'conc_col': self.concave_colors.get(),
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
