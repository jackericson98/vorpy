import tkinter as tk
from tkinter import ttk


class CustomExportWindow(tk.Toplevel):
    """Window for custom export settings."""

    def __init__(self, parent, group_name):
        import os
        super().__init__(parent)
        self.title(f"{group_name} Export Settings")
        self.geometry("360x380")
        self.resizable(False, False)
        if not isinstance(parent, tk.Tk):
            self.transient(parent)
            self.grab_set()

        # Add icon in the same way as main frame
        try:
            vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
            icon_path = os.path.join(vorpy_root, "assets", "VorpyIcon.ico")
            if os.path.exists(icon_path):
                self.iconbitmap(icon_path)
        except Exception:
            pass

        self.parent = parent
        self.group_name = group_name

        # Main frame
        main_frame = ttk.Frame(self)
        main_frame.pack(fill="both", expand=True, padx=5, pady=5)

        # Make the group header
        group_header = ttk.Label(main_frame, text=f"{group_name}", font=("Helvetica", 15, "bold"))
        group_header.pack(padx=2, pady=(0, 5))

        # Info Frame
        info_frame = ttk.LabelFrame(main_frame, text="Info")
        info_frame.pack(fill="x", padx=2, pady=(0, 5))

        # Info checkbuttons in horizontal layout
        self.logs_var = tk.BooleanVar(value=True)
        self.verts_var = tk.BooleanVar(value=False)
        self.info_var = tk.BooleanVar(value=False)

        ttk.Checkbutton(info_frame, text="Logs", variable=self.logs_var).pack(side="left", padx=10, pady=2)
        ttk.Checkbutton(info_frame, text="Verts", variable=self.verts_var).pack(side="left", padx=10, pady=2)
        ttk.Checkbutton(info_frame, text="Info", variable=self.info_var).pack(side="left", padx=10, pady=2)

        # Balls Frame
        balls_frame = ttk.LabelFrame(main_frame, text="Balls")
        balls_frame.pack(fill="x", padx=2, pady=(0, 5))

        # Grid for balls
        grid_frame = ttk.Frame(balls_frame)
        grid_frame.pack(fill="x", padx=5, pady=2)

        # Column headers
        ttk.Label(grid_frame, text="Group").grid(row=1, column=0, padx=5, sticky="w")
        ttk.Label(grid_frame, text="Surface Neighbors").grid(row=2, column=0, padx=5, sticky="w")
        ttk.Label(grid_frame, text="Surrounding Residues").grid(row=3, column=0, padx=5, sticky="w")

        # File format options
        self.formats = ['pdb', 'cif', 'mol', 'gro', 'xyz', 'txt']

        # Create variables for checkbuttons
        self.group_vars = {}
        self.surrounding_vars = {}
        self.surrounding_resid_vars = {}

        # Create grid of checkbuttons with format labels
        for i, fmt in enumerate(self.formats, start=1):
            # Format label
            ttk.Label(grid_frame, text=fmt).grid(row=0, column=i, padx=(0, 5))

            # Group balls checkbutton
            self.group_vars[fmt] = tk.BooleanVar(value=fmt == 'pdb')
            ttk.Checkbutton(grid_frame, variable=self.group_vars[fmt]).grid(row=1, column=i, padx=5)

            # Surrounding balls checkbutton
            self.surrounding_vars[fmt] = tk.BooleanVar(value=False)
            ttk.Checkbutton(grid_frame, variable=self.surrounding_vars[fmt]).grid(row=2, column=i, padx=5)

            # Residue-expanded surrounding atoms
            self.surrounding_resid_vars[fmt] = tk.BooleanVar(value=False)
            ttk.Checkbutton(grid_frame, variable=self.surrounding_resid_vars[fmt]).grid(row=3, column=i, padx=5)

        # Create a frame for the three similar sections
        three_sections_frame = ttk.Frame(main_frame)
        three_sections_frame.pack(fill="x", padx=2, pady=(0, 5))

        # Surfaces Frame
        surfs_frame = ttk.LabelFrame(three_sections_frame, text="Surfaces")
        surfs_frame.pack(side="left", fill="x", expand=True, padx=(0, 2))

        # Edges Frame
        edges_frame = ttk.LabelFrame(three_sections_frame, text="Edges")
        edges_frame.pack(side="left", fill="x", expand=True, padx=2)

        # Vertices Frame
        verts_frame = ttk.LabelFrame(three_sections_frame, text="Vertices")
        verts_frame.pack(side="left", fill="x", expand=True, padx=(2, 0))

        # Function to create checkbuttons in a frame
        def create_section_checkbuttons(frame, all_var, shell_var, separate_var, cell_var):
            ttk.Checkbutton(frame, text="All", variable=all_var).pack(anchor="w", padx=5, pady=1)
            ttk.Checkbutton(frame, text="Shell", variable=shell_var).pack(anchor="w", padx=5, pady=1)
            ttk.Checkbutton(frame, text="Separate", variable=separate_var).pack(anchor="w", padx=5, pady=1)
            ttk.Checkbutton(frame, text="Cell", variable=cell_var).pack(anchor="w", padx=5, pady=1)

        # Create checkbuttons for each section
        self.surfs_all_var = tk.BooleanVar(value=False)
        self.surfs_shell_var = tk.BooleanVar(value=True)
        self.surfs_separate_var = tk.BooleanVar(value=False)
        self.surfs_cell_var = tk.BooleanVar(value=False)
        create_section_checkbuttons(surfs_frame, self.surfs_all_var, self.surfs_shell_var, self.surfs_separate_var,
                                    self.surfs_cell_var)

        self.edges_all_var = tk.BooleanVar(value=False)
        self.edges_shell_var = tk.BooleanVar(value=True)
        self.edges_separate_var = tk.BooleanVar(value=False)
        self.edges_cell_var = tk.BooleanVar(value=False)
        create_section_checkbuttons(edges_frame, self.edges_all_var, self.edges_shell_var, self.edges_separate_var,
                                    self.edges_cell_var)

        self.verts_all_var = tk.BooleanVar(value=False)
        self.verts_shell_var = tk.BooleanVar(value=False)
        self.verts_separate_var = tk.BooleanVar(value=False)
        self.verts_cell_var = tk.BooleanVar(value=False)
        create_section_checkbuttons(verts_frame, self.verts_all_var, self.verts_shell_var, self.verts_separate_var,
                                    self.verts_cell_var)

        # Buttons frame
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill="x", pady=(5, 0))

        ttk.Button(button_frame, text="Apply", command=self._on_ok).pack(side="right", padx=(5, 0), pady=5)
        ttk.Button(button_frame, text="Cancel", command=self._on_cancel).pack(side="right", padx=(5, 0), pady=5)

        # Center the window on the parent
        self.update_idletasks()
        width = self.winfo_width()
        height = self.winfo_height()
        x = (self.winfo_screenwidth() // 2) - (width // 2)
        y = (self.winfo_screenheight() // 2) - (height // 2)
        self.geometry(f"{width}x{height}+{x}+{y}")

    def _on_ok(self):
        """Handle OK button click."""
        self.settings = {
            'logs': self.logs_var.get(),
            'verts': self.verts_var.get(),
            'info': self.info_var.get(),
            'surfs_all': self.surfs_all_var.get(),
            'surfs_shell': self.surfs_shell_var.get(),
            'surfs_separate': self.surfs_separate_var.get(),
            'surfs_cell': self.surfs_cell_var.get(),
            'edges_all': self.edges_all_var.get(),
            'edges_shell': self.edges_shell_var.get(),
            'edges_separate': self.edges_separate_var.get(),
            'edges_cell': self.edges_cell_var.get(),
            'verts_all': self.verts_all_var.get(),
            'verts_shell': self.verts_shell_var.get(),
            'verts_separate': self.verts_separate_var.get(),
            'verts_cell': self.verts_cell_var.get(),
            'group_vars': {fmt: self.group_vars[fmt].get() for fmt in self.formats},
            'surrounding_vars': {fmt: self.surrounding_vars[fmt].get() for fmt in self.formats},
            'surrounding_resid_vars': {
                fmt: self.surrounding_resid_vars[fmt].get() for fmt in self.formats
            },
        }
        # Update the export frame's settings directly
        self.parent.settings['custom_settings'] = self.settings
        self.destroy()

    def _on_cancel(self):
        """Handle Cancel button click."""
        self.destroy()
