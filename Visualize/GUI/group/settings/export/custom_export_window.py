import tkinter as tk
from tkinter import ttk


class CustomExportWindow(tk.Toplevel):
    """Window for custom export settings."""

    def __init__(self, parent):
        super().__init__(parent)
        self.title("Custom Export Settings")
        self.geometry("600x400")  # Adjusted size for new layout
        self.resizable(False, False)
        self.transient(parent)
        self.grab_set()

        # Main frame
        main_frame = ttk.Frame(self)
        main_frame.pack(fill="both", expand=True, padx=10, pady=10)

        # System frame
        system_frame = ttk.LabelFrame(main_frame, text="System")
        system_frame.pack(fill="x", padx=5, pady=5)

        # Ball Files section
        ttk.Label(system_frame, text="Ball Files:", font=('Helvetica', 10, 'bold')).grid(column=0, columnspan=5, row=0,
                                                                                         padx=5, pady=5, sticky="e")

        self.pdb_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(system_frame, text="PDB", variable=self.pdb_var).grid(column=0, row=1, padx=5, pady=1)

        self.cif_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(system_frame, text="CIF", variable=self.cif_var).grid(column=1, row=1, padx=5, pady=1)

        self.mol_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(system_frame, text="MOL", variable=self.mol_var).grid(column=2, row=1, padx=5, pady=1)

        self.gro_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(system_frame, text="GRO", variable=self.gro_var).grid(column=3, row=1, padx=5, pady=1)

        self.txt_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(system_frame, text="TXT", variable=self.txt_var).grid(column=4, row=1, padx=5, pady=1)

        # Other Files section
        ttk.Label(system_frame, text="Other Files:", font=('Helvetica', 10, 'bold')).grid(column=0, columnspan=5, row=2,
                                                                                          padx=5, pady=5, sticky="e")

        self.pymol_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(system_frame, text="Set Balls PyMOL", variable=self.pymol_var).grid(column=0, row=3, padx=5,
                                                                                            pady=1)

        self.vmd_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(system_frame, text="Set Balls VMD", variable=self.vmd_var).grid(column=1, row=3, padx=5, pady=1)

        self.info_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(system_frame, text="Information", variable=self.info_var).grid(column=2, row=3, padx=5, pady=1)

        # Group frame
        group_frame = ttk.LabelFrame(main_frame, text="Group")
        group_frame.pack(fill="x", padx=5, pady=5)

        self.group_balls_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(group_frame, text="Group Balls", variable=self.group_balls_var).pack(anchor="w", padx=5, pady=2)

        self.surrounding_balls_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(group_frame, text="Surrounding Balls", variable=self.surrounding_balls_var).pack(anchor="w",
                                                                                                         padx=5, pady=2)

        # Buttons frame
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill="x", pady=10)

        ttk.Button(button_frame, text="Apply", command=self._on_ok).pack(side="right", padx=5)
        ttk.Button(button_frame, text="Cancel", command=self._on_cancel).pack(side="right", padx=5)

        # Center the window on the parent
        self.update_idletasks()
        width = self.winfo_width()
        height = self.winfo_height()
        x = (self.winfo_screenwidth() // 2) - (width // 2)
        y = (self.winfo_screenheight() // 2) - (height // 2)
        self.geometry(f"{width}x{height}+{x}+{y}")

    def _on_ok(self):
        """Handle OK button click."""
        # Here you would collect all the checkbox values and update the export settings
        self.destroy()

    def _on_cancel(self):
        """Handle Cancel button click."""
        self.destroy()
