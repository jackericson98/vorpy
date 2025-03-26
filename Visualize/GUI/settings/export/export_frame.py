import os
import sys
from pathlib import Path

# Add the project root directory to the Python path
project_root = Path(__file__).resolve().parents[4]
sys.path.append(str(project_root))

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog


class ExportFrame(ttk.LabelFrame):
    """
    A frame for export settings configuration.
    """
    def __init__(self, parent, gui):
        super().__init__(parent, text="Export Settings", padding="5")
        self.gui = gui
        
        # Create and pack widgets
        self._create_widgets()
        
        # Initialize custom mode flag
        self.is_custom = False
        
    def _create_widgets(self):
        """Create and pack all widgets in the frame."""
        # Export Size Section
        ttk.Label(self, text="Export Size:").grid(row=0, column=0, sticky="w", padx=2, pady=(2, 0))
        
        # Radio buttons for export size
        self.export_size = tk.StringVar(value="Medium")
        ttk.Radiobutton(self, text="Small", variable=self.export_size, value="Small").grid(row=1, column=0, sticky="w", padx=2, pady=1)
        ttk.Radiobutton(self, text="Medium", variable=self.export_size, value="Medium").grid(row=1, column=1, sticky="w", padx=2, pady=1)
        ttk.Radiobutton(self, text="Large", variable=self.export_size, value="Large").grid(row=1, column=2, sticky="w", padx=2, pady=1)
        ttk.Radiobutton(self, text="All", variable=self.export_size, value="All").grid(row=1, column=3, sticky="w", padx=2, pady=1)
        
        # Custom Button
        self.custom_button = ttk.Button(self, text="Custom", command=self.toggle_custom, width=8)
        self.custom_button.grid(row=1, column=4, sticky="w", padx=2, pady=1)
        
        # Export Location
        ttk.Label(self, text="Location:").grid(row=2, column=0, sticky="e", padx=2, pady=2)
        self.export_location = ttk.Entry(self, width=25)
        self.export_location.grid(row=2, column=1, columnspan=3, sticky="w", padx=2, pady=2)
        self.export_location.insert(0, "Default Output Directory")
        
        # Browse Button
        browse_button = ttk.Button(self, text="Browse", command=self.choose_export_location, width=8)
        browse_button.grid(row=2, column=4, sticky="w", padx=2, pady=2)
    
    def toggle_custom(self):
        """Toggle custom mode and handle radio button selection."""
        self.is_custom = not self.is_custom
        if self.is_custom:
            self.export_size.set("")  # Deselect all radio buttons
            self.custom_button.state(['pressed'])  # Visual feedback
            self.open_custom_settings()
        else:
            self.export_size.set("Medium")  # Reset to default
            self.custom_button.state(['!pressed'])
        
    def choose_export_location(self):
        """Open a directory chooser dialog for export location."""
        directory = filedialog.askdirectory(title='Choose Export Location')
        if directory:
            self.export_location.delete(0, tk.END)
            self.export_location.insert(0, directory)
    
    def open_custom_settings(self):
        """Open the custom export settings window."""
        CustomExportWindow(self.gui)


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
        ttk.Label(system_frame, text="Ball Files:", font=('Helvetica', 10, 'bold')).grid(column=0, columnspan=5, row=0, padx=5, pady=5, sticky="e")
        
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
        ttk.Label(system_frame, text="Other Files:", font=('Helvetica', 10, 'bold')).grid(column=0, columnspan=5, row=2, padx=5, pady=5, sticky="e")
        
        self.pymol_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(system_frame, text="Set Balls PyMOL", variable=self.pymol_var).grid(column=0, row=3, padx=5, pady=1)
        
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
        ttk.Checkbutton(group_frame, text="Surrounding Balls", variable=self.surrounding_balls_var).pack(anchor="w", padx=5, pady=2)
        
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


if __name__ == "__main__":
    root = tk.Tk()
    export_frame = ExportFrame(root, None)
    export_frame.pack(fill="both", expand=True)
    root.mainloop()
