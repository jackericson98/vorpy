import tkinter as tk
from tkinter import ttk


class SystemExportsWindow:
    """
    Opens a new window to display system exports.
    """
    def __init__(self, gui):
        """
        Initializes the SystemExportsWindow.

        Args:
            gui: The main GUI application object.
        """
        self.gui = gui
        # Use gui directly as parent since it's a tk.Tk instance
        self.window = tk.Toplevel(self.gui)
        self.window.title("System Exports")
        self.window.geometry("250x270") 

        # Create main frame
        self.main_frame = ttk.Frame(self.window)
        self.main_frame.pack(fill="both", expand=True, padx=10, pady=5)

        # Create system files frame
        self.system_files_frame = ttk.LabelFrame(self.main_frame, text="System Files")
        self.system_files_frame.pack(fill="x", padx=5, pady=5)

        # Create a check box for the set atomic radii file
        self.set_radii_var = tk.BooleanVar(value=True)
        self.set_radii_check = ttk.Checkbutton(self.system_files_frame, text="Set Atomic Radii",
                                               variable=self.set_radii_var)
        self.set_radii_check.pack(padx=5, pady=5)

        # Create a check box for the info file
        self.info_var = tk.BooleanVar(value=True)
        self.info_check = ttk.Checkbutton(self.system_files_frame, text="Info File", variable=self.info_var)
        self.info_check.pack(padx=5, pady=5)

        # Create ball files frame
        self.ball_files_frame = ttk.LabelFrame(self.main_frame, text="Ball Files")
        self.ball_files_frame.pack(fill="x", padx=5, pady=5)

        # Create left and right columns for ball files
        self.left_column = ttk.Frame(self.ball_files_frame)
        self.left_column.pack(side="left", fill="both", expand=True, padx=5)
        
        self.right_column = ttk.Frame(self.ball_files_frame)
        self.right_column.pack(side="left", fill="both", expand=True, padx=5)

        # Create checkboxes for the different ball file types (pdb, mol, gro, xyz, cif, txt)
        self.pdb_var = tk.BooleanVar(value=True)
        self.pdb_check = ttk.Checkbutton(self.left_column, text="PDB", variable=self.pdb_var)
        self.pdb_check.pack(padx=5, pady=5)

        self.mol_var = tk.BooleanVar(value=False)
        self.mol_check = ttk.Checkbutton(self.left_column, text="MOL", variable=self.mol_var) 
        self.mol_check.pack(padx=5, pady=5)

        self.gro_var = tk.BooleanVar(value=False)
        self.gro_check = ttk.Checkbutton(self.left_column, text="GRO", variable=self.gro_var)
        self.gro_check.pack(padx=5, pady=5)    

        self.xyz_var = tk.BooleanVar(value=False)
        self.xyz_check = ttk.Checkbutton(self.right_column, text="XYZ", variable=self.xyz_var)
        self.xyz_check.pack(padx=5, pady=5)

        self.cif_var = tk.BooleanVar(value=False)
        self.cif_check = ttk.Checkbutton(self.right_column, text="CIF", variable=self.cif_var)
        self.cif_check.pack(padx=5, pady=5)

        self.txt_var = tk.BooleanVar(value=False)   
        self.txt_check = ttk.Checkbutton(self.right_column, text="TXT", variable=self.txt_var)
        self.txt_check.pack(padx=5, pady=5)

        # Create button frame
        self.button_frame = ttk.Frame(self.main_frame)
        self.button_frame.pack(fill="x", pady=10)

        # Create an apply button and a cancel button
        self.apply_button = ttk.Button(self.button_frame, text="Apply", command=self.apply_exports)
        self.apply_button.pack(side="right", padx=5)

        self.cancel_button = ttk.Button(self.button_frame, text="Cancel", command=self.cancel_exports)
        self.cancel_button.pack(side="right", padx=5)
        
    def apply_exports(self):
        """Apply the system exports. Update the gui.sys.exports dictionary."""
        self.gui.exports = {
            'set_radii': self.set_radii_var.get(),
            'info': self.info_var.get(),
            'pdb': self.pdb_var.get(),
            'mol': self.mol_var.get(),
            'gro': self.gro_var.get(),
            'xyz': self.xyz_var.get(),
            'cif': self.cif_var.get(),
            'txt': self.txt_var.get()
        }
        # Close the window
        self.window.destroy()

    def cancel_exports(self):
        """Cancel the system exports."""
        self.window.destroy()


if __name__ == "__main__":
    root = tk.Tk()
    SystemExportsWindow(root)
    root.mainloop()


