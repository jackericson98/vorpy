import tkinter as tk
from tkinter import ttk, filedialog
from Visualize.GUI.system.radii_adjustments.periodic_table_GUI import PeriodicTableGUI


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
        self.window = tk.Toplevel(self.gui.root)
        self.window.title("System Exports")
        self.window.geometry("300x200") 

        # Create a frame for the system exports
        self.exports_frame = ttk.Frame(self.window)
        self.exports_frame.pack(fill="both", padx=10, pady=5)

        # Create a check box for the set atomic radii file
        self.set_radii_var = tk.BooleanVar(value=False)
        self.set_radii_check = ttk.Checkbutton(self.exports_frame, text="Set Atomic Radii", variable=self.set_radii_var)
        self.set_radii_check.grid(row=0, column=0, padx=5, pady=5)

        # Create a check box for the info file
        self.info_var = tk.BooleanVar(value=False)
        self.info_check = ttk.Checkbutton(self.exports_frame, text="Info File", variable=self.info_var)
        self.info_check.grid(row=0, column=1, padx=5, pady=5)

        # Create checkboxes for the different ball file types (pdb, mol, gro, xyz, cif, txt)
        self.pdb_var = tk.BooleanVar(value=False)
        self.pdb_check = ttk.Checkbutton(self.exports_frame, text="PDB", variable=self.pdb_var)
        self.pdb_check.grid(row=0, column=2, padx=5, pady=5)

        self.mol_var = tk.BooleanVar(value=False)
        self.mol_check = ttk.Checkbutton(self.exports_frame, text="MOL", variable=self.mol_var) 
        self.mol_check.grid(row=0, column=3, padx=5, pady=5)

        self.gro_var = tk.BooleanVar(value=False)
        self.gro_check = ttk.Checkbutton(self.exports_frame, text="GRO", variable=self.gro_var)
        self.gro_check.grid(row=0, column=4, padx=5, pady=5)    

        self.xyz_var = tk.BooleanVar(value=False)
        self.xyz_check = ttk.Checkbutton(self.exports_frame, text="XYZ", variable=self.xyz_var)
        self.xyz_check.grid(row=0, column=5, padx=5, pady=5)

        self.cif_var = tk.BooleanVar(value=False)
        self.cif_check = ttk.Checkbutton(self.exports_frame, text="CIF", variable=self.cif_var)
        self.cif_check.grid(row=0, column=6, padx=5, pady=5)

        self.txt_var = tk.BooleanVar(value=False)   
        self.txt_check = ttk.Checkbutton(self.exports_frame, text="TXT", variable=self.txt_var)
        self.txt_check.grid(row=0, column=7, padx=5, pady=5)

        # Create an apply button and a cancel button
        self.apply_button = ttk.Button(self.exports_frame, text="Apply", command=self.apply_exports)
        self.apply_button.grid(row=0, column=8, padx=5, pady=5)

        self.cancel_button = ttk.Button(self.exports_frame, text="Cancel", command=self.cancel_exports)
        self.cancel_button.grid(row=0, column=9, padx=5, pady=5)
        
    def apply_exports(self):
        """Apply the system exports. Update the gui.sys.exports dictionary."""
        self.gui.sys.exports = {
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


class SystemFrame:
    """
    Builds the system information frame with the specified layout.

    Args:
        gui: The main GUI application object.
        parent: The parent frame to which this system frame will be added.
    """
    def __init__(self, gui, parent):
        """
        The frame that gets the file
        """
        self.gui = gui
        system_frame = ttk.LabelFrame(parent, text=" System ")
        system_frame.pack(fill="both", padx=10, pady=5)

        # System info Frame
        sys_info_frame = ttk.LabelFrame(system_frame, text="Files")
        sys_info_frame.grid(row=1, padx=10, pady=5, sticky="nsew")
        system_frame.grid_rowconfigure(1, weight=0)
        system_frame.grid_columnconfigure(0, weight=1)

        # Configure grid weights for centering
        sys_info_frame.grid_columnconfigure(0, weight=1)
        sys_info_frame.grid_columnconfigure(1, weight=2)
        sys_info_frame.grid_columnconfigure(2, weight=1)

        # System Name in the top center
        self.system_name = tk.StringVar(value="System Name" if gui is None else gui.sys.name)
        font = ('Helvetica', 12) if gui is None else gui.fonts['class 1']
        tk.Label(system_frame, textvariable=self.system_name, font=font).grid(row=0, column=0, columnspan=3, pady=2)

        # Input File Section
        (tk.Label(sys_info_frame, text="Input File:",
                 font=('Helvetica', 10) if gui is None else gui.fonts['class 2'])
         .grid(row=1, column=0, sticky="w", padx=5, pady=2))
        self.input_file_label = tk.Label(sys_info_frame, text="",
                                         font=('Helvetica', 10) if gui is None else gui.fonts['class 2'])
        self.input_file_label.grid(row=1, column=1, sticky='w')
        (ttk.Button(sys_info_frame, text="Browse", command=self.choose_ball_file)
         .grid(row=1, column=2, sticky="e", padx=5, pady=2))

        # Other Files Section
        (tk.Label(sys_info_frame, text="Other Files:",
                 font=('Helvetica', 10) if gui is None else gui.fonts['class 2'])
         .grid(row=2, column=0, sticky="w", padx=5, pady=2))
        
        # Create a frame for the file display and dropdown
        self.files_frame = ttk.Frame(sys_info_frame)
        self.files_frame.grid(row=2, column=1, sticky="w", pady=2)
        
        # Initialize the files list in gui if it doesn't exist
        if gui is not None and 'other_files' not in gui.files:
            gui.files['other_files'] = []
            
        # Create the file display widget
        self.file_display = ttk.Label(self.files_frame, text="",
                                      font=('Helvetica', 10) if gui is None else gui.fonts['class 2'])
        self.file_display.pack(side="left", fill="x", expand=True)
        
        # Create the dropdown (initially hidden)
        self.file_dropdown = ttk.Combobox(self.files_frame, state="readonly", width=50)
        self.file_dropdown.pack(side="left", fill="x", expand=True)
        self.file_dropdown.pack_forget()  # Hide initially
        
        # Update the display based on the number of files
        self._update_file_display()
        
        ttk.Button(sys_info_frame, text="Add", command=self._browse_other_files).grid(row=2, column=2, sticky="e",
                                                                                      padx=5, pady=2)

        # Output Directory Section
        (tk.Label(sys_info_frame, text="Output Directory:",
                  font=('Helvetica', 10) if gui is None else gui.fonts['class 2'])
         .grid(row=3, column=0, sticky="w", padx=5, pady=2))
        self.output_dir_label = tk.Label(sys_info_frame, text="None",
                                         font=('Helvetica', 10) if gui is None else gui.fonts['class 2'])
        self.output_dir_label.grid(row=3, column=1, sticky="w")
        ttk.Button(sys_info_frame, text="Browse", command=self.choose_output_directory).grid(row=3, column=2,
                                                                                             sticky="e", padx=5, pady=2)

        # Create a frame for buttons to change the atomic radii and the system exports
        self.radii_frame = ttk.LabelFrame(system_frame, text="Settings")
        self.radii_frame.grid(row=1, column=1, rowspan=4, sticky="nsew", pady=5, padx=5)
        ttk.Button(self.radii_frame, text="Radii", command=self.change_atomic_radii).grid(row=1, column=0, pady=2,
                                                                                          padx=5, sticky="s")
        ttk.Button(self.radii_frame, text="Exports", command=self.system_exports).grid(row=2, column=0, pady=2, padx=5,
                                                                                       sticky="s")
        ttk.Button(self.radii_frame, text="Reset", command=self.delete_system).grid(row=3, column=0, pady=2, padx=5,
                                                                                    sticky="s")

    def choose_ball_file(self):
        """Open file dialog to select a ball file."""
        filename = filedialog.askopenfilename(
            title="Select Ball File",
            filetypes=[("Ball files", "*.pdb"), ("All files", "*.*")]
        )
        if filename:
            self.gui.ball_file = filename
            self.gui.sys.ball_file = filename
            self.gui.sys.load_sys(filename)
            self.system_name.set(self.gui.sys.name.capitalize())  # Update the display
            
            # Truncate the filename display with ellipses in the middle
            if len(filename) > 50:
                truncated = filename[:20] + "..." + filename[-20:]
            else:
                truncated = filename
            self.input_file_label.config(text=truncated)

    def _browse_other_files(self):
        """Open file dialog to select other files."""
        filename = filedialog.askopenfilename(
            title="Select Other File",
            filetypes=[("All files", "*.*")]
        )
        if filename:
            if self.gui is not None:
                self.gui.files['other_files'].append(filename)
                self._update_file_display()

    def _update_file_display(self, file_string_len=46):
        """Update the display based on the number of files."""
        if self.gui is None or not self.gui.files['other_files']:
            self.file_display.config(text="None")
            self.file_dropdown.pack_forget()
            self.file_display.pack(side="left", fill="x", expand=True)
            return

        files = self.gui.files['other_files']
        if len(files) == 1:
            # Show first 100 characters of the single file
            self.file_display.config(text=files[0][:int(file_string_len / 2) - 2] + "..." +
                                          files[0][-(int(file_string_len / 2) - 2):]
                                          if len(files[0]) > file_string_len else files[0])
            self.file_dropdown.pack_forget()
            self.file_display.pack(side="left", fill="x", expand=True)
        else:
            # Show dropdown with all files
            self.file_display.pack_forget()
            # Create truncated versions of file paths for the dropdown
            truncated_files = [f[:int(file_string_len / 2) - 2] + "..." + f[-(int(file_string_len / 2) - 2):]
                               if len(f) > file_string_len else f for f in files]
            self.file_dropdown['values'] = truncated_files
            self.file_dropdown.set(truncated_files[0])  # Set to first file
            self.file_dropdown.pack(side="left", fill="x", expand=True)
            self.file_dropdown.bind('<<ComboboxSelected>>', self._on_file_selected)

    def _on_file_selected(self, event=None, file_string_len=30):
        """Handle file selection from dropdown."""
        selected = self.file_dropdown.get()
        if selected:
            self.file_display.config(text=selected[:int(file_string_len / 2) - 2] + "..." +
                                          selected[-(int(file_string_len / 2) - 2):]
                                          if len(selected) > file_string_len else selected)

    def choose_output_directory(self):
        """Open directory dialog to select output directory."""
        directory = filedialog.askdirectory(
            title="Select Output Directory"
        )
        if directory:
            if self.gui is not None:
                self.gui.output_dir = directory
                self.gui.sys.output_dir = directory
                
                # Truncate the directory path display with ellipses in the middle
                if len(directory) > 50:
                    truncated = directory[:23] + "..." + directory[-23:]
                else:
                    truncated = directory
                self.output_dir_label.config(text=truncated)

    def change_atomic_radii(self):
        """Open the atomic radii window."""
        PeriodicTableGUI(self.gui)

    def system_exports(self):
        """Open the system exports window."""
        SystemExportsWindow(self.gui)

    def delete_system(self):
        """Delete the system."""
        self.gui.sys.delete()


if __name__ == "__main__":
    root = tk.Tk()
    root.title("System Information")
    sys_info_frame = SystemFrame(None, root)
    root.mainloop()