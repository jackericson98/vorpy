import tkinter as tk
from tkinter import ttk, filedialog


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
        # System info Frame
        sys_info_frame = ttk.LabelFrame(parent, text=" System info ")
        sys_info_frame.pack(fill="both", padx=10, pady=5)

        # Configure grid weights for centering
        sys_info_frame.grid_columnconfigure(0, weight=1)
        sys_info_frame.grid_columnconfigure(1, weight=2)
        sys_info_frame.grid_columnconfigure(2, weight=1)

        # System Name in the top center
        system_name = "System Name" if gui is None else gui.sys.name
        font = ('Helvetica', 12) if gui is None else gui.fonts['class 1']
        tk.Label(sys_info_frame, text=system_name, font=font).grid(row=0, column=0, columnspan=3, pady=10)

        # Input File Section
        tk.Label(sys_info_frame, text="Input File:", font=('Helvetica', 10) if gui is None else gui.fonts['class 2']).grid(row=1, column=0, sticky="w", padx=5, pady=2)
        input_file_label = tk.Label(sys_info_frame, text="file", font=('Helvetica', 10) if gui is None else gui.fonts['class 2'])
        input_file_label.grid(row=1, column=1)
        ttk.Button(sys_info_frame, text="Browse", command=lambda: None if gui is None else gui.choose_ball_file).grid(row=1, column=2, sticky="e", padx=5, pady=2)

        # Other Files Section
        tk.Label(sys_info_frame, text="Other Files:", font=('Helvetica', 10) if gui is None else gui.fonts['class 2']).grid(row=2, column=0, sticky="w", padx=5, pady=2)
        self.other_files_label = tk.Label(sys_info_frame, text="None", font=('Helvetica', 10) if gui is None else gui.fonts['class 2'])
        self.other_files_label.grid(row=2, column=1)
        ttk.Button(sys_info_frame, text="Add", command=self._browse_other_files).grid(row=2, column=2, sticky="e", padx=5, pady=2)

        # Output Directory Section
        tk.Label(sys_info_frame, text="Output Directory:", font=('Helvetica', 10) if gui is None else gui.fonts['class 2']).grid(row=3, column=0, sticky="w", padx=5, pady=2)
        output_dir_label = tk.Label(sys_info_frame, text="None", font=('Helvetica', 10) if gui is None else gui.fonts['class 2'])
        output_dir_label.grid(row=3, column=1)
        ttk.Button(sys_info_frame, text="Browse", command=lambda: None if gui is None else gui.choose_output_directory).grid(row=3, column=2, sticky="e", padx=5, pady=2)

    def _browse_other_files(self):
        """Open file dialog to select other files."""
        filename = filedialog.askfilename(
            title="Select Other File",
            filetypes=[("All files", "*.*")]
        )
        if filename:
            self.other_files_label.config(text=filename)

    # Subframe for Specific info

if __name__ == "__main__":
    root = tk.Tk()
    root.title("System Information")
    sys_info_frame = SystemFrame(None, root)
    root.mainloop()