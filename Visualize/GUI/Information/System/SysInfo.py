import tkinter as tk
from tkinter import ttk


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
        # System Information Frame
        sys_info_frame = ttk.LabelFrame(parent, text=" System Information ")
        sys_info_frame.pack(fill="both", padx=10, pady=5)

        # System Name in the top center
        tk.Label(sys_info_frame, text=gui.sys.name, font=gui.fonts['class 1']).grid(row=0, column=0, columnspan=3, pady=5,
                                                                                    sticky="n")

        # Input File Section
        tk.Label(sys_info_frame, text="Input File:", font=gui.fonts['class 2']).grid(row=1, column=0, sticky="w", padx=10)
        input_file_label = tk.Label(sys_info_frame, text="file", font=gui.fonts['class 2'])
        input_file_label.grid(row=1, column=1, sticky="w")
        tk.Button(sys_info_frame, text="Select File", command=gui.choose_ball_file).grid(row=1, column=2, sticky="e",
                                                                                         padx=10)

        # Other Files Section
        tk.Label(sys_info_frame, text="Other Files:", font=gui.fonts['class 2']).grid(row=2, column=0, sticky="w", padx=10)
        other_files_label = tk.Label(sys_info_frame, text="None", font=gui.fonts['class 2'])
        other_files_label.grid(row=2, column=1, sticky="w")

        # Output Directory Section
        tk.Label(sys_info_frame, text="Output Directory:", font=gui.fonts['class 2']).grid(row=3, column=0, sticky="w",
                                                                                           padx=10)
        output_dir_label = tk.Label(sys_info_frame, text="None", font=gui.fonts['class 2'])
        output_dir_label.grid(row=3, column=1, sticky="w")
        tk.Button(sys_info_frame, text="Select Directory", command=gui.choose_output_directory).grid(row=3, column=2,
                                                                                                 sticky="e", padx=10)

    # Subframe for Specific Information

