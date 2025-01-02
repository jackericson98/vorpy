import tkinter as tk
from tkinter import filedialog
from System.system import System
from Visualize.GUI.Information.information_frame import create_information_section
from Visualize.GUI.settings.settings_frame import create_settings_section
import os

"""
This GUI operates the whole VorPy interface. Once running, the command line that it was run out of will inform you on 
your progress. The GUI is only for launching the program.

GUI Options:



GUI returns
"""


class VorPyGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        # Create a default system
        self.sys = System(simple=True)

        # Create a variable for the settings dictionary
        self.settings = None
        # Create a variable for the exports dictionary
        self.exports = None

        self.title("VorPy")

        # Font classes
        self.fonts = {'class 1': ("Arial", 16), 'class 2': ("Arial", 10), 'class 3': ("Arial", 12, "bold"),
                      'class 4': ("Arial", 14)}

        # Title Section
        title_label = tk.Label(self, text="VorPy", font=self.fonts['class 1'], pady=10)
        title_label.pack()

        # Main Frame to replace notebook
        main_frame = tk.Frame(self)
        main_frame.pack(expand=True, fill="both")

        # Information Section
        info_frame = tk.Frame(main_frame, height=300, width=500)
        info_frame.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")
        self.create_information_section(info_frame)

        # Settings Section
        settings_frame = tk.Frame(main_frame, height=300, width=500)
        settings_frame.grid(row=0, column=1, padx=10, pady=10, sticky="nsew")
        self.create_settings_section(settings_frame)

        # Configure equal sizing of columns and rows
        main_frame.grid_columnconfigure(0, weight=1)
        main_frame.grid_columnconfigure(1, weight=1)
        main_frame.grid_rowconfigure(0, weight=1)

        # Run and Cancel Buttons
        button_frame = tk.Frame(self, pady=10)
        button_frame.pack()

        run_button = tk.Button(button_frame, text="Run", command=self.run_program, font=self.fonts['class 2'])
        run_button.pack(side="left", padx=5)

        cancel_button = tk.Button(button_frame, text="Cancel", command=self.quit, font=self.fonts['class 2'])
        cancel_button.pack(side="left", padx=5)

    def create_information_section(self, frame):
        create_information_section(self, frame)

    def create_settings_section(self, frame):
        create_settings_section(self, frame)

    def open_surface_settings_gui(self):
        print("Opening Surface Settings GUI...")

    def run_program(self):
        print("Running the program...")

    def choose_ball_file(self):
        self.sys.files['base_file'] = filedialog.askopenfilename(title='Choose Base File')
        print(f"Ball file selected: {self.sys.files['base_file']}")

    def choose_output_directory(self):
        self.sys.files['dir'] = filedialog.askdirectory(title='Choose Output Directory')
        print(f"Output directory selected: {self.sys.files['dir']}")

    def add_group(self):
        print("Adding a new group...")


if __name__ == "__main__":
    os.chdir('../..')
    app = VorPyGUI()
    app.mainloop()

