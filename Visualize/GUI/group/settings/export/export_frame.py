import sys
from pathlib import Path

# Add the project root directory to the Python path
project_root = Path(__file__).resolve().parents[5]
sys.path.append(str(project_root))

from Visualize.GUI.group.settings.export.custom_export_window import CustomExportWindow
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog


class ExportFrame(ttk.LabelFrame):
    """
    A frame for export settings configuration.
    """
    def __init__(self, parent, gui, group_name_entry):
        super().__init__(parent, text="Export Settings", padding="5")
        self.gui = gui
        # Create and pack widgets
        self._create_widgets()
        self.group_name_entry = group_name_entry
        
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
        CustomExportWindow(self.gui, self.group_name_entry.get())


if __name__ == "__main__":
    root = tk.Tk()
    export_frame = ExportFrame(root, None)
    export_frame.pack(fill="both", expand=True)
    root.mainloop()
