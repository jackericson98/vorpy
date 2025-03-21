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
        self.is_custom = False
        
        # Export Location
        ttk.Label(self, text="Location:").grid(row=2, column=0, sticky="e", padx=2, pady=2)
        self.export_location = ttk.Entry(self, width=15)
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
        
        # Configure window
        self.title("Custom Export Settings")
        self.geometry("300x250")
        self.resizable(False, False)
        
        # Make window modal
        self.transient(parent)
        self.grab_set()
        
        # Create main frame
        main_frame = ttk.Frame(self, padding="5")
        main_frame.pack(fill="both", expand=True)
        
        # Add settings
        ttk.Label(main_frame, text="Customize Export Settings", font=("Arial", 10, "bold")).pack(pady=(0, 5))
        
        # Add your custom settings widgets here
        
        # Close button
        close_button = ttk.Button(main_frame, text="Close", command=self.destroy, width=8)
        close_button.pack(pady=5)
        
        # Center the window
        self.update_idletasks()
        width = self.winfo_width()
        height = self.winfo_height()
        x = (self.winfo_screenwidth() // 2) - (width // 2)
        y = (self.winfo_screenheight() // 2) - (height // 2)
        self.geometry(f"{width}x{height}+{x}+{y}")


if __name__ == "__main__":
    root = tk.Tk()
    export_frame = ExportFrame(root, None)
    export_frame.pack(fill="both", expand=True)
    root.mainloop()
