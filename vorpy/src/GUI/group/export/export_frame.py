import sys
from pathlib import Path

# Add the project root directory to the Python path
project_root = Path(__file__).resolve().parents[5]
sys.path.append(str(project_root))

from vorpy.src.GUI.group.export.custom_export_window import CustomExportWindow
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
        self.group_name_entry = group_name_entry
        
        # Initialize settings dictionary
        self.settings = {
            'size': 'Medium',
            'custom_settings': None,
            'directory': 'Default Output Directory',
            'round_to': 3
        }
        
        # Create and pack widgets
        self._create_widgets()
        
    def _create_widgets(self):
        """Create and pack all widgets in the frame."""
        # Export Size Section
        ttk.Label(self, text="Export Size:").grid(row=0, column=0, sticky="w", padx=2, pady=(2, 0))
        
        # Radio buttons for export size
        self.export_size = tk.StringVar(value=self.settings['size'])
        self.export_size.trace_add('write', self._on_size_change)
        
        (ttk.Radiobutton(self, text="Small", variable=self.export_size, value="Small", command=self._on_size_change)
         .grid(row=1, column=0, sticky="w", padx=2, pady=1))
        (ttk.Radiobutton(self, text="Medium", variable=self.export_size, value="Medium", command=self._on_size_change)
         .grid(row=1, column=1, sticky="w", padx=2, pady=1))
        (ttk.Radiobutton(self, text="Large", variable=self.export_size, value="Large", command=self._on_size_change)
         .grid(row=1, column=2, sticky="w", padx=2, pady=1))
        
        # Custom Button
        self.custom_button = ttk.Button(self, text="Custom", command=self._open_custom_settings, width=8)
        self.custom_button.grid(row=1, column=4, sticky="w", padx=2, pady=1)
        
        # Export Location
        ttk.Label(self, text="Location:").grid(row=2, column=0, padx=2, pady=5)
        self.export_location = ttk.Entry(self, width=25)
        self.export_location.grid(row=2, column=1, columnspan=3, sticky="w", padx=2, pady=5)
        self.export_location.insert(0, self.settings['directory'])
        
        # Browse Button
        browse_button = ttk.Button(self, text="Browse", command=self._choose_export_location, width=8)
        browse_button.grid(row=2, column=4, sticky="w", padx=2, pady=5)
        # Round To setting
        ttk.Label(self, text="Round To (Decimal Points):").grid(row=3, column=0, columnspan=2, padx=2, pady=5, sticky="w")

        self.round_to = tk.IntVar(value=self.settings["round_to"])

        round_spin = ttk.Spinbox(
            self,
            from_=0,
            to=12,
            textvariable=self.round_to,
            width=5,
            command=self._on_round_to_change
        )
        round_spin.grid(row=3, column=4, sticky="w", padx=2, pady=5)

        round_spin.bind("<FocusOut>", lambda e: self._on_round_to_change())
        round_spin.bind("<Return>", lambda e: self._on_round_to_change())

    def _on_round_to_change(self):
        try:
            value = int(self.round_to.get())
        except (ValueError, tk.TclError):
            value = 3
            self.round_to.set(value)

        value = max(0, min(12, value))

        self.settings["round_to"] = value

        if self.gui is not None and getattr(self.gui, "sys", None) is not None:
            self.gui.sys.round_to = value

    def _on_size_change(self, *args):
        """Handle changes to the export size selection."""
        size = self.export_size.get()
        if size != "Custom":
            # Update local settings
            self.settings['size'] = size
            self.settings['custom_settings'] = None
            self.custom_button.state(['!pressed'])

    def _choose_export_location(self):
        directory = filedialog.askdirectory(title='Choose Export Location')

        if not directory:
            return

        self.export_location.delete(0, tk.END)
        self.export_location.insert(0, directory)

        self.settings["directory"] = directory

        if self.gui is not None and getattr(self.gui, "sys", None) is not None:
            self.gui.sys.files["dir"] = directory

        print(f"[ExportFrame] directory selected = {directory}")
        print(f"[ExportFrame] entry now = {self.export_location.get()}")
        print(f"[ExportFrame] sys.files['dir'] now = {self.gui.sys.files['dir']}")
    
    def _open_custom_settings(self):
        """Open the custom export settings window."""
        custom_window = CustomExportWindow(self, self.group_name_entry.get())
        self.custom_button.state(['pressed'])
        self.export_size.set("")  # Deselect radio buttons
        
        # Update local settings
        self.settings['size'] = "Custom"
        
        # Wait for the window to close
        self.wait_window(custom_window)
        
        # Update settings if the window was closed with Apply
        if hasattr(custom_window, 'settings'):
            self.settings['size'] = "Custom"
            self.settings['custom_settings'] = custom_window.settings
        else:
            # If cancelled, revert to previous settings
            self.export_size.set(self.settings['size'])
            self.custom_button.state(['!pressed'])

    def get_settings(self):
        """Get the group's export settings."""
        directory = self.export_location.get()

        if directory in {"", "Default Output Directory"}:
            if self.gui is not None and getattr(self.gui, "sys", None) is not None:
                directory = self.gui.sys.files.get("dir", directory)

        self.settings["directory"] = directory

        if hasattr(self, "round_to"):
            self._on_round_to_change()

        print(f"[ExportFrame.get_settings] entry directory = {self.export_location.get()}")
        print(f"[ExportFrame.get_settings] resolved directory = {self.settings['directory']}")

        return self.settings.copy()
    
    def copy_settings_from(self, other_frame):
        """Copy settings from another export frame."""
        self.settings = other_frame.get_settings()

        self.export_size.set(self.settings['size'])
        self.export_location.delete(0, tk.END)
        self.export_location.insert(0, self.settings['directory'])

        self.round_to.set(int(self.settings.get("round_to", 3)))
        self._on_round_to_change()

        if self.settings['size'] == "Custom":
            self.custom_button.state(['pressed'])
        else:
            self.custom_button.state(['!pressed'])


if __name__ == "__main__":
    root = tk.Tk()
    export_frame = ExportFrame(root, None, None)
    export_frame.pack(fill="both", expand=True)
    root.mainloop()
