import tkinter as tk
from tkinter import ttk


class SurfaceSettingsWindow(tk.Toplevel):
    """
    A window for configuring surface settings.
    """
    def __init__(self, parent):
        super().__init__(parent)
        
        # Configure window
        self.title("Surface Settings")
        self.geometry("400x300")
        self.resizable(False, False)
        
        # Make window modal
        self.transient(parent)
        self.grab_set()
        
        # Create main frame
        main_frame = ttk.Frame(self, padding="10")
        main_frame.pack(fill="both", expand=True)
        
        # Surface Settings
        settings_frame = ttk.LabelFrame(main_frame, text="Surface Settings", padding="5")
        settings_frame.pack(fill="both", expand=True, padx=5, pady=5)
        
        # Create settings rows
        self.settings = {}
        self._create_setting_row(settings_frame, "Setting 1", 0)
        self._create_setting_row(settings_frame, "Setting 2", 1)
        self._create_setting_row(settings_frame, "Setting 3", 2)
        self._create_setting_row(settings_frame, "Setting 4", 3)
        
        # Buttons frame
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill="x", pady=10)
        
        # OK and Cancel buttons
        ok_button = ttk.Button(button_frame, text="OK", command=self._on_ok)
        ok_button.pack(side="right", padx=5)
        
        cancel_button = ttk.Button(button_frame, text="Cancel", command=self._on_cancel)
        cancel_button.pack(side="right", padx=5)
        
        # Center the window on the parent
        self.update_idletasks()
        width = self.winfo_width()
        height = self.winfo_height()
        x = (self.winfo_screenwidth() // 2) - (width // 2)
        y = (self.winfo_screenheight() // 2) - (height // 2)
        self.geometry(f"{width}x{height}+{x}+{y}")

    def _create_setting_row(self, parent, label, row):
        """Create a row for a setting with label and entry."""
        ttk.Label(parent, text=label).grid(row=row, column=0, padx=5, pady=5, sticky="w")
        entry = ttk.Entry(parent)
        entry.grid(row=row, column=1, padx=5, pady=5, sticky="ew")
        self.settings[label] = entry

    def _on_ok(self):
        """Handle OK button click."""
        # TODO: Save settings
        self.destroy()

    def _on_cancel(self):
        """Handle Cancel button click."""
        self.destroy() 