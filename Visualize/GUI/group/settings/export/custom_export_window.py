import tkinter as tk
from tkinter import ttk


class CustomExportWindow(tk.Toplevel):
    """Window for custom export settings."""

    def __init__(self, parent):
        super().__init__(parent)
        self.title("Custom Export Settings")
        self.geometry("800x600")
        self.resizable(False, False)
        if not isinstance(parent, tk.Tk):
            self.transient(parent)
            self.grab_set()

        # Main frame
        main_frame = ttk.Frame(self)
        main_frame.pack(fill="both", expand=True, padx=10, pady=10)

        # Info Frame
        info_frame = ttk.LabelFrame(main_frame, text="Info")
        info_frame.pack(fill="x", padx=5, pady=5)

        # Info checkbuttons
        self.logs_var = tk.BooleanVar(value=True)
        self.verts_var = tk.BooleanVar(value=False)
        self.info_var = tk.BooleanVar(value=False)

        ttk.Checkbutton(info_frame, text="Logs", variable=self.logs_var).pack(anchor="w", padx=5, pady=2)
        ttk.Checkbutton(info_frame, text="Verts", variable=self.verts_var).pack(anchor="w", padx=5, pady=2)
        ttk.Checkbutton(info_frame, text="Info", variable=self.info_var).pack(anchor="w", padx=5, pady=2)

        # Balls Frame
        balls_frame = ttk.LabelFrame(main_frame, text="Balls")
        balls_frame.pack(fill="x", padx=5, pady=5)

        # Create a frame for the column headers
        header_frame = ttk.Frame(balls_frame)
        header_frame.pack(fill="x", padx=5, pady=5)

        # Column headers
        ttk.Label(header_frame, text="Group\nBalls").grid(row=0, column=0, padx=20)
        ttk.Label(header_frame, text="Surrounding\nBalls").grid(row=0, column=1, padx=20)

        # Create a frame for the checkbuttons
        checkbuttons_frame = ttk.Frame(balls_frame)
        checkbuttons_frame.pack(fill="x", padx=5, pady=5)

        # File format options
        formats = ['.pdb', '.cif', '.mol', '.gro', '.xyz', '.txt']

        # Create variables for checkbuttons
        self.group_vars = {}
        self.surrounding_vars = {}

        # Create checkbuttons for each format
        for i, fmt in enumerate(formats):
            # Group balls checkbuttons
            self.group_vars[fmt] = tk.BooleanVar(value=False)
            ttk.Checkbutton(checkbuttons_frame, text=fmt, 
                          variable=self.group_vars[fmt]).grid(row=i, column=0, padx=20, pady=2)

            # Surrounding balls checkbuttons
            self.surrounding_vars[fmt] = tk.BooleanVar(value=False)
            ttk.Checkbutton(checkbuttons_frame, text=fmt, 
                          variable=self.surrounding_vars[fmt]).grid(row=i, column=1, padx=20, pady=2)

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


if __name__ == '__main__':
    root = tk.Tk()
    app = CustomExportWindow(root)
    root.mainloop()
