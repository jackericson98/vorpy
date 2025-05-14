import tkinter as tk
from tkinter import ttk, messagebox
from importlib import resources


class TestFileBrowserApp:
    def __init__(self, gui, sys):
        self.root = tk.Toplevel(gui)
        self.root.title("Select a Packaged Test File")
        self.sys = sys
        
        # Make window modal
        self.root.transient(gui)
        self.root.grab_set()
        
        # Set window size and position
        window_width = 275
        window_height = 120
        screen_width = self.root.winfo_screenwidth()
        screen_height = self.root.winfo_screenheight()
        x = (screen_width - window_width) // 2
        y = (screen_height - window_height) // 2
        self.root.geometry(f"{window_width}x{window_height}+{x}+{y}")

        # Fetch available data files
        self.files = self.get_data_files()
        if not self.files:
            messagebox.showerror("Error", "No test files found in vorpy.data")
            self.root.destroy()
            return

        # Create main frame
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.pack(fill="both", expand=True)

        # Selected file variable
        self.selected_file = tk.StringVar(value=self.files[0])

        # Dropdown label
        ttk.Label(main_frame, text="Select a test file:").pack(pady=(10, 0))

        # Dropdown menu
        self.dropdown = ttk.OptionMenu(main_frame, self.selected_file, self.files[0], *self.files)
        self.dropdown.pack(pady=5)

        # Preview button
        # Create a frame for the buttons
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(pady=5)

        # Add buttons side by side
        ttk.Button(button_frame, text="Preview", command=self.open_selected_file).pack(side="left", padx=5)
        ttk.Button(button_frame, text="Cancel", command=self.cancel).pack(side="left", padx=5)
        ttk.Button(button_frame, text="Save", command=self.save_selection).pack(side="left", padx=5)

    def get_data_files(self):
        # Return all file names in vorpy.data folder
        return [
            entry.name
            for entry in resources.files("vorpy.data").iterdir()
            if entry.is_file()
        ]

    def open_selected_file(self):
        file_name = self.selected_file.get()
        file_path = resources.files("vorpy.data").joinpath(file_name)

        try:
            with file_path.open("r") as f:
                content = f.read()
            self.show_preview(file_name, content)
        except Exception as e:
            messagebox.showerror("Error", f"Could not open file:\n{e}")

    def show_preview(self, filename, content):
        preview_window = tk.Toplevel(self.root)
        preview_window.title(f"Preview: {filename}")

        text_area = tk.Text(preview_window, wrap="word", height=20, width=80)
        text_area.insert("1.0", content[:20000])  # limit preview to 2000 chars
        text_area.config(state="disabled")
        text_area.pack(padx=10, pady=10)

    def save_selection(self):
        file_name = self.selected_file.get()
        file_path = resources.files("vorpy.data").joinpath(file_name)
        self.sys.files['ball_file'] = str(file_path)
        self.root.destroy()

    def cancel(self):
        self.sys.files['ball_file'] = None
        self.root.destroy()

if __name__ == "__main__":
    root = tk.Tk()
    app = TestFileBrowserApp(root)
    root.mainloop()
