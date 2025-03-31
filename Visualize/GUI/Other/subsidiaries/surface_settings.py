import tkinter as tk
from tkinter import ttk, messagebox
import warnings
from matplotlib import MatplotlibDeprecationWarning
import matplotlib as mpl

# Suppress MatplotlibDeprecationWarning
warnings.filterwarnings("ignore", category=MatplotlibDeprecationWarning)


# Helper function to validate color maps
def validate_colormap(colormap):
    try:
        mpl.colormaps.get_cmap(colormap)
        return True
    except ValueError:
        messagebox.showerror("Invalid Colormap", f"'{colormap}' is not a valid matplotlib colormap.")
        return False


class SurfaceOptionsGUI:

    def __init__(self):
        self.result = None
        # Create the main application window
        self.root = tk.Tk()

    # Function to handle submission
    def submit(self):
        try:
            resolution = float(self.resolution_var.get())
            if not (0.001 <= resolution <= 5):
                raise ValueError("Surface Resolution must be between 0.001 and 5.")

            colormap = self.colormap_var.get()
            if not validate_colormap(colormap):
                raise ValueError(f"Invalid colormap: {colormap}.")

            shading = self.shading_var.get()
            scale = self.scale_var.get()

            self.result = {
                "Surface Resolution": resolution,
                "Surface Color Spectrum": colormap,
                "Shading Scheme": shading,
                "Coloring Scale": scale
            }

            self.root.destroy()
        except ValueError as e:
            messagebox.showerror("Invalid Input", str(e))

    def run(self):

        self.root.title("Surface Options GUI")

        # Surface Resolution
        resolution_label = ttk.Label(self.root, text="Surface Resolution (0.001 - 5):")
        resolution_label.grid(row=0, column=0, padx=10, pady=10, sticky="w")
        self.resolution_var = tk.StringVar(value="0.2")
        resolution_entry = ttk.Entry(self.root, textvariable=self.resolution_var)
        resolution_entry.grid(row=0, column=1, padx=10, pady=10)

        # Surface Color Spectrum
        colormap_label = ttk.Label(self.root, text="Surface Color Spectrum:")
        colormap_label.grid(row=1, column=0, padx=10, pady=10, sticky="w")
        self.colormap_var = tk.StringVar()
        colormap_combobox = ttk.Combobox(self.root, textvariable=self.colormap_var)
        colormap_combobox["values"] = [
            "viridis",
            "plasma",
            "inferno",
            "magma",
            "cividis",
            "cool",
            "hot",
            "spring",
            "Other (must be a valid matplotlib colormap)"
        ]
        colormap_combobox.grid(row=1, column=1, padx=10, pady=10)
        colormap_combobox.set("viridis")

        verify_button = ttk.Button(self.root, text="Verify", command=validate_colormap(self.colormap_var.get()))
        verify_button.grid(row=1, column=2, padx=10, pady=10)

        # Shading Scheme
        shading_label = ttk.Label(self.root, text="Shading Scheme:")
        shading_label.grid(row=2, column=0, padx=10, pady=10, sticky="w")
        self.shading_var = tk.StringVar()
        shading_combobox = ttk.Combobox(self.root, textvariable=self.shading_var)
        shading_combobox["values"] = [
            "Mean Curvature",
            "Gaussian Curvature",
            "Overlapping",
            "NonOverlapping",
            "Distance from Ball"
        ]
        shading_combobox.grid(row=2, column=1, padx=10, pady=10)
        shading_combobox.set("Mean Curvature")

        # Coloring Scale
        scale_label = ttk.Label(self.root, text="Coloring Scale:")
        scale_label.grid(row=3, column=0, padx=10, pady=10, sticky="w")
        self.scale_var = tk.StringVar()
        scale_combobox = ttk.Combobox(self.root, textvariable=self.scale_var)
        scale_combobox["values"] = ["Log", "Linear", "Squared", "Cubed"]
        scale_combobox.grid(row=3, column=1, padx=10, pady=10)
        scale_combobox.set("Linear")

        # Submit Button
        submit_button = ttk.Button(self.root, text="Submit", command=self.submit)
        submit_button.grid(row=4, column=0, columnspan=2, pady=20)

        # Run the application
        self.root.mainloop()

        return self.result


# Example usage

if __name__ == "__main__":
    gui = SurfaceOptionsGUI()
    options = gui.run()
    if options:
        print("User Selections:", options)

