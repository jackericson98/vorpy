import tkinter as tk
from tkinter import ttk


class HelpWindow(tk.Toplevel):
    """
    A window containing help information for the GUI.
    """
    def __init__(self, parent):
        super().__init__(parent)
        
        # Configure window
        self.title("VorPy Help")
        self.geometry("800x600")  # Increased window size
        self.resizable(False, False)
        
        # Make window modal
        self.transient(parent)
        self.grab_set()
        
        # Create main frame
        main_frame = ttk.Frame(self, padding="10")
        main_frame.pack(fill="both", expand=True)
        
        # Add header
        header = ttk.Label(main_frame, text="VorPy Help", font=("Arial", 16, "bold"))
        header.pack(pady=(0, 10))
        
        # Create notebook for tabs
        notebook = ttk.Notebook(main_frame)
        notebook.pack(fill="both", expand=True)
        
        # Create tabs
        self._create_about_tab(notebook)
        self._create_system_info_tab(notebook)
        self._create_groups_tab(notebook)
        self._create_build_settings_tab(notebook)
        self._create_export_settings_tab(notebook)
        
        # Close button
        close_button = ttk.Button(main_frame, text="Close", command=self.destroy)
        close_button.pack(pady=10)
        
        # Center the window on the parent
        self.update_idletasks()
        width = self.winfo_width()
        height = self.winfo_height()
        x = (self.winfo_screenwidth() // 2) - (width // 2)
        y = (self.winfo_screenheight() // 2) - (height // 2)
        self.geometry(f"{width}x{height}+{x}+{y}")

    def _create_about_tab(self, notebook):
        """Create the About tab."""
        frame = ttk.Frame(notebook, padding="10")
        notebook.add(frame, text="About")
        
        # Configure grid weights
        frame.grid_rowconfigure(1, weight=1)
        frame.grid_columnconfigure(0, weight=1)
        
        # Title
        title = ttk.Label(frame, text="About VorPy", font=("Arial", 14, "bold"))
        title.grid(row=0, column=0, pady=(0, 10), sticky="w")
        
        # Content
        content = """
        VorPy is a comprehensive Voronoi diagram calculation tool designed for molecular analysis and network 
        generation.
        
        Features:
        • Multiple network types support:
          - Additively Weighted Voronoi
          - Power Diagram
          - Primitive (Delaunay) triangulation
        • Flexible group management for different molecular components
        • Customizable surface settings and parameters
        • Comprehensive export options for analysis results
        
        Purpose:
        This tool is designed to help researchers and scientists analyze molecular structures and generate various types 
        of Voronoi networks for their analysis. It provides a user-friendly interface for managing complex calculations 
        and visualizing results.
        
        Usage:
        1. Configure your system information
        2. Set up groups for different molecular components
        3. Adjust build settings for network generation
        4. Configure export settings for results
        5. Run the analysis
        
        The program will process your input and generate the requested networks, saving results according to your 
        specified settings.
        """
        text = tk.Text(frame, wrap="word", width=70)
        text.insert("1.0", content)
        text.config(state="disabled")
        text.grid(row=1, column=0, sticky="nsew", padx=5, pady=5)

    def _create_system_info_tab(self, notebook):
        """Create the System Information tab."""
        frame = ttk.Frame(notebook, padding="10")
        notebook.add(frame, text="System Information")
        
        # Configure grid weights
        frame.grid_rowconfigure(1, weight=1)
        frame.grid_columnconfigure(0, weight=1)
        
        # Title
        title = ttk.Label(frame, text="System Information Frame", font=("Arial", 12, "bold"))
        title.grid(row=0, column=0, pady=(0, 10), sticky="w")
        
        # Content
        content = """
        The System Information frame displays and manages the basic system configuration:
        
        • System Name: Shows the name of the current system
        • Input File: Select the main input file for the system
        • Other Files: Shows any additional files associated with the system
        • Output Directory: Choose where to save the results
        
        Use the "Select File" and "Select Directory" buttons to choose your input and output locations.
        """
        text = tk.Text(frame, wrap="word", width=70)
        text.insert("1.0", content)
        text.config(state="disabled")
        text.grid(row=1, column=0, sticky="nsew", padx=5, pady=5)

    def _create_groups_tab(self, notebook):
        """Create the Groups tab."""
        frame = ttk.Frame(notebook, padding="10")
        notebook.add(frame, text="Groups")
        
        # Configure grid weights
        frame.grid_rowconfigure(1, weight=1)
        frame.grid_columnconfigure(0, weight=1)
        
        # Title
        title = ttk.Label(frame, text="Groups Frame", font=("Arial", 12, "bold"))
        title.grid(row=0, column=0, pady=(0, 10), sticky="w")
        
        # Content
        content = """
        The Groups frame manages different groups of molecules in your system:
        
        • Network Type: Choose the type of network to generate
        • Maximum Vertex: Set the maximum number of vertices
        • Box Multiplier: Adjust the size of the bounding box
        • Atomic Radii/Masses: Configure atomic properties
        • Surface Settings: Adjust surface parameters
        
        Each group can have its own settings and will be processed independently.
        """
        text = tk.Text(frame, wrap="word", width=70)
        text.insert("1.0", content)
        text.config(state="disabled")
        text.grid(row=1, column=0, sticky="nsew", padx=5, pady=5)

    def _create_build_settings_tab(self, notebook):
        """Create the Build Settings tab."""
        frame = ttk.Frame(notebook, padding="10")
        notebook.add(frame, text="Build Settings")
        
        # Configure grid weights
        frame.grid_rowconfigure(1, weight=1)
        frame.grid_columnconfigure(0, weight=1)
        
        # Title
        title = ttk.Label(frame, text="Build Settings Frame", font=("Arial", 12, "bold"))
        title.grid(row=0, column=0, pady=(0, 10), sticky="w")
        
        # Content
        content = """
        The Build Settings frame controls how the networks are generated:
        
        • Max Vert Rad: Maximum vertex radius (0.01-500)
        • Max Box Multi: Maximum box multiplier (1-200)
        • Network Type: Choose from various network types:
          - Additively Weighted
          - Power Diagram
          - Primitive (Delaunay)
          - Compare options for different network types
        • Surface Settings: Configure surface parameters
        
        These settings affect how the networks are constructed and analyzed.
        """
        text = tk.Text(frame, wrap="word", width=70)
        text.insert("1.0", content)
        text.config(state="disabled")
        text.grid(row=1, column=0, sticky="nsew", padx=5, pady=5)

    def _create_export_settings_tab(self, notebook):
        """Create the Export Settings tab."""
        frame = ttk.Frame(notebook, padding="10")
        notebook.add(frame, text="Export Settings")
        
        # Configure grid weights
        frame.grid_rowconfigure(1, weight=1)
        frame.grid_columnconfigure(0, weight=1)
        
        # Title
        title = ttk.Label(frame, text="Export Settings Frame", font=("Arial", 12, "bold"))
        title.grid(row=0, column=0, pady=(0, 10), sticky="w")
        
        # Content
        content = """
        The Export Settings frame controls how results are saved:
        
        • Output Format: Choose the format for saving results
        • Data Selection: Select which data to export
        • File Options: Configure file naming and organization
        • Export Location: Choose where to save the results
        
        These settings determine how your analysis results are saved and organized.
        """
        text = tk.Text(frame, wrap="word", width=70)
        text.insert("1.0", content)
        text.config(state="disabled")
        text.grid(row=1, column=0, sticky="nsew", padx=5, pady=5)


if __name__ == "__main__":
    root = tk.Tk()
    help_window = HelpWindow(root)
    root.mainloop() 
