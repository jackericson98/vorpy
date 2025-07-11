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
        VorPy is a comprehensive Voronoi diagram calculation tool designed for molecular 
        analysis and network generation. 
        

        Features:

        • Multiple network types support:
          - Additively Weighted
          - Power
          - Primitive

        • Flexible group management for different molecular components

        • Customizable surface settings and parameters

        • Comprehensive export options for analysis results

        • Visualization of the network and its components
        

        Purpose:

        This tool is designed to help with the analysis of molecular structures and generate 
        various types of Voronoi networks for their analysis. It provides a user-friendly 
        interface for managing complex calculations and visualizing results.
        

        Usage:

        1. Run the program
        2. Select the input file from the pre-loaded list or from your device.
        3. Select the groupings of atoms (or balls) you want to solve for.
        4. Adjust build settings
        5. Configure export settings for results
        6. Run the analysis, find the outputs in the output directory
        

        The program will process your input and generate the requested networks, saving 
        results according to your specified settings.
        """
        text = tk.Text(frame, wrap="word", width=70)
        text.insert("1.0", content)
        text.config(state="disabled")
        text.grid(row=1, column=0, sticky="nsew", padx=5, pady=5)

    def _create_system_info_tab(self, notebook):
        """Create the System Information tab."""
        frame = ttk.Frame(notebook, padding="10")
        notebook.add(frame, text="System")
        
        # Configure grid weights
        frame.grid_rowconfigure(1, weight=1)
        frame.grid_columnconfigure(0, weight=1)
        
        # Title
        title = ttk.Label(frame, text="System Information", font=("Arial", 12, "bold"))
        title.grid(row=0, column=0, pady=(0, 10), sticky="w")
        
        # Content
        content = """
        The System Information frame provides an overview and management interface for your 
        system's core configuration:

        
        • System Name: Displays the current system's name, which is typically derived from 
          the selected input file.

        • Input File: Allows you to select the primary structure file (e.g., PDB) that 
          defines the molecular system to be analyzed.

        • Additional Files: Lists any supplementary files (such as radii, topology, or 
          parameter files) associated with the system, if applicable.

        • Output Directory: Lets you specify the folder where all generated results and 
          exports will be saved.

          
        To configure your system, use the "Select File" button to choose your main input
        file, and the "Select Directory" button to set the output location. The frame 
        will update to reflect your selections and display relevant file information. This 
        ensures that all subsequent operations and exports are performed using the correct 
        files and directories.
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
        The Groups frame provides a comprehensive interface for managing and configuring 
        distinct groups within your system. Groups allow you to organize atoms, residues, 
        chains, or molecules into logical sets, each of which can be processed independently 
        with its own settings.

        
        Key features of the Groups frame include:

        • Group Selection: Select specific atoms/balls, residues, chains, or molecules from 
          the input file to define the members of each group. This enables targeted analysis 
          or export of particular regions or components of your system.

        • Group Management: 
          - Add or delete groups as needed to organize your workflow.
          - Rename groups to provide meaningful identifiers for each set.
          - Use the provided tabs to easily track, switch between, and manage multiple 
            groups within your project.

        • Group-Specific Build Settings: Each group contains its own build settings (see the 
          "Build Settings" tab for details), allowing you to customize parameters such as 
          network type, maximum vertex count, box size, and surface calculation options on a 
          per-group basis.

        • Group-Specific Export Settings: Configure export options for each group 
          individually (see the "Export Settings" tab for more information). This includes 
          selecting which data to export, output formats, and destination directories.

        • Run Controls: 
          - Run calculations for all groups collectively, or execute only a specific group 
            as needed.
          - Each group can be processed independently, enabling flexible and efficient 
            workflows.

        • Visual Tracking: The Groups frame uses tabs or a list to help you keep track of 
          all defined groups, their settings, and their current status.

          
        The Groups frame is central to organizing your analysis, enabling you to tailor 
        network construction and export options for different parts of your system. By 
        leveraging group-specific settings, you can perform detailed, customized analyses 
        and exports for each region of interest.
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
        The Build Settings frame controls how network is generated for the given group. 
        Choosing the diffrent build settings can create drastically different networks 
        and calculation times. 
        
        
        The main three settings are:

        • Network Type: Select the type of network to construct:
          - Additively Weighted: Curved surfaces, calculates all points closer to one ball 
            than another and is considered the most accurate representation of the spatial 
            decomposition of spheres. 
          - Power Diagram: Flat surfaces, with weights for radii. Represents a happy medium 
            between accurracy 
            and speed. For most molecular applications (volume, surface area, etc.) this is 
            the best choice. 
          - Primitive (Delaunay): Flat surfaces, simple and fast. This diagram only 
            considers the centers of the input balls and therefore loses accuracy in the 
            representation of the spatial decomposition of spheres.

        • Probe Distance: This setting controls how far between balls the brobe reaches 
          before chosing the best vertex for network construction. The larger the probe 
          the more accurate the network, but the slower the generation. Default is 40 
          (Angstroms), with a minimum value of 0.01. 

        • Box Size: Adjusting this value controls the outer reach of the network. The 
          retining box is initially calculates and the multiplier for probe reach. This
          allows the probe to expand to a large size within the input balls without creating 
          surfaces that extend far away from them. The default value for this setting is 
          1.25x, which provides just enough reach to capture the entire system, while 
          limiting any vertices that reach too far out of the system, creating spiked 
          surfaces.
        
          
        For outputs and coloring options, use the Outputs Color Settings button. This will 
        open a seperate window in which the user can set the following settings:
        
        • Surface Resolution: This controls how large the triangles in the triangulation of 
          the output surface are. The smaller they are the more accurate the values. The 
          default is set to 0.2, because in molecular cases, the return on accuracy is not 
          worth the time of construction. For better visualizations a setting of 0.01 to 
          0.05 provides a seemless blend between the triangles in the construction.

        • Color Scheme: This allows the user to choose what specific scheme is being used 
          for the coloring of the surfaces.
          This setting decides what the colors on the surface will show, with the options 
          being:
          - Mean Curvature: Shows the mean curvature of the surface.
          - Gaussian Curvature: Shows the Gaussian curvature of the surface.
          - Distance: Shows the distance from the center of the ball.
          - Overlapping: Shows the overlapping of the surface with other surfaces.
          - No Scheme: Shows the surface in a single color.

        • Surface Colorway: This controls the output colorway of the surfaces. The options 
          for this come from matplotlib's colormaps. Some common examples of colormaps are: 
          Viridis, Plasma, Inferno, Blues, Greens, Reds, etc. A complete list can be found 
          at: https://matplotlib.org/stable/users/explain/colors/colormaps.html

        • Surface Coloring Factor: This controls how the values are represented on the 
          surface coloring. Because the linear representation of the values is not always 
          the best way to visualize them, this setting allows the user to choose 
          a logarithmic, exponential, squared or cubed representation.

        • Concave Colors: This allows the user to choose if the concave regions of the 
          surface should be colored inversely. This is usefule for showing shells of atoms, 
          residues, chains or molecules and shows the full color range from the perspective 
          of the item that the shell is surrounding ranging from the most negative to the most 
          positive curvature values. When selected only shell like objects will be represented 
          this way, since the choice of which perspective would be arbitrary otherwise.

        • Vertex Color: This allows the user to choose the color of the vertices of the 
          network.
          
        • Edge Color: This allows the user to choose the color of the edges of the network.
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
