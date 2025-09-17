import os
import sys
import matplotlib as mpl
from matplotlib._api.deprecation import MatplotlibDeprecationWarning as MPLDepWarn

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.system import System

try:
    from vorpy.src.GUI.vorpy_gui import VorPyGUI
    HAS_GUI = True
except ImportError:
    HAS_GUI = False


class Run:
    def __init__(self, file=None, load_files=None, settings=None, groups=None, exports=None):
        """
        Initialize the Run class for running the VorPy program. If no file is loaded, the GUI will be launched.. 

        Args:
            file (str): The file to run the VorPy program on.
            load_files (list): The files to load into the program.
            settings (dict): The settings for the VorPy program.
            groups (list): The groups to build.
            exports (list): The exports to run.

        Returns:
            None
        """

        self.file = file
        self.load_files = load_files
        self.settings = settings
        self.groups = groups
        self.exports = exports

        # Check if no file is loaded, if so launch the GUI
        if self.file is None:
            if HAS_GUI:
                gui = VorPyGUI()
                gui.mainloop()
                return
            else:
                print("GUI not available (tkinter not installed). Please provide a file to process.")
                return
        
        # Get the settings
        self.get_settings()

        # Get the groups
        self.get_groups()

        # Create the system
        self.system = System(self.file)
    
    def get_settings(self):
        """
        Get the settings for the VorPy program.
        """
        # Initialize the settings dictionary
        settings = {'net_type': None, 'surf_res': None, 'box_size': None, 'max_vert': None, 'build_type': None, 'net_file': None, 'net_dir': None, 'net_name': None, 'net_format': None, 'net_format_file': None, 'net_format_dir': None, 'net_format_name': None}

        # Loop through the settings to be able to set the settings
        for key, value in self.settings.items():
            # Check for the net type in the settings
            if key.lower() in {'net type', 'nt', 'network type', 'net_type', 'network_type'}:
                # Check that the value is a valid net type
                nt_dict = {
                    **{_: 'aw' for _ in {'voronoi', 'vor', 'additively weighted', 'aw', 'additively_weighted', 'apolonius'}},
                    **{_: 'prm' for _ in {'primitive', 'prm', 'delaunay', 'del'}},
                    **{_: 'pow' for _ in {'power', 'pow', 'power diagram', 'power_diagram'}},
                }
                if value.lower() in nt_dict.keys():
                    settings['net_type'] = nt_dict[value.lower()]
                else:
                    print(f"Error: {value} is not a valid net type.")
                    continue
            # Check for the surface resolution in the settings
            elif key.lower() in {'surf_res', 'surface resolution', 'sr', 'surface_resolution', 'surf_res'}:
                try:
                    settings['surf_res'] = float(value)
                except ValueError:
                    print(f"Error: {value} is not a valid float for the surface resolution.")
            # Check for the box size in the settings
            elif key.lower() in {'box size', 'bs', 'box_size', 'box_size'}:
                try:
                    settings['box_size'] = float(value)
                except ValueError:
                    print(f"Error: {value} is not a valid float for the box size.")
            # Check for the max vert in the settings
            elif key.lower() in {'max vert', 'mv', 'max_vert', 'max_vert', 'probe distance', 'pd', 'probe_distance', 'probe_dist', 'probe dist'}:
                try:
                    settings['max_vert'] = float(value)
                except ValueError:
                    print(f"Error: {value} is not a valid integer for the max vert.")
            # Check for the coloring scheme in the settings
            elif key.lower() in {'coloring scheme', 'cs', 'color_scheme', 'color_scheme', 'col_scheme', 'col_scheme', 'color_scheme', 'color_scheme'}:
                # Check that the value is a valid color scheme
                # Set up the color map
                try:
                    my_cmap = mpl.colormaps.get_cmap(value)
                except MPLDepWarn:
                    my_cmap = mpl.cm.get_cmap(value)
                except AttributeError:
                    my_cmap = mpl.cm.get_cmap(value)
                except ValueError:
                    my_cmap = mpl.cm.get_cmap(value.capitalize())
                except Exception as e:
                    print(f"Error: {value} is not a valid color scheme.")
                    print(f"Error: {e}")
                    continue
                # Set up the color scheme
                settings['surf_col'] = value
            # Check for the coloring scheme factor in the settings
            elif key.lower() in {'coloring scheme factor', 'cs_factor', 'color_scheme_factor', 'color_scheme_factor', 'col_scheme_factor', 'col_scheme_factor', 'color_scheme_factor', 'color_scheme_factor'}:
                csf_dict = {
                    **{_: 'log' for _ in {'log', 'logarithmic', 'log10', 'log2', 'ln', 'natural log'}},
                    **{_: 'sqr' for _ in {'sqr', 'square', 'squared', 'square root'}},
                    **{_: 'cub' for _ in {'cub', 'cubic', 'cubed', 'cubic root'}},
                    **{_: 'lin' for _ in {'lin', 'linear', 'linear scale', 'linear scale', 'linear_scale', 'linear_scale'}},
                    **{_: 'exp' for _ in {'exp', 'exponential', 'exponential scale', 'exponential scale', 'exponential_scale', 'exponential_scale'}}
                }
                if value.lower() in csf_dict.keys():
                    settings['scheme_factor'] = csf_dict[value.lower()]
                else:
                    print(f"Error: {value} is not a valid color scheme factor.")
                    continue
            # Check for the surface Scheme in the settings
            elif key.lower() in {'surface scheme', 'ss', 'surface_scheme', 'surface_scheme', 'surf_scheme', 'surf scheme'}:
                # Check that the value is a valid surface scheme
                ss_dict = {
                    **{_: 'mean' for _ in {'mean', 'mean curv', 'mean_curv', 'mean curvature', 'mean_curvature'}},
                    **{_: 'gauss' for _ in {'gauss', 'gaussian', 'gaussian curv', 'gaussian_curv', 'gaussian curvature', 'gaussian_curvature'}},
                    **{_: 'dist' for _ in {'dist', 'distance'}},
                    **{_: 'in_out' for _ in {'in_out', 'in out', 'in_out', 'in out', 'in_out', 'in out'}}
                }
                if value.lower() in ss_dict.keys():
                    settings['surf_scheme'] = ss_dict[value.lower()]
                else:
                    print(f"Error: {value} is not a valid surface scheme.")
                    continue
        # Set the settings
        self.settings = settings

    def get_groups(self):
        """
        Get the groups for the VorPy program.
        """
        # Initialize the groups dictionary
        groups = []

        # Check if the group is a list or a dictionary
        if isinstance(self.groups, dict):
            # Make it a list so it can be iterated through
            self.groups = list(self.groups)
        
        for my_group in self.groups:
            # Initialize the group dictionary
            group = {}
            # Loop through the groups to be able to set the groups
            for key, value in my_group.items():
                try:
                    # Check for commas in the value
                    if ',' in value:
                        value = value.split(',')
                    else:
                        value = [value]
                    # Create the list of numbers to add to the group
                    ndx_list = []
                    # Loop through the value to be able to fix the indexes
                    for ndx_set in value:
                        # Fix the value so that the indexes are a list
                        if '-' in ndx_set:
                            start, end = ndx_set.split('-')
                            ndx_list.extend(list(range(int(start), int(end) + 1)))
                        else:
                            ndx_list.append(int(ndx_set))
                    # Make sure there are no duplicates in the list
                    ndx_list = list(set(ndx_list))

                except ValueError:
                    print(f"Error: {value} is not a valid integer for the group.")
                    continue

                # Check if the key is an atom key
                if key.lower() in {'atoms', 'atom', 'balls', 'ball'}:
                    # Add to the atoms list
                    group['atoms'] = ndx_list
                elif key.lower() in {'molecules', 'molecule', 'mols', 'mol'}:
                    # Add to the molecules list
                    group['molecules'] = ndx_list
                elif key.lower() in {'chains', 'chain', 'chns', 'chn'}:
                    # Add to the chains list
                    group['chains'] = ndx_list
                elif key.lower() in {'residues', 'residue', 'rsds', 'rsd'}:
                    # Add to the residues list
                    group['residues'] = ndx_list
                else:
                    print(f"Error: {key} is not a valid group key.")
                    continue

            # Add the group to the groups list
            groups.append(group)

        # Set the groups
        self.groups = groups


if __name__ == "__main__":
    # Initialize the run
    run = Run()