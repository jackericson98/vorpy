import tkinter as tk
from tkinter import ttk


"""
This file updates the group info. If you add more than one group, the gui will update it.

Options include:
    Defaults:
        Centered Text saying no molecule is loaded, but spaced out for options
    Defaults Once a file is loaded:
        1. Molecule: No Sol - this means as soon as a molecule file is loaded anything marked SOL,
           HOH, WAT, etc... is excluded
        2. Foam:


"""


class GroupsFrame:
    """
    Builds the groups configuration frame with the specified layout.
    Handles group management, network type selection, and various group settings.

    Args:
        gui: The main GUI application object.
        parent: The parent frame to which this groups frame will be added.
    """
    def __init__(self, gui, parent):
        """
        Initialize the groups frame with all necessary components.
        """
        self.gui = gui
        
        # Main Groups Frame
        self.group_settings_frame = ttk.LabelFrame(parent, text=" Groups ")
        self.group_settings_frame.pack(fill="both", padx=10, pady=5)

        # Add tabs for groups
        self.group_notebook = ttk.Notebook(self.group_settings_frame)
        self.group_notebook.pack(fill="both", expand=True)

        # Initialize the frame based on system state
        self._initialize_groups()

    def _initialize_groups(self):
        """
        Initialize the groups frame based on the system state.
        Creates appropriate tabs and settings for each group.
        """
        # Default: If no file is loaded
        if not hasattr(self.gui.sys, 'groups') or self.gui.sys.groups is None:
            self._create_no_groups_tab("No file is loaded or no groups available.")
            return

        # Default: There are no groups
        if len(self.gui.sys.groups) == 0:
            self._create_no_groups_tab("System Must Be Loaded")
            return

        # Create tabs for each group
        for group in self.gui.sys.groups:
            self._create_group_tab(group)

    def _create_no_groups_tab(self, message):
        """
        Create a tab for when no groups are available.
        
        Args:
            message: The message to display in the tab.
        """
        no_groups_frame = ttk.Frame(self.group_notebook)
        self.group_notebook.add(no_groups_frame, text="No Groups")
        tk.Label(no_groups_frame, text=message, font=self.gui.fonts['class 2']).pack(
            padx=20, pady=20)

    def _create_group_tab(self, group):
        """
        Create a tab for a specific group with all its settings.
        
        Args:
            group: The group object to create settings for.
        """
        group_frame = ttk.Frame(self.group_notebook)
        self.group_notebook.add(group_frame, text=group.name)

