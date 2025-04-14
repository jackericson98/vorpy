import tkinter as tk
from tkinter import ttk, messagebox
from Visualize.GUI.group.build.build_frame import BuildFrame
from Visualize.GUI.group.export.export_frame import ExportFrame


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


class GroupsFrame(ttk.Frame):
    def __init__(self, parent, gui, settings):
        super().__init__(parent)
        self.gui = gui
        self.settings = settings

        # Create notebook for tabs
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill="both", expand=True, padx=5, pady=5)
        
        # Add initial group tab
        self.add_group_tab("Group 1")
        
        # Add button frame at bottom
        button_frame = ttk.Frame(self)
        button_frame.pack(fill="x", pady=5)
        
        # Create a run button for the current group
        run_button = ttk.Button(button_frame, text="Run", command=lambda: self.run_current_group())
        run_button.pack(side="right", padx=5)

        # Add button for current group        
        add_button = ttk.Button(button_frame, text="Add", command=self.add_group_tab)
        add_button.pack(side="right", padx=5)
        
        # Add duplicate button
        duplicate_button = ttk.Button(button_frame, text="Duplicate",
                                    command=lambda: self.duplicate_current_group())
        duplicate_button.pack(side="right", padx=5)

        # Add delete button for current group
        delete_button = ttk.Button(button_frame, text="Delete", 
                                 command=lambda: self.delete_current_group())
        delete_button.pack(side="right", padx=5)
    
    def add_group_tab(self, group_name=None):
        """Add a new group tab with build and export settings."""
        if group_name is None:
            group_name = f"Group {len(self.settings) + 1}"
        
        # Create tab frame
        tab_frame = ttk.Frame(self.notebook)
        self.notebook.add(tab_frame, text=group_name)
        
        # Select the newly created tab
        self.notebook.select(tab_frame)
        
        # Create group name entry frame at the top of the tab
        name_frame = ttk.Frame(tab_frame)
        name_frame.pack(fill="x", padx=5, pady=(5, 10))
        
        # Center the name frame
        name_frame.grid_columnconfigure(0, weight=1)
        name_frame.grid_columnconfigure(2, weight=1)
        
        # Group name label and entry
        ttk.Label(name_frame, text="Group Name:").grid(row=0, column=0, padx=5, sticky='w')
        group_name_entry = ttk.Entry(name_frame, width=40)
        group_name_entry.grid(row=0, column=1, columnspan=2, padx=5, sticky='w')
        group_name_entry.insert(0, group_name)
        
        # Save button
        save_button = ttk.Button(name_frame, text="Save Name", 
                               command=lambda: self.save_group_name(group_name, group_name_entry.get()))
        save_button.grid(row=0, column=3, padx=5)
        
        # Create main content frame with two columns
        content_frame = ttk.Frame(tab_frame)
        content_frame.pack(fill="both", expand=True, padx=5, pady=5)
        
        # Configure grid weights for the content frame
        content_frame.grid_columnconfigure(1, weight=1)  # Settings column takes most space
        
        # Create group selection frame (left column)
        selection_frame = ttk.LabelFrame(content_frame, text="Group Selection", padding="5")
        selection_frame.grid(row=0, column=0, columnspan=2, sticky="nsew", padx=(0, 10))
        
        # Create tracking frame at the top
        tracking_frame = ttk.LabelFrame(selection_frame, text="Current Selections", padding="5")
        tracking_frame.pack(fill="x", pady=(0, 5))
        
        # Create text widget for tracking
        tracking_text = tk.Text(tracking_frame, height=4, width=30, wrap=tk.WORD)
        tracking_text.pack(fill="x", padx=5, pady=5)
        tracking_text.config(state='disabled')
        
        # Create frame for dropdown and entry boxes
        selection_options_frame = ttk.Frame(selection_frame)
        selection_options_frame.pack(fill="x", pady=(0, 5))
        
        # Dropdown menu for selection type
        selection_type = tk.StringVar(value="Atoms")
        selection_dropdown = ttk.Combobox(selection_options_frame, 
                                        textvariable=selection_type,
                                        values=["Atom(s) (or Ball(s))", "Residue(s)", "Chain(s)", "Molecule(s)"],
                                        state="readonly",
                                        width=10)
        selection_dropdown.grid(row=1, column=0, padx=(0, 5))
        
        # Create frame for entry boxes and their labels
        entries_frame = ttk.Frame(selection_options_frame)
        entries_frame.grid(row=1, column=1, padx=5)
        
        # Entry boxes frame
        entry_boxes_frame = ttk.Frame(entries_frame)
        entry_boxes_frame.grid(row=0, column=0, padx=5)
        
        # Create frames for each entry box and its label
        start_frame = ttk.Frame(entry_boxes_frame)
        start_frame.grid(row=0, column=0, padx=5)
        
        end_frame = ttk.Frame(entry_boxes_frame)    
        end_frame.grid(row=0, column=2, padx=5)
        
        # Index label above first entry
        ttk.Label(selection_options_frame, text="Index").grid(row=0, column=1, padx=5)

        # Selection label
        ttk.Label(selection_options_frame, text="Selection").grid(row=0, column=0, padx=5)
        
        # Entry box for start value
        start_entry = ttk.Entry(selection_options_frame, width=5)
        start_entry.grid(row=1, column=1, padx=5)
        
        # "to" label
        ttk.Label(selection_options_frame, text="to").grid(row=1, column=2, padx=5)
        
        # Range label above second entry
        ttk.Label(selection_options_frame, text="Range").grid(row=0, column=3, padx=5)
        
        # Entry box for end value
        end_entry = ttk.Entry(selection_options_frame, width=5)
        end_entry.grid(row=1, column=3, padx=5)
        
        # Create button frame for Add and Remove buttons
        button_frame = ttk.Frame(selection_frame)
        button_frame.pack(fill="x", pady=2)
        
        # Add and Remove buttons side by side
        ttk.Button(button_frame, text="Remove", command=lambda: self.delete_group(group_name)).pack(side="left", expand=True, padx=2)
        ttk.Button(button_frame, text="Add", command=lambda: self.add_selection(group_name, selection_type.get(), start_entry.get(), end_entry.get(), tracking_text)).pack(side="right", expand=True, padx=2)
        
        # Create settings container (right column)
        settings_container = ttk.Frame(content_frame)
        settings_container.grid(row=0, column=2, sticky="nsew")
        
        # Create build settings frame
        build_frame = BuildFrame(settings_container, self.gui)
        build_frame.pack(fill="x", pady=(0, 10))
        
        # Create export settings frame
        export_frame = ExportFrame(settings_container, self.gui, group_name_entry)
        export_frame.pack(fill="x")
        
        # Store settings for this group
        self.settings[group_name] = {
            'build_settings': build_frame.get_settings(),  # Store frame reference
            'export_settings': export_frame.get_settings(),  # Store frame reference
            'name_entry': group_name_entry.get(),
            'tracking_text': tracking_text.get(1.0, tk.END),
            'selections': []  # List to store all selections
        }
    
    def add_selection(self, group_name, selection_type, start_val, end_val, tracking_text):
        """Add a new selection to the group and update the tracking display."""
        try:
            start = int(start_val)
            
            # If end_val is empty, treat as single index
            if not end_val.strip():
                end = start
            else:
                end = int(end_val)
            
            if start > end:
                messagebox.showerror("Invalid Range", "Start value must be less than or equal to end value.")
                return
            
            # Create new selection
            new_selection = {
                'type': selection_type,
                'start': start,
                'end': end
            }
            
            # Get current selections
            selections = self.settings[group_name]['selections']
            
            # Check for overlaps
            for selection in selections:
                if (selection['type'] == selection_type and
                    ((start <= selection['end'] and end >= selection['start']) or
                     (selection['start'] <= end and selection['end'] >= start))):
                    # Merge overlapping ranges
                    selection['start'] = min(selection['start'], start)
                    selection['end'] = max(selection['end'], end)
                    break
            else:
                # No overlap found, add new selection
                selections.append(new_selection)
            
            # Update tracking display
            tracking_text.config(state='normal')
            tracking_text.delete(1.0, tk.END)
            
            for selection in selections:
                if selection['start'] == selection['end']:
                    tracking_text.insert(tk.END, f"{selection['type']}: {selection['start']}\n")
                else:
                    tracking_text.insert(tk.END, f"{selection['type']}: {selection['start']}-{selection['end']}\n")
            
            tracking_text.config(state='disabled')

        except ValueError:
            messagebox.showerror("Invalid Input", "Please enter valid numbers for the range.")
    
    def delete_group(self, group_name):
        """Delete a group and its settings."""
        # Don't allow deleting the last group
        if len(self.settings) <= 1:
            messagebox.showwarning("Cannot Delete", "Cannot delete the last group.")
            return
        
        # Ask for confirmation
        if messagebox.askyesno("Confirm Delete", f"Are you sure you want to delete {group_name}?"):
            # Find the tab index
            for i in range(self.notebook.index('end')):
                if self.notebook.tab(i, "text") == group_name:
                    # Remove the tab
                    self.notebook.forget(i)
                    break
            
            # Remove the group settings
            del self.settings[group_name]
            
            # Select the first remaining tab
            self.notebook.select(0)
    
    def save_group_name(self, old_name, new_name):
        """Update the group name and tab text."""
        if new_name and new_name != old_name:
            # Update the tab text
            current_tab = self.notebook.select()
            self.notebook.tab(current_tab, text=new_name)
            
            # Update the settings dictionary
            if old_name in self.settings:
                self.settings[new_name] = self.settings.pop(old_name)
    
    def get_current_group_settings(self):
        """Get the settings for the currently selected group."""
        current_tab = self.notebook.select()
        group_name = self.notebook.tab(current_tab, "text")
        group_data = self.settings[group_name]
        return {
            'name': group_name,
            'build_settings': group_data['build_settings'].get_settings(),
            'export_settings': group_data['export_settings'].get_settings(),
            'selections': group_data['selections']
        }
    
    def get_all_group_settings(self):
        """Get settings for all groups."""
        return {
            group_name: {
                'name': group_name,
                'build_settings': data['build_settings'].get_settings(),
                'export_settings': data['export_settings'].get_settings(),
                'selections': data['selections']
            }
            for group_name, data in self.settings.items()
        }

    def delete_current_group(self):
        """Delete the currently selected group."""
        current_tab = self.notebook.select()
        if current_tab:  # If there is a selected tab
            group_name = self.notebook.tab(current_tab, "text")
            self.delete_group(group_name)

    def duplicate_current_group(self):
        """Duplicate the currently selected group with all its settings."""
        current_tab = self.notebook.select()
        if current_tab:
            # Get current group name
            current_name = self.notebook.tab(current_tab, "text")
            
            # Find next available number for duplicate
            base_name = current_name.split(" (")[0]  # Remove any existing (n) suffix
            counter = 1
            while f"{base_name} ({counter})" in self.settings:
                counter += 1
            new_name = f"{base_name} ({counter})"
            
            # Create new tab
            self.add_group_tab(new_name)
            
            # Get current group's settings
            current_settings = self.settings[current_name]
            new_settings = self.settings[new_name]
            
            # Copy build and export settings using the new methods
            new_settings['build_settings'].copy_settings_from(current_settings['build_settings'])
            new_settings['export_settings'].copy_settings_from(current_settings['export_settings'])
            
            # Copy selections
            new_settings['selections'] = current_settings['selections'].copy()
            
            # Update tracking text
            new_settings['tracking_text'].config(state='normal')
            new_settings['tracking_text'].delete(1.0, tk.END)
            for selection in new_settings['selections']:
                if selection['start'] == selection['end']:
                    new_settings['tracking_text'].insert(tk.END, f"{selection['type']}: {selection['start']}\n")
                else:
                    new_settings['tracking_text'].insert(tk.END, f"{selection['type']}: {selection['start']}-{selection['end']}\n")
            new_settings['tracking_text'].config(state='disabled')
            
            # Select the new tab
            self.notebook.select(self.notebook.index('end')-1)

    def run_current_group(self):
        """Run the currently selected group."""
        current_tab = self.notebook.select()
        if current_tab:
            group_name = self.notebook.tab(current_tab, "text")
            settings = self.get_current_group_settings()
            if self.gui and self.gui.sys:
                # Create a new group in the system
                self.gui.sys.create_group()
                # Apply the settings to the group
                current_group = self.gui.sys.groups[-1]
                current_group.name = settings['name']
                
                # Get settings from the frames
                current_group.build_settings = settings['build_settings']
                current_group.export_settings = settings['export_settings']
                current_group.selections = settings['selections']
                print(f"Running group: {group_name} with settings:", settings)

