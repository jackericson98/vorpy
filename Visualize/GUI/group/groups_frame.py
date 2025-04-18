import tkinter as tk
from tkinter import ttk, messagebox
from Visualize.GUI.group.build.build_frame import BuildFrame
from Visualize.GUI.group.export.export_frame import ExportFrame
from Visualize.GUI.group.selection.selection_frame import SelectionFrame


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
        self.selection_frame = SelectionFrame(content_frame, self.gui, group_name_entry)
        self.selection_frame.grid(row=0, column=0, sticky="nsew", padx=5)

        # Create settings container (right column)
        settings_container = ttk.Frame(content_frame)
        settings_container.grid(row=0, column=2, sticky="nsew")
        
        # Create build settings frame
        self.build_frame = BuildFrame(settings_container, self.gui)
        self.build_frame.pack(fill="x")
        
        # Create export settings frame
        self.export_frame = ExportFrame(settings_container, self.gui, group_name_entry)
        self.export_frame.pack(fill="x")
        
        # Store settings for this group
        self.settings[group_name] = {
            'build_settings': self.build_frame,  # Store frame reference
            'export_settings': self.export_frame,  # Store frame reference
            'name_entry': group_name_entry,
            'selections': self.selection_frame  # List to store all selections
        }
        
        # Create a button frame for the buttons
        button_frame = ttk.Frame(content_frame)
        button_frame.grid(row=1, column=1, columnspan=2, sticky="ew", padx=5, pady=5)

        # Create a run button for the current group
        run_button = ttk.Button(button_frame, text="Run", command=self.run_current_group)
        run_button.pack(side="right", padx=5)

        # Add button for current group        
        add_button = ttk.Button(button_frame, text="Add", command=self.add_group_tab)
        add_button.pack(side="right", padx=5)

        # Add delete button for current group
        delete_button = ttk.Button(button_frame, text="Delete", 
                                 command=lambda: self.delete_current_group())
        delete_button.pack(side="right", padx=5)
    
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
        return {
            'name': group_name,
            'build_settings': self.settings[group_name]['build_settings'].get_settings(),
            'export_settings': self.settings[group_name]['export_settings'].get_settings(),
            'selections': self.settings[group_name]['selections'].selections
        }
    
    def get_all_group_settings(self):
        """Get settings for all groups."""
        return {
            group_name: {
                'name': group_name,
                'build_settings': data['build_settings'].get_settings(),
                'export_settings': data['export_settings'].get_settings(),
                'selections': data['selections'].selections
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
            self.gui.run_group(group_name)
