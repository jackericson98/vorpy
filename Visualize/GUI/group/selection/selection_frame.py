import tkinter as tk
from tkinter import ttk

class SelectionFrame(ttk.LabelFrame):
    def __init__(self, parent, gui, group_name_entry):
        super().__init__(parent, text="Selection", padding=(5, 5))
        self.parent = parent
        self.gui = gui
        self.group_name_entry = group_name_entry

        self.selections = {'balls': None, 'residues': None, 'chains': None, 'molecules': None}
        self.tracking = {'balls': '', 'residues': '', 'chains': '', 'molecules': ''}
        self.undo_stack = []
        self.selection_labels = {}  # Store references to the selection labels

        self._create_widgets()

    def _create_widgets(self):
        # # Create a centered and underlined label for the selection frame
        # label = ttk.Label(self, text="Choose Group Members", font=('TkDefaultFont', 10, 'underline'))
        # label.grid(row=0, column=0, columnspan=2, sticky='ew', pady=5)

        # Create frame for dropdown and entry boxes
        selection_options_frame = ttk.Frame(self)
        selection_options_frame.grid(row=1, column=0, columnspan=2, sticky='ew', pady=(0, 5))

        # Create a label for the added selections
        ttk.Label(selection_options_frame, text="Add/Remove Selections", font=('TkDefaultFont', 10, 'underline')).grid(row=0, column=0, columnspan=4, sticky='n', pady=(5, 10))

        # Selection label
        ttk.Label(selection_options_frame, text="Selection").grid(row=1, column=0, padx=5)

        # Index label above first entry
        ttk.Label(selection_options_frame, text="Index").grid(row=1, column=1, padx=5)
                
        # Or indicator for user to know that they can either do a single index or a range
        ttk.Label(selection_options_frame, text="or").grid(row=1, column=2, padx=5)
        
        # Range label above second entry
        ttk.Label(selection_options_frame, text="Range").grid(row=1, column=3, padx=5)

        # Dropdown menu for selection type
        selection_type = tk.StringVar(value="Atoms/Balls")
        selection_dropdown = ttk.Combobox(selection_options_frame, 
                                        textvariable=selection_type,
                                        values=["Atoms/Balls", "Residues", "Chains", "Molecules"],
                                        state="readonly",
                                        width=20, justify="center")
        selection_dropdown.grid(row=2, column=0, padx=(5, 10))
        
        # Create frame for entry boxes and their labels
        entries_frame = ttk.Frame(selection_options_frame)
        entries_frame.grid(row=2, column=1, padx=5)

        # Entry box for start value
        start_entry = ttk.Entry(selection_options_frame, width=5)
        start_entry.grid(row=2, column=1, padx=5)
        
        # "to" label
        ttk.Label(selection_options_frame, text="to").grid(row=2, column=2, padx=5)
        
        # Entry box for end value
        end_entry = ttk.Entry(selection_options_frame, width=5)
        end_entry.grid(row=2, column=3, padx=5)
                
        # Create button frame for Add and Remove buttons
        button_frame = ttk.Frame(self)
        button_frame.grid(row=2, column=0, columnspan=2, sticky='ew', pady=2)
        
        # Add and Remove buttons side by side
        ttk.Button(button_frame, text="Add", command=lambda: self.add_selection(selection_type.get(), start_entry.get(), end_entry.get(), undo_command=False)).grid(row=0, column=3)
        ttk.Button(button_frame, text="Remove", command=lambda: self.delete_selection(selection_type.get(), start_entry.get(), end_entry.get(), undo_command=False)).grid(row=0, column=2)
        # ttk.Button(button_frame, text="Undo", command=lambda: self.undo_selections()).grid(row=0, column=1, padx=2)
        ttk.Button(button_frame, text="Clear", command=lambda: self.clear_selections()).grid(row=0, column=0)

        # Configure button frame columns
        button_frame.grid_columnconfigure((0,1,2,3), weight=1)

        # Create selection display frames
        selections_container = ttk.Frame(self)
        selections_container.grid(row=3, column=0, columnspan=2, sticky='nsew', padx=5, pady=5)
        selections_container.grid_columnconfigure((0, 1), weight=1)
        selections_container.grid_rowconfigure(0, weight=1)  # Allow vertical expansion

        # Create a label for the added selections
        ttk.Label(selections_container, text="Selection Tracker", font=('TkDefaultFont', 10, 'underline')).grid(row=0, column=0, columnspan=2, padx=5, pady=5, sticky='n')

        # Create label frames for each selection type
        selection_types = ['Atoms/Balls', 'Residues', 'Chains', 'Molecules']
        
        # Create two columns of selection frames
        for idx, sel_type in enumerate(selection_types):
            # Calculate row and column for 2x2 grid
            row = idx // 2 + 1
            col = idx % 2
            
            # Create frame for this selection type
            type_frame = ttk.LabelFrame(selections_container, text=selection_types[idx])
            type_frame.grid(row=row, column=col, sticky='nsew', pady=1, padx=2)
            type_frame.grid_columnconfigure(row, weight=1)
            type_frame.grid_rowconfigure(col, weight=1)
            
            # Label for selections
            selection_label = ttk.Label(type_frame, text="None", wraplength=150, justify="center")
            selection_label.grid(row=0, column=0, padx=5, pady=1, sticky='w')
            
            # Store reference to label
            key = sel_type.lower().replace('/', '_')
            self.selection_labels[key] = selection_label

    def create_selection_string(self, list_of_ndxs):
        """Takes a list of indexes and returns a string of the form 'index1-index2, index3-index4, etc.'"""
        # Copy the list of indexes
        sorted_list_of_ndxs = list_of_ndxs.copy()
        # Sort the list of indexes
        sorted_list_of_ndxs.sort()
        if len(sorted_list_of_ndxs) == 0:
            return ''
        if len(sorted_list_of_ndxs) == 1:
            return str(sorted_list_of_ndxs[0])
        # Create a list of ranges
        ranges = [[sorted_list_of_ndxs[0]]]
        for i in range(1, len(sorted_list_of_ndxs)):
            if sorted_list_of_ndxs[i] == ranges[-1][-1] + 1:
                ranges[-1].append(sorted_list_of_ndxs[i])
            else:
                ranges.append([sorted_list_of_ndxs[i]])
        return ', '.join([f'{r[0]}-{r[-1]}' if len(r) > 1 else str(r[0]) for r in ranges])
    
    def update_tracking_text(self):
        """Updates the selection labels to show the current selections"""
        # Update each selection type's label
        self.selection_labels['atoms_balls'].config(text=self.tracking['balls'] or "None")
        self.selection_labels['residues'].config(text=self.tracking['residues'] or "None")
        self.selection_labels['chains'].config(text=self.tracking['chains'] or "None")
        self.selection_labels['molecules'].config(text=self.tracking['molecules'] or "None")

    def add_selection(self, selection_type, start, end=None, undo_command=False):
        """Goes through the selections added, checks the self.selections dictionary and only adds the selected indexes if they are not already in the 
        dictionary. Once added, the text showing the current selections will be updated to show any index ranges, minimizing the amount of text listed"""
        if start == '' or start is None or not start.isdigit():
            return
        # Check if the start and end indices are strings
        if isinstance(start, str):
            start = int(start)
            if start < 0:
                return
        if end is not None and end != '' and isinstance(end, str):
            end = int(end)
            if end < 0:
                return
        # Get a list of new indexes to add to whatever the selection is
        if end is None or end == '':
            new_ndxs = [start]
        else:
            new_ndxs = list(range(start, end+1))
        # Add the given selection to the self.selections dictionary and make sure duplicates are filtered out
        if selection_type == 'Atoms/Balls':
            if self.selections['balls'] is None:
                self.selections['balls'] = new_ndxs
            else:
                self.selections['balls'].extend(new_ndxs)
                self.selections['balls'] = list(set(self.selections['balls']))
            self.tracking['balls'] = self.create_selection_string(self.selections['balls'])
        elif selection_type == 'Residues':
            if self.selections['residues'] is None:
                self.selections['residues'] = new_ndxs
            else:
                self.selections['residues'].extend(new_ndxs)
                self.selections['residues'] = list(set(self.selections['residues']))
            self.tracking['residues'] = self.create_selection_string(self.selections['residues'])
        elif selection_type == 'Chains':
            if self.selections['chains'] is None:
                self.selections['chains'] = new_ndxs
            else:
                self.selections['chains'].extend(new_ndxs)
                self.selections['chains'] = list(set(self.selections['chains']))
            self.tracking['chains'] = self.create_selection_string(self.selections['chains'])
        elif selection_type == 'Molecules':
            if self.selections['molecules'] is None:
                self.selections['molecules'] = new_ndxs
            else:
                self.selections['molecules'].extend(new_ndxs)
                self.selections['molecules'] = list(set(self.selections['molecules']))
            self.tracking['molecules'] = self.create_selection_string(self.selections['molecules'])
        
        # Update the tracking text
        self.update_tracking_text()
        # Add the current state to the undo stack
        if undo_command:
            self.undo_stack.append({
                'selections': selection_type,
                'start': start,
                'end': end,
                'action': 'add'
            })
    
    def delete_selection(self, selection_type, start, end=None, undo_command=False):
        """Deletes a selection from the self.selections dictionary"""
        # Check if the start and end indices are strings
        if isinstance(start, str):
            start = int(start)
        if end is not None and end != '' and isinstance(end, str):
            end = int(end)
        # Get a list of indexes to remove
        if end is None or end == '':
            removal_ndxs = [start]
        else:
            removal_ndxs = list(range(start, end+1))
        if selection_type == 'Atoms/Balls':
            self.selections['balls'] = [ndx for ndx in self.selections['balls'] if ndx not in removal_ndxs]
            self.tracking['balls'] = self.create_selection_string(self.selections['balls'])
        elif selection_type == 'Residues':
            self.selections['residues'] = [ndx for ndx in self.selections['residues'] if ndx not in removal_ndxs]
            self.tracking['residues'] = self.create_selection_string(self.selections['residues'])
        elif selection_type == 'Chains':
            self.selections['chains'] = [ndx for ndx in self.selections['chains'] if ndx not in removal_ndxs]
            self.tracking['chains'] = self.create_selection_string(self.selections['chains'])
        elif selection_type == 'Molecules':
            self.selections['molecules'] = [ndx for ndx in self.selections['molecules'] if ndx not in removal_ndxs]
            self.tracking['molecules'] = self.create_selection_string(self.selections['molecules'])
        
        # Update the tracking text
        self.update_tracking_text()
        # Add the current state to the undo stack
        if not undo_command:
            self.undo_stack.append({
                'selections': selection_type,
                'start': start,
                'end': end,
                'action': 'remove'
            })

    def clear_selections(self):
        """Clears all the selections in the self.selections dictionary"""
        self.undo_stack.append({
            'selections': self.selections,
            'start': None,
            'end': None,
            'action': 'clear',
            'tracking': self.tracking
        })
        self.selections = {'balls': [], 'residues': [], 'chains': [], 'molecules': []}
        self.tracking = {'balls': '', 'residues': '', 'chains': '', 'molecules': ''}
        self.update_tracking_text()
        

    def undo_selections(self):
        """Undoes the last action in the undo stack"""
        if self.undo_stack:
            last_action = self.undo_stack.pop()
            if last_action['action'] == 'add':
                self.delete_selection(last_action['selections'], last_action['start'], last_action['end'], undo_command=True)
            elif last_action['action'] == 'remove':
                self.add_selection(last_action['selections'], last_action['start'], last_action['end'], undo_command=True)
            elif last_action['action'] == 'clear':
                self.selections = last_action['selections']
                self.tracking = last_action['tracking']
            self.update_tracking_text()



