import tkinter as tk
from tkinter import ttk
from vorpy.src.GUI.system.radii_adjustments.element_classifications import elements, special_radii


class ElementDialog(tk.Toplevel):
    def __init__(self, parent, element, callback):
        super().__init__(parent)
        self.element = element
        self.callback = callback
        self.title(f"Edit Properties - {element['name']}")
        
        # Configure window
        self.resizable(False, False)
        self.transient(parent)
        self.grab_set()
        
        # Create main frame with padding
        main_frame = ttk.Frame(self, padding="10")
        main_frame.pack(fill="both", expand=True)
        
        # Element Info Frame
        info_frame = ttk.LabelFrame(main_frame, text="Element Information", padding="5")
        info_frame.pack(fill="x", padx=5, pady=5)
        
        # Element details in a grid
        ttk.Label(info_frame, text=f"{element['name']}", font=('Arial', 12, 'bold')).grid(row=0, column=0, padx=5,
                                                                                          pady=2)
        ttk.Label(info_frame, text=f"Atomic Number: {element['number']}").grid(row=1, column=0, sticky="w", padx=5,
                                                                               pady=2)
        ttk.Label(info_frame, text=f"Group: {element['group']}").grid(row=2, column=0, sticky="w", padx=5, pady=2)
        
        # Properties Frame
        props_frame = ttk.LabelFrame(main_frame, text="Edit Properties", padding="5")
        props_frame.pack(fill="x", padx=5, pady=5)
        
        # Mass entry
        ttk.Label(props_frame, text=f"Mass ({element['mass']} u):").grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.mass_entry = ttk.Entry(props_frame, width=15)
        self.mass_entry.insert(0, str(element['mass']))
        self.mass_entry.grid(row=0, column=1, sticky="w", padx=5, pady=2)
        
        # Radius entry
        ttk.Label(props_frame, text=f"Radius ({element['radius']} Å):").grid(row=1, column=0, sticky="w", padx=5,
                                                                             pady=2)
        self.radius_entry = ttk.Entry(props_frame, width=15)
        self.radius_entry.insert(0, str(element['radius']))
        self.radius_entry.grid(row=1, column=1, sticky="w", padx=5, pady=2)
        
        # Button frame
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill="x", pady=(10, 0))
        
        # Buttons
        self.apply_button = ttk.Button(button_frame, text="Apply", command=self.apply)
        self.apply_button.pack(side="right", padx=5)
        
        self.cancel_button = ttk.Button(button_frame, text="Cancel", command=self.destroy)
        self.cancel_button.pack(side="right", padx=5)
        
        # Center the window on the parent
        self.update_idletasks()
        width = self.winfo_width()
        height = self.winfo_height()
        x = (self.winfo_screenwidth() // 2) - (width // 2)
        y = (self.winfo_screenheight() // 2) - (height // 2)
        self.geometry(f"{width}x{height}+{x}+{y}")

    def apply(self):
        try:
            new_mass = float(self.mass_entry.get())
            new_radius = float(self.radius_entry.get())
            self.callback(new_mass, new_radius)
        except ValueError:
            pass  # Handle incorrect input gracefully
        self.destroy()


class PeriodicTableGUI(tk.Toplevel):
    def __init__(self, parent):
        super().__init__(parent)
        
        # Configure window
        self.title("Editable Periodic Table")
        self.parent = parent
        
        # Make window modal
        self.transient(parent)
        self.grab_set()
        
        self.buttons = {}
        self.color_scheme = {
            'Nonmetal': '#FFD700',  # Gold
            'Noble Gas': '#FFC0CB',  # Pink
            'Alkali Metal': '#F08080',  # Light Coral
            'Alkaline Earth Metal': '#00BFFF',  # Deep Sky Blue
            'Metalloid': '#ADFF2F',  # Green Yellow
            'Halogens': '#FFA500',  # Orange
            'Post-transition Metal': '#20B2AA',  # Light Sea Green
            'Transition Metal': '#B0C4DE',  # Light Steel Blue
        }

        # Create main container frame
        main_container = ttk.Frame(self)
        main_container.pack(fill='both', expand=True, padx=10, pady=10)

        # Create periodic table frame
        self.periodic_table_frame = ttk.Frame(main_container)
        self.periodic_table_frame.pack(fill='both', expand=True)

        # Create button frame at bottom
        self.button_frame = ttk.Frame(main_container)
        self.button_frame.pack(fill='x')
        
        # Add OK and Cancel buttons
        ttk.Button(self.button_frame, text="Save", command=self.apply_and_close).pack(side='right', padx=5)
        ttk.Button(self.button_frame, text="Cancel", command=self.destroy).pack(side='right', padx=5)

        self._create_widgets()

    def apply_and_close(self):
        """Apply changes and close the window."""
        # Here you can add code to save/apply any changes made to the periodic table
        self.destroy()

    def _update_properties(self, element, frame):
        def apply_changes(new_mass, new_radius):
            # Create or update the element changes list in the parent GUI
            if not hasattr(self.parent, 'element_changes'):
                self.parent.element_changes = []
            
            # Only add to the list if values have changed
            if new_mass != element['mass'] or new_radius != element['radius']:
                # Update the parent GUI's radii_changes list
                if not hasattr(self.parent, 'radii_changes'):
                    self.parent.radii_changes = []
                self.parent.radii_changes.append({
                    'element': element['name'],
                    'radius': new_radius,
                    'mass': new_mass
                })
                
            # Update the element properties
            element['mass'] = new_mass
            element['radius'] = new_radius
            
            # Update the display
            for widget in frame.winfo_children():
                if isinstance(widget, tk.Frame):
                    for label in widget.winfo_children():
                        if label.winfo_x() < frame.winfo_width() // 2:
                            label.config(text=f"{new_radius} Å")
                        else:
                            label.config(text=f"{new_mass} u")

        ElementDialog(self, element, apply_changes)

    def _create_element_frame(self, symbol, element):
        frame = tk.Frame(self.periodic_table_frame, bg=self.color_scheme[element['group']])
        frame.grid(row=element['row'], column=element['column'], sticky='nsew', padx=1, pady=1)
        
        atomic_num = tk.Label(frame, text=str(element['number']), 
                              font=('Arial', 8), bg=self.color_scheme[element['group']])
        atomic_num.pack(anchor='nw')
        
        symbol_label = tk.Label(frame, text=symbol, 
                                font=('Arial', 15, 'bold'), bg=self.color_scheme[element['group']])
        symbol_label.pack()
        
        name_label = tk.Label(frame, text=element['name'], 
                              font=('Arial', 8), bg=self.color_scheme[element['group']])
        name_label.pack()
        
        bottom_frame = tk.Frame(frame, bg=self.color_scheme[element['group']])
        bottom_frame.pack(fill='x')
        
        radius_label = tk.Label(bottom_frame, text=f"{round(element['radius'], 2)} Å", 
                                font=('Arial', 8), bg=self.color_scheme[element['group']])
        radius_label.pack(side='left')
        
        mass_label = tk.Label(bottom_frame, text=f"{round(element['mass'], 2)} u", 
                              font=('Arial', 8), bg=self.color_scheme[element['group']])
        mass_label.pack(side='right')
        
        frame.bind('<Button-1>', lambda e, el=element, f=frame: self._update_properties(el, f))
        for child in frame.winfo_children():
            child.bind('<Button-1>', lambda e, el=element, f=frame: self._update_properties(el, f))
        
        self.buttons[element['name']] = frame
        return frame

    def _create_widgets(self):
        max_row = 0
        for symbol, element in elements.items():
            max_row = max(max_row, element['row'])
            self._create_element_frame(symbol, element)
        # Create a frame for special radii information below the periodic table
        special_radii_frame = tk.LabelFrame(self.periodic_table_frame, text="Special Radii", 
                                          font=('Arial', 10, 'bold'))
        special_radii_frame.grid(row=max_row + 1, column=0, columnspan=19, sticky='nsew', padx=5, pady=(15, 5))
                
        
        def update_atom_dropdown(*args):
            selected_residue = residue_var.get()
            if selected_residue in special_radii:
                atom_dropdown['values'] = list(special_radii[selected_residue].keys())
                if atom_dropdown['values']:
                    atom_dropdown.set(atom_dropdown['values'][0])
                else:
                    atom_dropdown.set('')
        labels_frame = tk.Frame(special_radii_frame)
        labels_frame.pack(fill='x', padx=10, pady=5)

        labels_frame.grid_columnconfigure((0, 1, 2, 3, 4), weight=1)
        # Residue label
        ttk.Label(labels_frame, text="Residue:").grid(row=0, column=0, padx=5, pady=5)
        # Atom label
        ttk.Label(labels_frame, text="Atom:").grid(row=0, column=1, padx=5, pady=5)
        # Radius label
        ttk.Label(labels_frame, text="Radius:").grid(row=0, column=2, padx=5, pady=5)
        # Mass label
        ttk.Label(labels_frame, text="Mass:").grid(row=0, column=3, padx=5, pady=5)
        
        
        # Residue dropdown
        residue_var = tk.StringVar()
        residue_dropdown = ttk.Combobox(labels_frame, textvariable=residue_var, state="readonly", width=15)
        residue_dropdown['values'] = list(special_radii.keys())
        residue_dropdown.grid(row=1, column=0, padx=5, pady=5)
        
        # Atom dropdown
        atom_var = tk.StringVar()
        atom_dropdown = ttk.Combobox(labels_frame, textvariable=atom_var, state="readonly", width=15)
        atom_dropdown.grid(row=1, column=1, padx=5, pady=5)
        
        # Radius entry
        radius_var = tk.StringVar()
        radius_entry = tk.Entry(labels_frame, textvariable=radius_var, width=10)
        radius_entry.grid(row=1, column=2, padx=5, pady=5)
        
        # Mass entry
        mass_var = tk.StringVar()
        mass_entry = tk.Entry(labels_frame, textvariable=mass_var, width=10)
        mass_entry.grid(row=1, column=3, padx=5, pady=5)
        
        # Change button
        change_button = ttk.Button(labels_frame, text="Change", width=10)
        change_button.grid(row=1, column=4, padx=5, pady=5)
        
        def update_radius_and_mass(*args):
            selected_residue = residue_var.get()
            selected_atom = atom_var.get()
            if selected_residue in special_radii and selected_atom in special_radii[selected_residue]:
                radius_var.set(str(special_radii[selected_residue][selected_atom]['radius']))
                mass_var.set(str(special_radii[selected_residue][selected_atom]['mass']))

        # Bind the residue dropdown to update the atom dropdown
        residue_var.trace_add('write', update_atom_dropdown)
        # Bind the atom dropdown to update radius and mass entries
        atom_var.trace_add('write', update_radius_and_mass)
        
        # Configure grid weights for the special radii frame
        self.periodic_table_frame.grid_rowconfigure(max_row + 1, weight=1)
        for i in range(18):
            self.periodic_table_frame.grid_columnconfigure(i, weight=1)


if __name__ == "__main__":
    root = tk.Tk()
    root.withdraw()  # Hide the main window
    app = PeriodicTableGUI(root)
    root.mainloop()
