import tkinter as tk
from tkinter import simpledialog

# Define basic properties of elements (symbol, atomic number, atomic mass, atomic radius)
elements = {
        'H': {'name': 'Hydrogen', 'number': 1, 'mass': 1.008, 'radius': 1.30, 'row': 1, 'column': 1,
              'group': 'Nonmetal'},
        'He': {'name': 'Helium', 'number': 2, 'mass': 4.003, 'radius': 1.40, 'row': 1, 'column': 18,
               'group': 'Noble Gas'},
        'Li': {'name': 'Lithium', 'number': 3, 'mass': 6.941, 'radius': 0.76, 'row': 2, 'column': 1,
               'group': 'Alkali Metal'},
        'Be': {'name': 'Beryllium', 'number': 4, 'mass': 9.012, 'radius': 0.45, 'row': 2, 'column': 2,
               'group': 'Alkaline Earth Metal'},
        'B': {'name': 'Boron', 'number': 5, 'mass': 10.811, 'radius': 1.92, 'row': 2, 'column': 13,
              'group': 'Metalloid'},
        'C': {'name': 'Carbon', 'number': 6, 'mass': 12.011, 'radius': 1.80, 'row': 2, 'column': 14,
              'group': 'Nonmetal'},
        'N': {'name': 'Nitrogen', 'number': 7, 'mass': 14.007, 'radius': 1.60, 'row': 2, 'column': 15,
              'group': 'Nonmetal'},
        'O': {'name': 'Oxygen', 'number': 8, 'mass': 15.999, 'radius': 1.50, 'row': 2, 'column': 16,
              'group': 'Nonmetal'},
        'F': {'name': 'Fluorine', 'number': 9, 'mass': 18.998, 'radius': 1.33, 'row': 2, 'column': 17,
              'group': 'Halogens'},
        'Ne': {'name': 'Neon', 'number': 10, 'mass': 20.180, 'radius': 1.54, 'row': 2, 'column': 18,
               'group': 'Noble Gas'},
        'Na': {'name': 'Sodium', 'number': 11, 'mass': 22.990, 'radius': 1.02, 'row': 3, 'column': 1,
               'group': 'Alkali Metal'},
        'Mg': {'name': 'Magnesium', 'number': 12, 'mass': 24.305, 'radius': 0.72, 'row': 3, 'column': 2,
               'group': 'Alkaline Earth Metal'},
        'Al': {'name': 'Aluminum', 'number': 13, 'mass': 26.982, 'radius': 0.60, 'row': 3, 'column': 13,
               'group': 'Post-transition Metal'},
        'Si': {'name': 'Silicon', 'number': 14, 'mass': 28.086, 'radius': 2.10, 'row': 3, 'column': 14,
               'group': 'Metalloid'},
        'P': {'name': 'Phosphorus', 'number': 15, 'mass': 30.974, 'radius': 1.90, 'row': 3, 'column': 15,
              'group': 'Nonmetal'},
        'S': {'name': 'Sulfur', 'number': 16, 'mass': 32.066, 'radius': 1.90, 'row': 3, 'column': 16,
              'group': 'Nonmetal'},
        'Cl': {'name': 'Chlorine', 'number': 17, 'mass': 35.453, 'radius': 1.81, 'row': 3, 'column': 17,
               'group': 'Halogens'},
        'Ar': {'name': 'Argon', 'number': 18, 'mass': 39.948, 'radius': 1.88, 'row': 3, 'column': 18,
               'group': 'Noble Gas'},
        'K': {'name': 'Potassium', 'number': 19, 'mass': 39.098, 'radius': 1.38, 'row': 4, 'column': 1,
              'group': 'Alkali Metal'},
        'Ca': {'name': 'Calcium', 'number': 20, 'mass': 40.078, 'radius': 1.00, 'row': 4, 'column': 2,
               'group': 'Alkaline Earth Metal'},
        'Ga': {'name': 'Gallium', 'number': 31, 'mass': 69.723, 'radius': 0.62, 'row': 4, 'column': 13,
               'group': 'Post-transition Metal'},
        'Ge': {'name': 'Germanium', 'number': 32, 'mass': 72.631, 'radius': 0.73, 'row': 4, 'column': 14,
               'group': 'Metalloid'},
        'As': {'name': 'Arsenic', 'number': 33, 'mass': 74.922, 'radius': 0.58, 'row': 4, 'column': 15,
               'group': 'Metalloid'},
        'Se': {'name': 'Selenium', 'number': 34, 'mass': 78.971, 'radius': 1.90, 'row': 4, 'column': 16,
               'group': 'Nonmetal'},
        'Br': {'name': 'Bromine', 'number': 35, 'mass': 79.904, 'radius': 1.83, 'row': 4, 'column': 17,
               'group': 'Halogens'},
        'Kr': {'name': 'Krypton', 'number': 36, 'mass': 83.798, 'radius': 2.02, 'row': 4, 'column': 18,
               'group': 'Noble Gas'},
        'Rb': {'name': 'Rubidium', 'number': 37, 'mass': 85.468, 'radius': 1.52, 'row': 5, 'column': 1,
               'group': 'Alkali Metal'},
        'Sr': {'name': 'Strontium', 'number': 38, 'mass': 87.62, 'radius': 1.18, 'row': 5, 'column': 2,
               'group': 'Alkaline Earth Metal'},
        'In': {'name': 'Indium', 'number': 49, 'mass': 114.818, 'radius': 1.93, 'row': 5, 'column': 13,
               'group': 'Post-transition Metal'},
        'Sn': {'name': 'Tin', 'number': 50, 'mass': 118.711, 'radius': 2.17, 'row': 5, 'column': 14,
               'group': 'Post-transition Metal'},
        'Sb': {'name': 'Antimony', 'number': 51, 'mass': 121.760, 'radius': 2.06, 'row': 5, 'column': 15,
               'group': 'Metalloid'},
        'Te': {'name': 'Tellurium', 'number': 52, 'mass': 127.6, 'radius': 2.06, 'row': 5, 'column': 16,
               'group': 'Metalloid'},
        'I': {'name': 'Iodine', 'number': 53, 'mass': 126.904, 'radius': 2.20, 'row': 5, 'column': 17,
              'group': 'Halogens'},
        'Xe': {'name': 'Xenon', 'number': 54, 'mass': 131.293, 'radius': 2.16, 'row': 5, 'column': 18,
               'group': 'Noble Gas'},
        'Cs': {'name': 'Cesium', 'number': 55, 'mass': 132.905, 'radius': 1.67, 'row': 6, 'column': 1,
               'group': 'Alkali Metal'},
        'Ba': {'name': 'Barium', 'number': 56, 'mass': 137.328, 'radius': 1.35, 'row': 6, 'column': 2,
               'group': 'Alkaline Earth Metal'},
        'Tl': {'name': 'Thallium', 'number': 81, 'mass': 204.383, 'radius': 1.96, 'row': 6, 'column': 13,
               'group': 'Post-transition Metal'},
        'Pb': {'name': 'Lead', 'number': 82, 'mass': 207.2, 'radius': 2.02, 'row': 6, 'column': 14,
               'group': 'Post-transition Metal'},
        'Bi': {'name': 'Bismuth', 'number': 83, 'mass': 208.980, 'radius': 2.07, 'row': 6, 'column': 15,
               'group': 'Post-transition Metal'},
        'Po': {'name': 'Polonium', 'number': 84, 'mass': 208.982, 'radius': 1.97, 'row': 6, 'column': 16,
               'group': 'Metalloid'},
        'At': {'name': 'Astatine', 'number': 85, 'mass': 209.987, 'radius': 2.02, 'row': 6, 'column': 17,
               'group': 'Halogens'},
        'Rn': {'name': 'Radon', 'number': 86, 'mass': 222.018, 'radius': 2.20, 'row': 6, 'column': 18,
               'group': 'Noble Gas'},
        'Fr': {'name': 'Francium', 'number': 87, 'mass': 223.020, 'radius': 3.48, 'row': 7, 'column': 1,
               'group': 'Alkali Metal'},
        'Ra': {'name': 'Radium', 'number': 88, 'mass': 226.025, 'radius': 2.83, 'row': 7, 'column': 2,
               'group': 'Alkaline Earth Metal'},
        'Zn': {'name': 'Zinc', 'number': 30, 'mass': 65.38, 'radius': 1.39, 'row': 4, 'column': 12,
               'group': 'Transition Metal'}
    }


class ElementDialog(tk.Toplevel):
    def __init__(self, parent, element, callback):
        super().__init__(parent)
        self.element = element
        self.callback = callback
        self.title(f"Edit Properties - {element['name']}")

        tk.Label(self, text=f"Mass ({element['mass']} u):").grid(row=0, column=0)
        self.mass_entry = tk.Entry(self)
        self.mass_entry.insert(0, str(element['mass']))
        self.mass_entry.grid(row=0, column=1)

        tk.Label(self, text=f"Radius ({element['radius']} Å):").grid(row=1, column=0)
        self.radius_entry = tk.Entry(self)
        self.radius_entry.insert(0, str(element['radius']))
        self.radius_entry.grid(row=1, column=1)

        self.apply_button = tk.Button(self, text="Apply", command=self.apply)
        self.apply_button.grid(row=2, column=0)

        self.cancel_button = tk.Button(self, text="Cancel", command=self.destroy)
        self.cancel_button.grid(row=2, column=1)

    def apply(self):
        try:
            new_mass = float(self.mass_entry.get())
            new_radius = float(self.radius_entry.get())
            self.callback(new_mass, new_radius)
        except ValueError:
            pass  # Handle incorrect input gracefully
        self.destroy()


def update_properties(element, button, root):
    def apply_changes(new_mass, new_radius):
        element['mass'] = new_mass
        element['radius'] = new_radius
        button.config(text=f"{element['name']}\n{element['mass']} u\n{element['radius']} Å")

    ElementDialog(root, element, apply_changes)


def create_button(root, element, color, buttons):
    button = tk.Button(root, text=f"{element['name']}\n{element['mass']} u\n{element['radius']} Å",
                       command=lambda: update_properties(element, button, root), padx=10, pady=10, bg=color)
    return button


def periodic_table():
    root = tk.Tk()
    root.title("Editable Periodic Table")

    color_scheme = {
        'Nonmetal': '#FFD700',  # Gold
        'Noble Gas': '#FFC0CB',  # Pink
        'Alkali Metal': '#F08080',  # Light Coral
        'Alkaline Earth Metal': '#00BFFF',  # Deep Sky Blue
        'Metalloid': '#ADFF2F',  # Green Yellow
        'Halogens': '#FFA500',  # Orange
        'Post-transition Metal': '#20B2AA',  # Light Sea Green
        'Transition Metal': '#B0C4DE',  # Light Steel Blue
    }

    buttons = {}
    max_row = 0
    for symbol, element in elements.items():
        btn = create_button(root, element, color_scheme[element['group']], buttons)
        buttons[element['name']] = btn
        row, col = element['row'], element['column']
        btn.grid(row=row, column=col, sticky='nsew')
        max_row = max(max_row, element['row'])

    for i in range(18):
        root.grid_columnconfigure(i, weight=1)
        root.grid_rowconfigure(i, weight=1)

    tk.Button(root, text='Cancel', command=root.destroy).grid(row=max_row + 1, column=17)
    tk.Button(root, text='Apply', command=root.destroy).grid(row=max_row + 1, column=18)

    root.mainloop()


if __name__ == '__main__':
    periodic_table()
