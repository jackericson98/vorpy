import tkinter as tk
from tkinter import ttk

data = False


def settings_gui():
    # Function to collect values and print dictionary
    global data

    def apply_values():
        global data
        data = {
            "bubble size": float(avg_bubble_size_var.get()),
            "bubble sd": float(std_deviation_var.get()),
            "bubble num": int(num_bubbles_var.get()),
            "bubble density": float(density_var.get()),
            "open cell": open_cell_var.get(),
            "distribution": distribution_var.get().lower()
        }
        root.destroy()
        return data  # Or replace this with return data if using within another function or script

    # Function to handle cancel
    def cancel():
        root.destroy()

    # Main window
    root = tk.Tk()
    root.title("Foam Settings")

    # Variables for storing input
    avg_bubble_size_var = tk.StringVar(value='1.0')
    std_deviation_var = tk.StringVar(value='0.1')
    num_bubbles_var = tk.StringVar(value='100')
    density_var = tk.StringVar(value='0.25')
    open_cell_var = tk.BooleanVar(value=False)
    distribution_var = tk.StringVar(value='Lognormal')

    # Average Bubble Size Entry
    tk.Label(root, text="Avg Bubble Size").grid(row=0, column=0)
    tk.Entry(root, textvariable=avg_bubble_size_var).grid(row=0, column=1)
    tk.Label(root, text='(0.001 - 1000 \u212B)').grid(row=0, column=2)

    # Standard Deviation Entry
    tk.Label(root, text="Standard Deviation").grid(row=1, column=0)
    tk.Entry(root, textvariable=std_deviation_var).grid(row=1, column=1)
    tk.Label(root, text='(0 - 1000)').grid(row=1, column=2)

    # Number of bubbles entry
    tk.Label(root, text="Number of Bubbles").grid(row=2, column=0)
    tk.Entry(root, textvariable=num_bubbles_var).grid(row=2, column=1)
    tk.Label(root, text='(10 - 100,000)').grid(row=2, column=2)

    #
    tk.Label(root, text="Density").grid(row=3, column=0)
    tk.Entry(root, textvariable=density_var).grid(row=3, column=1)
    tk.Label(root, text='(0.001 - 0.6)').grid(row=3, column=2)

    # Checkbox
    tk.Checkbutton(root, text="Open Cell?", variable=open_cell_var).grid(row=4, column=0, columnspan=3)

    # Dropdown
    tk.Label(root, text="Distribution").grid(row=5, column=0)
    distribution_menu = ttk.Combobox(root, textvariable=distribution_var, values=["Lognormal", "Gamma", "Real"])
    distribution_menu.grid(row=5, column=1)
    distribution_menu.current(0)  # Set the default value

    # Buttons
    tk.Button(root, text="Apply", command=apply_values).grid(row=6, column=1)
    tk.Button(root, text="Cancel", command=cancel).grid(row=6, column=2)

    # Run the GUI
    root.mainloop()

    if not data:
        data = {'bubble size': 1, 'bubble sd': 0.1, 'bubble num': 100, 'bubble density': 0.25,
                'open cell': False, 'distribution': 'lognormal'}
    return data


