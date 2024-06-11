from tkinterdnd2 import DND_FILES, TkinterDnD
import tkinter as tk
from tkinter import ttk


setting = False


def my_GUI():
    global settings

    def drop(event):
        file_path.set(event.data)  # Set the file path in the Label

    def update_variables():
        global settings
        settings = {
            # Files
            'base file': file_path.get(),
            # Settings
            # 'mv':
        }
        root.destroy()
        return settings

    root = TkinterDnD.Tk()
    root.title('Vorpy')

    # Variable to store the path of the dragged file
    file_path = tk.StringVar()
    file_path.set('Drag and drop your file here.')

    # Create a label to display the file path or instructions
    label = tk.Label(root, textvariable=file_path, width=80, height=4, relief='sunken')
    label.grid(row=0, column=0, columnspan=5, padx=10, pady=10)

    # Configure the label to accept dragged files
    label.drop_target_register(DND_FILES)
    label.dnd_bind('<<Drop>>', drop)

    # Add a run button
    tk.Button(root, command=update_variables, text='Run').grid(row=4, column=4, padx=10, pady=10)

    # Add a cancel function and button
    def cancel():
        root.destroy()

    tk.Button(root, command=cancel, text='Cancel').grid(row=4, column=3, padx=10, pady=10)

    # Network Type
    net_type = tk.StringVar(value='Additively Weighted')
    (ttk.Combobox(root, textvariable=net_type, values=['Additively Weighted', 'Power', 'Primitive'])
     .grid(row=2, column=3, padx=10, pady=10))


    # Maximum Vertex Size
    max_vert = tk.StringVar(value='40')
    tk.Entry(root, textvariable=max_vert).grid(row=1, column=0, padx=10, pady=10)

    # Surface Resolution
    surf_res = tk.StringVar(value='0.2')
    tk.Entry(root, textvariable=surf_res).grid(row=2, column=0, padx=10, pady=10)

    # Export type
    exp_type = tk.StringVar(value='Large')
    ttk.Combobox(root, textvariable=exp_type, values=['Logs', 'Small', 'Medium', 'Large', 'All']).grid(padx=10, pady=10)


    # End the main loop
    root.mainloop()
    return settings


if __name__ == '__main__':
    my_GUI()