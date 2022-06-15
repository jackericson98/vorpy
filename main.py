import tkinter as tk
from tkinter import filedialog, Button

# Initial dialog for vorpy
root = tk.Tk()
root.title('vorpy')


# Grab the file
def grab_file():
    root.withdraw()
    file_path = filedialog.askopenfilename()


head = tk.Label(text="VorPy", font=('Helvetical bold', 40))
greeting = tk.Label(text="Please select an input file(.pdb, .gro, .mol):")
browse = tk.Button(text="browse.. ", command=grab_file)
cancel = tk.Button(text="cancel", command=root.destroy)

head.grid(row=0)
greeting.grid(row=1)
browse.grid(row=2, column=1)
cancel.grid(row=2, column=0)


root.mainloop()

