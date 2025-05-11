import os

import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
import tkinter as tk
from tkinter import filedialog
from vorpy.src.system.system import System
from System.Group.group import Group
#
# # Points and labels
# points = [1, 2, 4, 5, 6, 6.5, 7, 7.1, 8.5, 9.2, 10]
# labels = ['one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine', 'ten', 'eleven']
#
# # Plotting the points
# plt.figure(figsize=(8, 4))  # Adjust the figure size as needed
# plt.plot(points, [0] * len(points), 'ro')  # 'ro' means red color, circle markers
#
# pad = 0.1
# # Adding labels for each point with diagonal rotation
# for point, label in zip(points, labels):
#     text = plt.text(point + pad, 0.01, f'{label}', ha='right', va='bottom', rotation=90)
#
#     # # Adjust the starting position of the labels
#     # trans = Affine2D().translate(-12, 0) + plt.gca().transData
#     # text.set_transform(trans)
#
# # Drawing lines from points to labels
# for point, label in zip(points, labels):
#     plt.plot([point, point], [0, 0.01], 'k-', linewidth=0.7)
#
# # Setting up plot limits and removing y-axis
# plt.xlim(min(points) - 1, max(points) + 1)
# plt.ylim(-0.1, 0.1)
# plt.yticks([])  # Remove y-axis ticks and labels
#
# # Display the plot
# plt.show()


if __name__ == '__main__':
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    folder = filedialog.askdirectory()
    systems = []
    for root, dir, files in os.walk(folder):
        for file in files:
           if file[-3:] == 'pdb':
              systems.append(System(file=folder + '/' + file))
    num_atoms = [len(_.atoms) for _ in systems]
    systems = [x for _, x in sorted(zip(num_atoms, systems))]
    print('name', 'atoms', 'residues', 'chains', 'non sol atoms')
    for my_sys in systems:
        my_group = Group(sys=my_sys, residues=my_sys.residues)
        print(my_sys.name, len(my_sys.atoms), len(my_sys.residues), len(my_sys.chains), len(my_group.atoms))

