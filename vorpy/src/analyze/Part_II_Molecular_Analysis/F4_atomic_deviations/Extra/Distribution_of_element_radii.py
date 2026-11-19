import os
from matplotlib import pyplot as plt
import tkinter as tk
from tkinter import filedialog

import os
import sys

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.system.system import System
from vorpy.src.analyze.tools.plot_templates.histogram import histogram


root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)
my_pdb = filedialog.askopenfilename()

# Load the system into the system using the system object
os.chdir('../../../../..')
my_sys = System(my_pdb)


# Check the radii
rads = my_sys.balls['rad'].tolist()
print(rads)

# Plot the distribution
histogram(rads)
plt.show()


