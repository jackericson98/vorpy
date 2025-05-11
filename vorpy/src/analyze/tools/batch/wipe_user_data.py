import tkinter as tk
from tkinter import filedialog
import os
import shutil

root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)
user_data = filedialog.askdirectory(title="Choose a user directory")

removal_directories = []
for roott, folders, files in os.walk(user_data):
    for folder in folders:
        for new_roott, sub_folders, filess in os.walk(roott + '/' + folder):
            if new_roott[-3:] == 'vor' or new_roott[-3:] == 'pow' or new_roott[-2:] == 'aw':
                removal_directories.append(new_roott)
            for file in filess:
                if 'vert' in file:
                    removal_directories.append(roott + '/' + folder + '/' + file)
                elif 'info.txt' in file:
                    removal_directories.append(roott + '/' + folder + '/' + file)

for dinkleschmitt in removal_directories:
    try:
        shutil.rmtree(dinkleschmitt + '/')
    except Exception:
        os.remove(dinkleschmitt)
    except FileNotFoundError:
        print(dinkleschmitt)
        pass
