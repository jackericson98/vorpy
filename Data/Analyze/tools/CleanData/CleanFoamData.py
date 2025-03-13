import os
import sys

# Add the project root directory to the Python path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../"))
sys.path.append(project_root)

import os.path
import tkinter as tk
from tkinter import filedialog
import csv


def clean_foam_data():

    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)

    files = []

    while True:
        new_file = filedialog.askopenfilename(title='Get Another Foam File? >>>   ')
        print(new_file)

        if new_file == '':
            break
        files.append(new_file)
        my_answer = input("Grab another file?")
        if my_answer.lower().strip() not in {'', 'y', 'yes'}:
            break

    # File structure: sorted by cv density
    info = {}
    for file in files:
        with open(file, 'r') as read_file:
            read_csv = csv.reader(read_file)
            for line in read_csv:
                print(line)
                try:
                    file_name, cv, density = line[0], line[3], line[5]
                except IndexError:
                    continue
                if file_name in info:
                    continue
                else:
                    info[file_name] = line
    with open(os.path.dirname(files[0]) + '/all_foam_data.csv', 'w') as foam_data:
        my_writer = csv.writer(foam_data)
        for _ in info:
            line = [entry.strip("\n") for entry in info[_]]
            line[0] = line[0].split('/')[-1]
            my_writer.writerow(line)



if __name__ == '__main__':
    clean_foam_data()
