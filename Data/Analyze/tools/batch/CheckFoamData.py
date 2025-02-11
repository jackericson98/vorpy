import os
import tkinter as tk
from tkinter import filedialog
from Data.Analyze.tools.batch.get_files import get_files


def completeness_check(cv_vals, density_vals, number_of_files=20, folder=None):
    # Get the folder to loop through
    if folder is None:
        folder = filedialog.askdirectory()

    # Get some information on the folder
    # Create a checklist for the subfolders
    checklist = {}
    extra_files = {}
    # Loop through the subfolders
    for subfolder in os.listdir(folder):
        # Get the full path for the subfolder
        full_path = os.path.join(folder, subfolder)

        # Get the files from the subfolder
        pdb_fl, aw_fl, pow_fl = get_files(full_path)
        # print(pdb_fl, aw_fl, pow_fl)

        # Get the subfolder information
        sub_info = subfolder.split('_')
        cv, den = float(sub_info[1]), float(sub_info[3])

        # Get the number for the file
        try:
            num = int(sub_info[-1])
        except ValueError:
            num = 0

        # If the number is too high
        if num > 19:
            extra_files[(cv, den, num)] = {'aw': aw_fl is not None, 'pow': pow_fl is not None,
                                           'pdb': pdb_fl is not None, 'exists': True,
                                           'complete': not (aw_fl is None or pow_fl is None or pdb_fl is None)}
        else:
            # Checklist
            checklist[(cv, den, num)] = {'aw': aw_fl is not None, 'pow': pow_fl is not None, 'pdb': pdb_fl is not None,
                                         'exists': True,
                                         'complete': not (aw_fl is None or pow_fl is None or pdb_fl is None)}

    total_count = len(cv_vals) * len(density_vals) * number_of_files
    num_complete, foam_done, incomplete = 0, 0, 0
    foam_makes = {}
    vorpy_solves = {}
    # Print the missing values
    for cv in cv_vals:
        for den in density_vals:
            for i in range(number_of_files):
                # Check if it is in the checklist
                if (cv, den, i) in checklist:
                    if checklist[(cv, den, i)]['complete']:
                        num_complete += 1
                    else:
                        if (cv, den) in vorpy_solves:
                            vorpy_solves[(cv, den)][i] = checklist[(cv, den, i)]
                        else:
                            vorpy_solves[(cv, den)] = {i: checklist[(cv, den, i)]}
                        foam_done += 1
                else:
                    # Add to the foam solves
                    if (cv, den) in foam_makes:
                        foam_makes[(cv, den)].append(i)
                    else:
                        foam_makes[(cv, den)] = [i]
                    incomplete += 1

    # Print the missing foam numbers from the data
    print("Missing Foams:\n")
    for cv in cv_vals:
        print(cv, ": ", *[f"{den} - {foam_makes[(cv, den)]} | " if (cv, den) in foam_makes else f"{den} - Complete" for den in density_vals])

    # Print the missing foam numbers from the data
    print("Missing Foam Solves:\n")
    for cv in cv_vals:
        print(cv, *[f"{den} - {vorpy_solves[(cv, den)].keys()} | " if (cv, den) in vorpy_solves else f"{den} - Complete" for den in density_vals])

    # Print the full data information
    print(f"Number complete = {num_complete}/{total_count}\nFoam Complete = {num_complete + foam_done}/{total_count}\n"
          f"Number not made = {incomplete} / {total_count}")

    # Return the information
    return foam_makes, vorpy_solves, extra_files

# def create_foam_scripts()


if __name__ == '__main__':
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    my_folder = filedialog.askdirectory()

    foams, vorpys, extras = completeness_check(cv_vals=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
                                               density_vals=[0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5],
                                               folder=my_folder)


