import os
import sys

# Add the project root directory to the Python path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../"))
sys.path.append(project_root)

import tkinter as tk
from tkinter import filedialog
from Data.Analyze.tools.batch.get_files import get_files
from Data.Analyze.tools.CleanData.VerifyFiles import verify_pdb, verify_logs


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
        # Quick verification
        pdb_fl = pdb_fl if pdb_fl is not None and verify_pdb(pdb_fl, simple=True) else None
        aw_fl = aw_fl if aw_fl is not None and verify_logs(aw_fl, simple=True) else None
        pow_fl = pow_fl if pow_fl is not None and verify_logs(pow_fl, simple=True) else None
        if aw_fl is None:
            if os.path.exists(folder + '/' + subfolder + '/aw/aw_logs.pdb'):
                os.rename(folder + '/' + subfolder + '/aw/aw_logs.pdb', folder + '/' + subfolder + '/aw/aw_logs.csv')
                aw_fl = folder + '/' + subfolder + '/aw/aw_logs.csv'
        if pow_fl is None:
            if os.path.exists(folder + '/' + subfolder + '/pow/pow_logs.pdb'):
                os.rename(folder + '/' + subfolder + '/pow/pow_logs.pdb', folder + '/' + subfolder + '/pow/pow_logs.csv')
                pow_fl = folder + '/' + subfolder + '/pow/pow_logs.csv'

        # Get the subfolder information
        sub_info = subfolder.split('_')
        try:
            cv, den = float(sub_info[1]), float(sub_info[3])
        except ValueError:
            print(sub_info)
            continue
        except IndexError:
            print(sub_info)
            continue

        # Get the number for the file
        try:
            num = int(sub_info[-1])
        except ValueError:
            num = 0

        # If the number is too high
        if num > 19:
            while (cv, den, num) in extra_files:
                num += 1
            extra_files[(cv, den, num)] = {'aw': aw_fl is not None, 'pow': pow_fl is not None,
                                           'pdb': pdb_fl is not None, 'exists': True,
                                           'complete': not (aw_fl is None or pow_fl is None or pdb_fl is None),
                                           'subfolder': subfolder}
        # Check for repeats
        elif (cv, den, num) in checklist and pdb_fl is not None:
            while (cv, den, num) in checklist:
                num += 1
            if num > 19:
                while (cv, den, num) in checklist:
                    num += 1
                extra_files[(cv, den, num)] = {'aw': aw_fl is not None, 'pow': pow_fl is not None,
                                               'pdb': pdb_fl is not None, 'exists': True,
                                               'complete': not (aw_fl is None or pow_fl is None or pdb_fl is None),
                                               'subfolder': subfolder}
            else:
                while (cv, den, num) in checklist:
                    num += 1
                # Checklist
                checklist[(cv, den, num)] = {'aw': aw_fl is not None, 'pow': pow_fl is not None, 'pdb': pdb_fl is not None,
                                             'exists': True,
                                             'complete': not (aw_fl is None or pow_fl is None or pdb_fl is None),
                                             'subfolder': subfolder}
        else:
            while (cv, den, num) in checklist:
                num += 1
            # Checklist
            checklist[(cv, den, num)] = {'aw': aw_fl is not None, 'pow': pow_fl is not None, 'pdb': pdb_fl is not None,
                                         'exists': True,
                                         'complete': not (aw_fl is None or pow_fl is None or pdb_fl is None),
                                         'subfolder': subfolder}

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
                        print(checklist[(cv, den, i)])
                        if (cv, den) in vorpy_solves:
                            vorpy_solves[(cv, den)][i] = checklist[(cv, den, i)]
                        else:
                            vorpy_solves[(cv, den)] = {i: checklist[(cv, den, i)]}
                        foam_done += 1
                else:
                    print()
                    # Add to the foam solves
                    if (cv, den) in foam_makes:
                        foam_makes[(cv, den)].append(i)
                    else:
                        foam_makes[(cv, den)] = [i]
                    incomplete += 1

    # Print the missing foam numbers from the data
    print("Missing Foams:\n")
    for cv in cv_vals:
        print(cv, ": ", *[f"{den} - {len(foam_makes[(cv, den)]):02d} | " if (cv, den) in foam_makes else f"{den} - 00 | " for den in density_vals])

    # Print the missing foam numbers from the data
    print("\n\nMissing Foam Solves:\n")
    for cv in cv_vals:
        print(cv, ": ", *[f"{den} - {len(vorpy_solves[(cv, den)]):02d} | " if (cv, den) in vorpy_solves else f"{den} - 00 | " for den in density_vals])

    # Print the full data information
    print(f"\n\nNumber complete = {num_complete}/{total_count}\nFoam Complete = {num_complete + foam_done}/{total_count}\n"
          f"Number not made = {incomplete} / {total_count}\nNumber Extra = {len(extra_files)}")

    # Create a ledger of the missing foam solves
    # for

    # Return the information
    return foam_makes, vorpy_solves, extra_files

# def create_foam_scripts()


if __name__ == '__main__':
    os.chdir('../../../..')
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    my_folder = filedialog.askdirectory()
    print(my_folder)

    foams, vorpys, extras = completeness_check(cv_vals=[0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
                                               density_vals=[0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5],
                                               folder=my_folder)


