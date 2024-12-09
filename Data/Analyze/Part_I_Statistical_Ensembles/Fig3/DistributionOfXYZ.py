from Data.Analyze.tools.batch.compile_logs import get_logs_and_pdbs
from System.system import System
import csv


pdb_files = get_logs_and_pdbs(False)

# Get the number of files
log_length = len(pdb_files)
# Loop through the loggys
for j, loggy in enumerate(pdb_files):
    # print the progress
    print("\rReading Log {}/{} - {}%".format(j+1, log_length, 100 * (j+1) / log_length), end="")
    try:
        # Get the pdb file address
        my_pdb = pdb_files[loggy]['pdb']
    except KeyError:
        continue
    # Create the simple system
    try:
        my_sys = System(my_pdb, simple=True)
    except ValueError:
        continue
    # split the loggy value
    loggy_list = loggy.split('_')
    # Ge the cv and density
    cv, density = float(loggy_list[1]), float(loggy_list[3])
    # Check to see if the cv and density have been used as a key yet
    # Get the box dimensions
    box_dimensions = float(my_sys.data[0][2])
    my_data = []
    # Loop through the balls in the system
    for i, ball in my_sys.balls.iterrows():
        if i > 999:
            continue
        # Get the adjusted x, y, and z coordinates
        x, y, z = [_ / box_dimensions for _ in ball['loc']]
        # Add the data to the dictionary
        my_data.append((x, y, z))

    # add the data to the file
    with open('location_data.csv', 'a') as locy_data:
        # Create the csv writer
        loc_writer = csv.writer(locy_data)
        # Write the line
        loc_writer.writerow([cv, density] + my_data)




