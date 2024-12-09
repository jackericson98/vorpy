import csv
import os

from Data.Analyze.tools.batch.compile_logs import get_logs_and_pdbs
from Data.Analyze.tools.compare.read_logs2 import read_logs2
import pandas as pd
from System.sys_funcs.calcs.calcs import calc_dist
from System.system import System

# First gather the log files
my_logs = get_logs_and_pdbs(False)
# Create the overlaps directory
olaps = {}
log_length = len(my_logs)
# Loop through each of the logs files looking for the overlap data
for j, loggy in enumerate(my_logs):
    # print the progress
    print("\rReading Log {}/{} - {}%".format(j+1, log_length, 100 * (j+1) / log_length), end="")
    # split the loggy value
    loggy_list = loggy.split('_')
    # Ge the cv and density
    cv, density = float(loggy_list[1]), float(loggy_list[3])
    # Get the files
    pdb_file, aw_logs_file, pow_logs_file = [my_logs[loggy][_] for _ in ['pdb', 'aw', 'pow']]
    # Get the logs dictionaries
    aw_logs = read_logs2(aw_logs_file, True, all_=False, balls=True)
    # System of balls
    my_sys = System(pdb_file, simple=True)
    # Create the dataframe for the logs
    aw_dataframe = pd.DataFrame(aw_logs['atoms'])
    # pow_logs = read_logs2(pow_logs_file, True, all_=False, balls=True)
    for i, ball in aw_dataframe.iterrows():
        # Loop through the balls
        overlaps = []
        # Loop through each of the neighbors looking for size
        for neighbor in ball['Neighbors']:
            # Get the neighbors rad and loc
            nrad, nloc = my_sys.balls['rad'][neighbor], my_sys.balls['loc'][neighbor]
            # Check the size
            if nrad > ball['rad']:
                # Add the overlap percentage
                overlaps.append(round(max(nrad + ball['rad'] - calc_dist(ball['loc'], nloc), 0) / ball['rad'], 5))
        # # Write the overlap data to match the overlaps file we already have
        with open('overlaps.csv', 'a') as my_overlap_file:
            # Create the csv_writer
            olap_csv = csv.writer(my_overlap_file)
            # Write the new line
            olap_csv.writerow([loggy, i] + overlaps)
