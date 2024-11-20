from Data.Analyze.tools.batch.compile_logs import get_logs_and_pdbs
from Data.Analyze.tools.compare.read_logs2 import read_logs2
import pandas as pd
from System.sys_funcs.calcs.calcs import calc_dist


my_logs = get_logs_and_pdbs(False)
# Create the overlaps directory
olaps = {}
# Loop through each of the logs files looking for the overlap data
for loggy in my_logs:
    # Get the
    # Ge the cv and density
    # Get the files
    pdb_file, aw_logs_file, pow_logs_file = [my_logs[loggy][_] for _ in ['pdb', 'aw', 'pow']]
    # Get the logs dictionaries
    aw_logs = read_logs2([aw_logs_file], True, all_=False, balls=True)
    # Create the dataframe for the logs
    aw_dataframe = pd.DataFrame(aw_logs[loggy]['atoms'])
    # pow_logs = read_logs2(pow_logs_file, True, all_=False, balls=True)
    # Loop through the balls
    overlaps = []
    for i, ball in aw_dataframe.iterrows():
        # Loop through each of the neighbors looking for size
        for neighbor in ball['Neighbors']:
            # Get the neighbors rad and loc
            nrad, nloc = aw_dataframe['rad'][neighbor], aw_dataframe['loc'][neighbor]
            # Check the size
            if nrad > ball['rad']:
                # Add the overlap percentage
                overlaps.append(max(nrad + ball['rad'] - calc_dist(ball['loc'], nloc), 0) / ball['rad'])

