from Data.Analyze.tools.batch.compile_logs import get_logs_and_pdbs
from System.system import System


pdb_files = get_logs_and_pdbs(False)

# Get the number of files
log_length = len(pdb_files)
# Loop through the loggys
for j, loggy in enumerate(pdb_files):
    # print the progress
    print("\rReading Log {}/{} - {}%".format(j+1, log_length, 100 * (j+1) / log_length), end="")
    # split the loggy value
    loggy_list = loggy.split('_')
    # Ge the cv and density
    cv, density = float(loggy_list[1]), float(loggy_list[3])
    # Get the pdb file address
    my_pdb = pdb_files[loggy]['pdb']
    # Create the simple system
    my_sys = System(my_pdb, simple=True)
    # Get the box dimensions
    box_dimensions = my_sys.data
    print(box_dimensions)
    # Loop through the balls in the system
    for i, ball in my_sys.balls.iterrows():
        # Get the box dimensions
        pass
        # Get the x, y, and z coordinates



