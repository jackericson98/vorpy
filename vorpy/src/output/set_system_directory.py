import os


def set_sys_dir(sys, dir_name=None):
    """
    Creates and sets up an output directory for system data with automatic naming and collision handling.
    
    This function manages the creation of output directories for system data, handling several cases:
    - Creates a user_data directory if it doesn't exist
    - Generates a unique directory name by appending incrementing numbers if needed
    - Supports both custom directory paths and default locations
    
    :param sys: System object that needs an output directory
    :param dir_name: Optional custom directory name/path. If None, generates a default path
    :return: None
    """
    # Check to see if the system has a name
    if sys.name is None:
        return
    # Make sure a user_data path exists
    if sys.files['vpy_dir'] is not None and not os.path.exists(sys.files['vpy_dir'] + "/output"):
        os.mkdir(sys.files['vpy_dir'] + "/output")

    # If no outer directory was specified use the directory outside the current one
    if dir_name is None:
        if sys.files['dir'] is not None:
            dir_name = sys.files['dir'] + '/' + sys.name
        elif sys.files['vpy_dir'] is not None:
            dir_name = sys.files['vpy_dir'] + "/output/" + sys.name
        else:
            dir_name = os.getcwd() + "/output/" + sys.name
    # Catch for existing directories. Keep trying out directories until one doesn't exist
    i = 0
    while True:
        # Try creating the directory with the system name + the current i_string
        try:
            # Create a string variable for the incrementing variable
            i_str = '_' + str(i)
            # If no file with the system name exists change the string to empty
            if i == 0:
                i_str = ""
            # Try to create the directory
            os.mkdir(dir_name + i_str)
            break
        # If the file exists increment the counter and try creating the directory again
        except FileExistsError:
            i += 1
    # Set the output directory for the system
    sys.files['dir'] = dir_name + i_str
