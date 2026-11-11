import os
from pathlib import Path


def default_output_root(sys):
    """Return the configured project output root, independent of cwd."""
    files = getattr(sys, 'files', None) or {}
    project_root = files.get('vpy_dir')
    if project_root is None:
        project_root = Path(__file__).resolve().parents[3]
    return Path(project_root).expanduser().resolve() / 'output'


def set_sys_dir(sys, dir_name=None, *, fail_if_exists=False):
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
    # Resolve defaults from the configured project root, not the process cwd.
    if dir_name is None:
        if sys.files['dir'] is not None:
            dir_name = str(Path(sys.files['dir']).expanduser() / sys.name)
        else:
            dir_name = str(default_output_root(sys) / sys.name)
    if fail_if_exists:
        try:
            os.makedirs(dir_name, exist_ok=False)
        except FileExistsError as error:
            raise FileExistsError(
                f'Refusing to reuse existing output directory: {os.path.abspath(dir_name)}'
            ) from error
        sys.files['dir'] = dir_name
        return

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
            os.makedirs(dir_name + i_str, exist_ok=False)
            break
        # If the file exists increment the counter and try creating the directory again
        except FileExistsError:
            i += 1
    # Set the output directory for the system
    sys.files['dir'] = dir_name + i_str
