import os
from os import path


def get_files(folder):
    """Traverse the folder and extract required PDB and log files.
    Logs are identified based on subfolder names containing 'aw' or 'pow'.
    """
    pdb_file, aw_logs, pow_logs = None, None, None

    for root_dir, sub_dirs, files in os.walk(folder):
        for file in files:
            if file.endswith('.pdb') and 'atoms' not in file and 'diff' not in file:
                pdb_file = path.join(root_dir, file)

        # Check sub-subfolders for logs
        for sub_dir in sub_dirs:
            if 'aw' in sub_dir.lower():
                for file in os.listdir(path.join(root_dir, sub_dir)):
                    if file.endswith('.csv'):
                        aw_logs = path.join(root_dir, sub_dir, file)
            elif 'pow' in sub_dir.lower():
                for file in os.listdir(path.join(root_dir, sub_dir)):
                    if file.endswith('.csv'):
                        pow_logs = path.join(root_dir, sub_dir, file)
    return pdb_file, aw_logs, pow_logs
