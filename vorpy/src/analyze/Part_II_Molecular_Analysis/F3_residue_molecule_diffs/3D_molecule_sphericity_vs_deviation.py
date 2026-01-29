import os
import sys
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning, module=r"PIL|matplotlib\.backends\._backend_tk")

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.plot_templates.bar import bar
from vorpy.src.analyze.tools.batch.get_files import get_all_files
from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2
from vorpy.src.calculations.calcs import calc_sphericity



def plot_data(plotting='Vol', diff='tot', exclude_keys=[], ylim=None):
    """Plots the data for the given plotting type and difference type"""
    # Get the files
    files = get_all_files()
    # Get the totals
    for key, value in files.items():
        if key in exclude_keys:
            continue
        # get the aw, pow, and prm volumes
        for log in ['aw', 'pow', 'prm']:
            # Read the logs
            logs = read_logs2(value[log], all_=False)
            vol = logs['group data']['Volume']
            sa = logs['group data']['Surface Area']
            files[key][log + ' vol'] = vol
            files[key][log + ' sa'] = sa
            files[key][log + ' sphericity'] = calc_sphericity(vol, sa)

    # Sort the files dictionary by key
    files = dict(sorted(files.items()))

    # Get the data and the labels


    # Get the