import os
import sys
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning, module=r"PIL|matplotlib\.backends\._backend_tk")

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.plot_templates.bar import bar
from vorpy.src.analyze.tools.batch.get_files import get_all_files
from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


POWER_COLOR = '#d62728'
PRIMITIVE_COLOR = '#7f3fbf'


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

    # Sort the files dictionary by key
    files = dict(sorted(files.items()))

    # 
    bar(
        data=[[100 * (files[key]['pow vol'] - files[key]['aw vol']) / files[key]['aw vol'] for key in files if key not in exclude_keys],
              [100 * (files[key]['prm vol'] - files[key]['aw vol']) / files[key]['aw vol'] for key in files if key not in exclude_keys]],
        x_names=[key for key in files if key not in exclude_keys],
        Show=True,
        y_axis_title='% Difference',
        x_axis_title='Model',
        print_vals_on_bars=False,
        legend_orientation='Vertical',
        y_range=ylim,
        xlabel_size=30,
        ylabel_size=30,
        tick_width=2,
        tick_length=12,
        xtick_label_size=25,
        ytick_label_size=25,
        x_tick_rotation=0,
        colors=[POWER_COLOR, PRIMITIVE_COLOR],
        legend_names=["Pow vs AW", "Prm vs AW"],
        legend_loc='lower right'
    )

    bar(
        data=[[100 * (files[key]['pow sa'] - files[key]['aw sa']) / files[key]['aw sa'] for key in files if key not in exclude_keys],
              [100 * (files[key]['prm sa'] - files[key]['aw sa']) / files[key]['aw sa'] for key in files if key not in exclude_keys]],
        x_names=[key for key in files if key not in exclude_keys],
        Show=True,
        y_axis_title='% Difference',
        x_axis_title='Model',
        print_vals_on_bars=False,
        legend_orientation='Vertical',
        y_range=ylim,
        xlabel_size=30,
        ylabel_size=30,
        tick_width=2,
        tick_length=12,
        xtick_label_size=25,
        ytick_label_size=25,
        x_tick_rotation=0,
        colors=[POWER_COLOR, PRIMITIVE_COLOR],
        legend_names=["Pow vs AW", "Prm vs AW"],
        legend_loc='lower right'
    )


if __name__ == '__main__':
    plot_data(exclude_keys=['A', 'B', 'C'], ylim=[-4.5, 3.5])


