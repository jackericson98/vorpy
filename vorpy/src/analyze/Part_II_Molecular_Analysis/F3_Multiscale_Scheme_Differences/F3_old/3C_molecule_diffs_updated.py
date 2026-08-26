import os
import sys
import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning, module=r'PIL|matplotlib\.backends\._backend_tk')

vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '', '../..', '..', '..', '..'))
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.plot_templates.bar import bar
from vorpy.src.analyze.tools.batch.get_files import get_all_files
from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2

POWER_COLOR = '#d62728'
PRIMITIVE_COLOR = '#7f3fbf'


def _res_key(atom):
    return f"{atom['Residue']}_{atom['Residue Sequence']}"


def get_sol_neighbor_count(logs, sol_name='SOL'):
    """Count unique SOL residues sharing at least one surface with a non-SOL molecule atom."""
    atom_info = {int(atom['Index']): (str(atom['Residue']).strip(), _res_key(atom)) for _, atom in logs['atoms'].iterrows()}
    sol_neighbors = set()

    for _, surf in logs['surfs'].iterrows():
        b1, b2 = (int(_) for _ in surf['Balls'])
        a1, a2 = atom_info.get(b1), atom_info.get(b2)
        if a1 is None or a2 is None:
            continue
        r1, key1 = a1
        r2, key2 = a2
        if r1.upper() == sol_name.upper() and r2.upper() != sol_name.upper():
            sol_neighbors.add(key1)
        elif r2.upper() == sol_name.upper() and r1.upper() != sol_name.upper():
            sol_neighbors.add(key2)

    return len(sol_neighbors)


def _pct(value, reference):
    return 100.0 * (value - reference) / reference if reference != 0 else 0.0


def get_molecule_data(exclude_keys=None, sol_name='SOL'):
    """Collect whole-molecule volume, SA, and number of touching SOL residues for AW/POW/PRM."""
    if exclude_keys is None:
        exclude_keys = []
    files = get_all_files()
    data = {}

    for key, value in files.items():
        if key in exclude_keys:
            continue
        data[key] = {}
        for scheme in ['aw', 'pow', 'prm']:
            logs = read_logs2(value[scheme], all_=False, balls=True, surfs=True)
            data[key][scheme] = {
                'vol': float(logs['group data']['Volume']),
                'sa': float(logs['group data']['Surface Area']),
                'neighbors': get_sol_neighbor_count(logs, sol_name=sol_name)
            }

    return dict(sorted(data.items()))


def plot_data(data, ylim=None, absolute=True, plot_scheme='both'):
    x_names = list(data.keys())
    prefix = 'Abs %' if absolute else '%'
    titles = {
        'vol': 'Molecular Volume Difference',
        'sa': 'Molecular Surface Area Difference',
        'neighbors': 'Molecular Solvent Neighbor Count Difference'
    }
    for metric, ylabel in [('vol', f'{prefix} Volume Difference'), ('sa', f'{prefix} Surface Area Difference'),
                           ('neighbors', f'{prefix} SOL Neighbor Count Difference')]:
        pow_diff = [_pct(data[k]['pow'][metric], data[k]['aw'][metric]) for k in x_names]
        prm_diff = [_pct(data[k]['prm'][metric], data[k]['aw'][metric]) for k in x_names]
        title = titles[metric]
        if absolute:
            pow_diff = [abs(_) for _ in pow_diff]
            prm_diff = [abs(_) for _ in prm_diff]

        if absolute and ylim is None:
            max_y = max(pow_diff + prm_diff, default=0.0)
            metric_ylim = [0, max_y * 1.1 if max_y > 0 else 1.0]
        else:
            metric_ylim = ylim

        common = dict(x_names=x_names, Show=True, y_axis_title=ylabel, x_axis_title='Model', print_vals_on_bars=False,
                      legend_orientation='Vertical', y_range=metric_ylim, xlabel_size=30, ylabel_size=30, tick_width=2,
                      tick_length=12, xtick_label_size=25, ytick_label_size=25, x_tick_rotation=0, legend_loc='lower right')

        if plot_scheme == 'both':
            bar([pow_diff, prm_diff], colors=[POWER_COLOR, PRIMITIVE_COLOR], legend_names=['Pow vs AW', 'Prm vs AW'],
                title=title, **common)
        elif plot_scheme == 'pow':
            bar([pow_diff], colors=[POWER_COLOR], legend_names=['Pow vs AW'], title=title + 'pow', **common)
        elif plot_scheme == 'prm':
            bar([prm_diff], colors=[PRIMITIVE_COLOR], legend_names=['Prm vs AW'], title=title + 'prm', **common)
        elif plot_scheme == 'separate':
            bar([pow_diff], colors=[POWER_COLOR], legend_names=['Pow vs AW'], title=title + 'pow', **common)
            bar([prm_diff], colors=[PRIMITIVE_COLOR], legend_names=['Prm vs AW'], title=title + 'prm', **common)
        else:
            raise ValueError("plot_scheme must be 'both', 'pow', 'prm', or 'separate'")


if __name__ == '__main__':
    absolute = False
    plot_scheme = 'both'
    molecule_data = get_molecule_data(exclude_keys=['A', 'B', 'C'], sol_name='SOL')
    plot_data(molecule_data, ylim=None, absolute=absolute, plot_scheme=plot_scheme)
