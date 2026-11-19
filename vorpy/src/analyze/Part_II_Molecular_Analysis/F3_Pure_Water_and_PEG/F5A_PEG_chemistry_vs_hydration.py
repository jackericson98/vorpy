"""Plot matched PEG atom volume and surface-area changes against AW contacts.

Both signed percentages are 100 * (Power - AW) / AW, using atomic volume and
surface area respectively. Observations are pooled atom/frame pairs, not
independent molecular replicates. F5 readers check PDB identity, surface
partitions, and complete, valid cells; no geometry is calculated here.
"""
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

if not __package__:
    import sys
    root = next(parent for parent in Path(__file__).resolve().parents
                if (parent / 'pyproject.toml').is_file())
    sys.path.insert(0, str(root))
    __package__ = ('vorpy.src.analyze.Part_II_Molecular_Analysis.F3_Pure_Water_and_PEG'
                   if 'F3_Pure_Water_and_PEG' in Path(__file__).parts else
                   'vorpy.src.analyze.Part_II_Molecular_Analysis.F3_SimpleModels')

from .F5A_exposed_vs_buried import assign_repeat_types
from vorpy.src.analyze.Part_II_Molecular_Analysis.Shared_Utilities.simple_models.F5_common import (
    load_scheme, match_schemes, read_pdb, save_figure, setup_plots,
)
from vorpy.src.analyze.Part_II_Molecular_Analysis.Shared_Utilities.simple_models.F5_config import (
    ATOM_KEY, LEGACY_OUTPUT_ROOT, Config,
)
try:
    from .support.F5_peg_dataset import screened_runs, paths, layer1_atoms
except ModuleNotFoundError:
    from .F5_peg_dataset import screened_runs, paths, layer1_atoms


# The AW/Power outputs are polymer scoped; AW surface endpoints are mapped to
# water molecules in each full-frame PDB to calculate distinct-water contacts.
INPUT_ROOT = Path(
    r'C:\Users\Ryzen-3600\Versiti Dropbox\John Ericson\Versiti lab\Jack'
    r'\Ericson Voronoi DNA\P2\NewModels\PEG\Reruns\layers_4'
)
OUTPUT_DIR = LEGACY_OUTPUT_ROOT / 'PEG' / 'F5A'
MODELS = ('EG', 'PEG300', 'PEG1000', 'PEG2000', 'PEG3000', 'PEG4000', 'PEG6000')
FIGURE_SIZE = (6.5, 4.5)
DPI = 200
FONT_SIZE = 12
CPK_COLORS = {'C': '#333333', 'H': '#f2f2f2', 'O': '#ff0000'}
POINT_SIZE = 10
POINT_ALPHA = 0.5
TERMINAL_POINT_SIZE = 23
TERMINAL_ALPHA = 0.75
MEAN_SIZE = 67.5
TERMINAL_MEAN_SIZE = 90
MEAN_MIN_N = 10
ZERO_LINE_COLOR = '#555555'
Y_LIMITS = (-25, 25)
TYPE_STYLES = {
    'Carbon': {'element': 'C', 'marker': 'o'},
    'Hydrogen': {'element': 'H', 'marker': 'o'},
    'Ether oxygen': {'element': 'O', 'marker': 'o'},
    'Terminal oxygen': {'element': 'O', 'marker': 'D'},
}


def load_frame(model, frame, config):
    """Load/match a frame, retaining AW distinct-water coordination per atom."""
    selected = screened_runs()
    selected = selected[selected.molecule.eq(model) & selected.frame.eq(frame)]
    if len(selected) != 3:
        raise ValueError(f'Unscreened PEG frame: {model}/{frame}')
    aw_row = next(selected[selected.scheme.eq('aw')].itertuples(index=False))
    power_row = next(selected[selected.scheme.eq('pow')].itertuples(index=False))
    aw_log, pdb_path = paths(aw_row)
    power_log, _ = paths(power_row)
    required = (aw_log, power_log, pdb_path)
    if not all(path.is_file() for path in required):
        raise ValueError(f'Screened AW/Power log or frame PDB missing: {model}/{frame}')

    pdb = read_pdb(pdb_path, 'PEG', config)
    metadata = {'polymer': 'PEG', 'model': model, 'frame': frame}
    aw = load_scheme(aw_log, pdb, config, dict(metadata, scheme='aw'))['atoms']
    power = load_scheme(power_log, pdb, config, dict(metadata, scheme='pow'))['atoms']
    invalid = {
        'aw': int((~aw.primary_valid.fillna(False)).sum()),
        'pow': int((~power.primary_valid.fillna(False)).sum()),
    }
    invalid['layer1_waters'] = len(layer1_atoms(aw_row, pdb)) // 3

    # The established F5 PEG naming and covalent-distance checks distinguish
    # terminal oxygens from ether oxygens and exclude O-H hydrogens.
    atoms = pd.concat([aw, power], ignore_index=True)
    aw_types = assign_repeat_types(atoms[atoms.scheme.eq('aw')])
    atoms = atoms.drop(columns='atom_class').merge(
        aw_types[ATOM_KEY + ['atom_class']], on=ATOM_KEY,
        how='left', validate='many_to_one',
    )
    # This shared matcher checks atom index, residue/name/coordinates, and cell
    # validity in both schemes while carrying AW water_coordination forward.
    pairs, _ = match_schemes(atoms, config)
    return pairs, '', invalid


def plot_data(pairs, metric):
    """Calculate one signed percentage difference and count its exclusions."""
    classes = pairs.atom_class_aw.isin(TYPE_STYLES)
    aw_value = pd.to_numeric(pairs[f'{metric}_aw'], errors='coerce')
    power_value = pd.to_numeric(pairs[f'{metric}_s'], errors='coerce')
    contacts = pd.to_numeric(pairs.water_coordination_aw, errors='coerce')
    denominator_valid = np.isfinite(aw_value) & aw_value.gt(0)
    power_valid = np.isfinite(power_value)
    contacts_valid = np.isfinite(contacts) & contacts.ge(0) & contacts.eq(contacts.round())
    valid = pairs.pair_valid & classes & denominator_valid & power_valid & contacts_valid
    points = pairs.loc[valid, ATOM_KEY + ['atom_class_aw']].copy()
    points['water_contacts'] = contacts.loc[valid].astype(int)
    points[f'delta_{metric}_pct'] = 100 * (power_value.loc[valid] - aw_value.loc[valid]) / aw_value.loc[valid]
    excluded = {
        'invalid_or_unmatched_scheme_pairs': int((~pairs.pair_valid).sum()),
        'valid_pairs_with_zero_or_invalid_AW_denominator': int((pairs.pair_valid & ~denominator_valid).sum()),
        'valid_pairs_with_nonfinite_Power_value': int((pairs.pair_valid & denominator_valid & ~power_valid).sum()),
        'valid_measurements_with_invalid_AW_contacts': int((pairs.pair_valid & denominator_valid & power_valid & ~contacts_valid).sum()),
        'valid_measurements_outside_four_atom_classes': int((pairs.pair_valid & denominator_valid & power_valid & contacts_valid & ~classes).sum()),
    }
    return points, excluded


def main():
    config = Config()
    observations = {'volume': [], 'surface_area': []}
    exclusions = {metric: {} for metric in observations}
    invalid_cells = {'aw': 0, 'pow': 0}
    runs = screened_runs()
    coverage = []

    for model in MODELS:
        model_frames = sorted(runs.loc[runs.molecule.eq(model), 'frame'].unique())
        model_counts = {metric: 0 for metric in observations}
        skipped = 0
        waters = 0
        for frame in model_frames:
            # A screened input failure is an error, not a silently smaller cohort.
            pair, reason, invalid = load_frame(model, frame, config)
            if pair is None:
                raise ValueError(f'Screened frame unavailable: {model}/{frame}: {reason}')
            for scheme in invalid_cells:
                invalid_cells[scheme] += invalid[scheme]
            waters += invalid['layer1_waters']
            coverage.append(dict(model=model, frame=frame, matched_atoms=int(pair.pair_valid.sum()),
                                 layer1_water_molecules=invalid['layer1_waters']))
            for metric in observations:
                points, frame_exclusions = plot_data(pair, metric)
                observations[metric].append(points)
                model_counts[metric] += len(points)
                for key, value in frame_exclusions.items():
                    exclusions[metric][key] = exclusions[metric].get(key, 0) + value

        print(f'{model}: {len(model_frames)} frames, {waters} layer-1 water molecules; '
              f'volume {model_counts["volume"]}, surface area {model_counts["surface_area"]} '
              'included atom/frame observations; '
              f'{skipped} frames skipped', flush=True)

    if not any(observations.values()):
        raise ValueError('No matched PEG AW/Power observations were available.')
    print(f'Excluded incomplete/invalid cells: AW {invalid_cells["aw"]}; '
          f'Power {invalid_cells["pow"]}', flush=True)
    summaries = {}
    metric_labels = {'volume': 'volume', 'surface_area': 'surface area'}
    for metric in observations:
        if not observations[metric]:
            continue
        points = pd.concat(observations[metric], ignore_index=True)
        print(f'Included {metric_labels[metric]} observations by atom type:', flush=True)
        for atom_type in TYPE_STYLES:
            print(f'  {atom_type}: {int(points.atom_class_aw.eq(atom_type).sum())}', flush=True)
        for key, count in exclusions[metric].items():
            print(f'  Excluded {key.replace("_", " ")}: {count}', flush=True)

        value_column = f'delta_{metric}_pct'
        summary = (points.groupby(['atom_class_aw', 'water_contacts'], observed=True)
            [value_column].agg(['mean', 'std', 'count']).reset_index()
            .rename(columns={'atom_class_aw': 'atom_type', 'mean': f'mean_delta_{metric}_pct',
                             'std': f'std_delta_{metric}_pct', 'count': f'n_{metric}'}))
        if metric == 'volume':
            summary = summary.rename(columns={'n_volume': 'n'})
        summaries[metric] = summary
        plot_metric(points, metric)

    # Preserve the original volume summary columns and append surface-area stats.
    summary = summaries['volume'].merge(
        summaries['surface_area'], on=['atom_type', 'water_contacts'], how='outer',
    )
    summary = summary[['atom_type', 'water_contacts', 'mean_delta_volume_pct',
                       'std_delta_volume_pct', 'n', 'mean_delta_surface_area_pct',
                       'std_delta_surface_area_pct', 'n_surface_area']]
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUTPUT_DIR / 'F5A_PEG_chemistry_vs_hydration.csv'
    summary.to_csv(csv_path, index=False)
    pd.DataFrame(coverage).to_csv(OUTPUT_DIR / 'F5A_coverage.csv', index=False)
    print(f'Saved {csv_path}', flush=True)


def plot_metric(points, metric):
    """Match the original F5A scatter style and save one measure-specific plot."""
    config = Config(figure_size=FIGURE_SIZE, dpi=DPI, font_size=FONT_SIZE,
                    export_formats=('png', 'svg'))
    setup_plots(config)
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    ax.tick_params(labelsize=config.font_size)
    y_column = f'delta_{metric}_pct'
    for atom_type in ['Hydrogen', 'Carbon', 'Ether oxygen', 'Terminal oxygen']:
        atom_points = points.loc[points.atom_class_aw.eq(atom_type)]
        if atom_points.empty:
            continue
        style = TYPE_STYLES[atom_type]
        color = CPK_COLORS[style['element']]
        terminal = atom_type == 'Terminal oxygen'
        ax.scatter(atom_points.water_contacts, atom_points[y_column],
                   s=TERMINAL_POINT_SIZE if terminal else POINT_SIZE,
                   alpha=TERMINAL_ALPHA if terminal else POINT_ALPHA,
                   linewidths=.3, edgecolors='#666666', color=color,
                   marker=style['marker'], label=atom_type)
        contact_stats = atom_points.groupby('water_contacts')[y_column].agg(['mean', 'count'])
        supported_means = contact_stats.loc[contact_stats['count'].ge(MEAN_MIN_N), 'mean']
        if not supported_means.empty:
            ax.scatter(supported_means.index, supported_means,
                       s=TERMINAL_MEAN_SIZE if terminal else MEAN_SIZE,
                       linewidths=1.4,
                       edgecolors='white' if style['element'] == 'C' else '#222222',
                       color=color, marker=style['marker'], zorder=4)

    ax.set_xlabel('Solvent contacts in AW (water molecules)')
    ax.set_ylabel(f'Signed {"volume" if metric == "volume" else "surface area"} difference from AW (%)')
    ax.set_ylim(Y_LIMITS)
    ax.set_title('PEG atoms: Power versus AW\nLarge markers show means for n ≥ 10')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.axhline(0, color=ZERO_LINE_COLOR, linewidth=.8)
    legend = ax.legend(title='PEG atom type', frameon=False, markerscale=2,
                       loc='upper left', bbox_to_anchor=(1.01, 1))
    for handle in legend.legend_handles:
        handle.set_alpha(1)
    stem = f'F5A_PEG_{metric}_vs_water_contacts'
    save_figure(fig, OUTPUT_DIR / stem, config)
    if metric == 'volume':
        # Keep the previous volume-figure filename from becoming a stale result.
        for extension in config.export_formats:
            legacy = OUTPUT_DIR / f'F5A_PEG_chemistry_vs_hydration.{extension}'
            if legacy.is_file():
                shutil.copy2(OUTPUT_DIR / f'{stem}.{extension}', legacy)
    print(f'Saved PNG/SVG {OUTPUT_DIR / stem}', flush=True)


if __name__ == '__main__':
    main()
