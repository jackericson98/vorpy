"""Compare matched PEG Power and AW cells using the cached F5 atom table.

Run from the repository root:
python -m vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5A_exposed_vs_buried
Direct-file execution, including PyCharm's Run button, is also supported.

Inputs and calculations, in the order used by main():
1. Read atoms.csv, classify atoms using the AW structure, and match AW to Power
   by polymer, model, frame and atom index. The shared F5 matcher also checks
   atom identity, coordinates and cell validity.
2. Compute 100 * (Power - AW) / AW for volume, surface area and logged neighbor
   count. Absolute plots take the magnitude of each atom's change before
   averaging. Missing measurements and zero AW denominators are excluded.
3. Group bars by atom type and AW water surface-area fraction. Exposed means
   strictly above EXPOSURE_THRESHOLD. Error bars show sample SD (ddof=1).
4. Use only exposed atoms in the surface-fraction scatter. The water-contact
   scatter instead uses cached AW water_coordination: one count per distinct
   water residue with a qualifying Voronoi face, including zero-contact atoms.
   Its large markers average the points at each atom type and contact count
   with at least 10 observations.

Each observation is one atom in one frame. Pooled means give larger models
more weight because they contribute more atoms. Terminal O-H hydrogens are
listed in the atom-type audit CSV and excluded from the plotted types.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

if not __package__:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[5]))
    __package__ = 'vorpy.src.analyze.Part_II_Molecular_Analysis.F5'

from .F5_config import Config, ATOM_KEY, KEY, METRICS as F5_METRICS
from .F5_common import match_schemes, setup_plots, save_figure, write_csv


# Configuration: exposure is strictly greater than this AW atomic water fraction.
EXPOSURE_THRESHOLD = 0.01
INPUT_TABLE = Path(__file__).parent / 'output/summaries/atoms.csv'
OUTPUT_DIR = Path(__file__).parent / 'output/F5A'
FIGURE_SIZE = (6.5, 4.5)
DPI = 200
CPK_COLORS = {'C': '#333333', 'H': '#f2f2f2', 'O': '#ff0000'}
TYPE_STYLES = {
    'Carbon': {'element': 'C', 'marker': 'o', 'hatch': ''},
    'Hydrogen': {'element': 'H', 'marker': 'o', 'hatch': ''},
    'Ether oxygen': {'element': 'O', 'marker': 'o', 'hatch': ''},
    'Terminal oxygen': {'element': 'O', 'marker': 'D', 'hatch': '///'},
}
METRICS = {'volume': 'Atomic volume', 'surface_area': 'Atomic surface area',
           'contacts': 'Number of contacts (Number of Neighbors)'}
CATEGORIES = ['Exposed', 'Not Exposed']
COORDINATE_COLUMNS = ['x', 'y', 'z']


def assign_repeat_types(atoms):
    """Classify each AW PEG chain and record reasons for unassigned atoms."""
    roles = atoms[ATOM_KEY].copy()
    roles['atom_class'] = None
    roles['attached_heavy_atom_index'] = np.nan
    roles['type_exclusion'] = 'Unvalidated PEG backbone'
    chain_keys = KEY + ['residue_id']
    for _, chain in atoms.groupby(chain_keys, sort=False, dropna=False):
        if not chain.order_valid.eq(True).all():
            continue

        is_hydrogen = chain.element.isin(['H', 'D'])
        heavy_atoms = chain.loc[~is_hydrogen].sort_values('atom_index')
        repeat_count = int(chain.repeat_count.iloc[0])
        expected_elements = ['O'] + ['C', 'C', 'O'] * repeat_count
        expected_names = [
            f'{element}{position}' for position, element in enumerate(expected_elements, start=1)
        ]
        backbone_matches = (
            repeat_count >= 1
            and heavy_atoms.element.tolist() == expected_elements
            and heavy_atoms.atom_name.tolist() == expected_names
        )
        if not backbone_matches:
            continue

        heavy_types = heavy_atoms.element.map({'C': 'Carbon', 'O': 'Ether oxygen'}).to_numpy(dtype=object)
        heavy_types[[0, -1]] = 'Terminal oxygen'
        roles.loc[heavy_atoms.index, 'atom_class'] = heavy_types
        roles.loc[heavy_atoms.index, 'type_exclusion'] = ''

        hydrogen_atoms = chain.loc[is_hydrogen]
        if hydrogen_atoms.empty:
            continue

        hydrogen_xyz = hydrogen_atoms[COORDINATE_COLUMNS].to_numpy()
        heavy_xyz = heavy_atoms[COORDINATE_COLUMNS].to_numpy()
        # Rows are hydrogens, columns are heavy atoms in the same chain. Distances are in Å.
        offsets = hydrogen_xyz[:, None, :] - heavy_xyz[None, :, :]
        distances = np.linalg.norm(offsets, axis=2)
        nearest_heavy_positions = distances.argmin(axis=1)
        nearest_distances = distances.min(axis=1)
        second_nearest_distances = np.partition(distances, 1, axis=1)[:, 1]
        # Require a plausible covalent distance and a 0.05 Å margin over the next candidate.
        valid_attachment = (
            (nearest_distances >= .6)
            & (nearest_distances <= 1.35)
            & (second_nearest_distances - nearest_distances > .05)
        )
        for row, hydrogen_index in enumerate(hydrogen_atoms.index):
            if not valid_attachment[row]:
                roles.loc[hydrogen_index, 'type_exclusion'] = 'Unresolved hydrogen attachment'
                continue
            parent_position = nearest_heavy_positions[row]
            roles.loc[hydrogen_index, 'attached_heavy_atom_index'] = heavy_atoms.atom_index.iloc[parent_position]
            if heavy_types[parent_position] != 'Carbon':
                roles.loc[hydrogen_index, 'type_exclusion'] = 'O-H hydrogen'
                continue
            roles.loc[hydrogen_index, 'atom_class'] = 'Hydrogen'
            roles.loc[hydrogen_index, 'type_exclusion'] = ''
    return roles


def load_pairs(input_table=INPUT_TABLE, config=None):
    """Read the cached atoms, classify from AW and attach the raw neighbor counts."""
    if not input_table.is_file():
        raise FileNotFoundError(f'{input_table} is missing; run F5_run_all first.')
    columns = list(dict.fromkeys(ATOM_KEY + F5_METRICS + [
        'scheme', 'primary_valid', 'neighbors', 'atom_name', 'residue_name',
        'residue_sequence', 'chain', 'x', 'y', 'z', 'Number of Neighbors',
        'element', 'residue_id', 'order_valid', 'repeat_count']))
    atoms = pd.read_csv(input_table, usecols=columns)
    atoms = atoms[atoms.scheme.isin(['aw', 'pow']) & atoms.polymer.eq('PEG')].copy()
    if not {'aw', 'pow'}.issubset(set(atoms.scheme)):
        raise ValueError('The F5 atom table must contain both AW and Power observations.')
    atoms['neighbors'] = atoms.neighbors.map(lambda value: frozenset(json.loads(value)))
    aw_atoms = atoms.loc[atoms.scheme.eq('aw')]
    atom_types = assign_repeat_types(aw_atoms)
    atoms = atoms.merge(atom_types[ATOM_KEY + ['atom_class']], on=ATOM_KEY, how='left', validate='many_to_one')
    # The shared matcher checks atom/residue identity, coordinates, and F5 validity.
    pairs, _ = match_schemes(atoms, config or Config())
    # Bars for contact count compare total logged neighbors. The water scatter
    # uses water_coordination_aw from the shared matcher instead.
    neighbor_counts = atoms.pivot(index=ATOM_KEY, columns='scheme', values='Number of Neighbors')
    neighbor_counts = neighbor_counts.rename(columns={'aw': 'contacts_aw', 'pow': 'contacts_s'})
    pairs = pairs.merge(neighbor_counts.reset_index(), on=ATOM_KEY, how='left', validate='one_to_one')
    audit_columns = ATOM_KEY + ['attached_heavy_atom_index', 'type_exclusion']
    return pairs.merge(atom_types[audit_columns], on=ATOM_KEY, how='left', validate='one_to_one')


def percentage_differences(aw_values, power_values, *, signed):
    """Calculate per-atom changes for finite values with nonzero AW denominators."""
    signed_change = 100 * (power_values - aw_values) / aw_values
    return signed_change if signed else signed_change.abs()


def summarize(pairs, threshold=EXPOSURE_THRESHOLD, *, signed=False):
    """Return pooled means, sample SDs and sample counts for each bar."""
    if not 0 <= threshold <= 1:
        raise ValueError('EXPOSURE_THRESHOLD must lie between zero and one.')
    aw_water_fraction = pd.to_numeric(pairs.water_fraction_aw, errors='coerce')
    valid_exposure = np.isfinite(aw_water_fraction) & aw_water_fraction.between(0, 1)
    classified_atom = pairs.atom_class_aw.isin(TYPE_STYLES)
    eligible = pairs.loc[pairs.pair_valid & valid_exposure & classified_atom].copy()
    exposed = aw_water_fraction.loc[eligible.index].gt(threshold)
    eligible['exposure'] = np.where(exposed, 'Exposed', 'Not Exposed')
    exclusions = {'invalid_or_unmatched_pairs': int((~pairs.pair_valid).sum()),
                  'invalid_aw_exposure': int((pairs.pair_valid & ~valid_exposure).sum()),
                  'outside_plot_types': int((pairs.pair_valid & valid_exposure & ~classified_atom).sum())}
    difference = 'signed' if signed else 'absolute'
    rows = []
    for metric, label in METRICS.items():
        for category in CATEGORIES:
            for atom_type in TYPE_STYLES:
                group = eligible[eligible.atom_class_aw.eq(atom_type) & eligible.exposure.eq(category)]
                aw_values = pd.to_numeric(group[f'{metric}_aw'], errors='coerce')
                power_values = pd.to_numeric(group[f'{metric}_s'], errors='coerce')
                finite = np.isfinite(aw_values) & np.isfinite(power_values)
                zero_aw = aw_values.eq(0)
                valid = finite & ~zero_aw
                differences = percentage_differences(aw_values[valid], power_values[valid], signed=signed)
                rows.append({
                    'metric': label,
                    'metric_key': metric,
                    'exposure': category,
                    'polymer': 'PEG',
                    'atom_type': atom_type,
                    f'mean_{difference}_percentage_difference': differences.mean(),
                    f'sd_{difference}_percentage_difference': differences.std(ddof=1),
                    'n_valid': int(valid.sum()),
                    'n_candidates': len(group),
                    'n_aw_zero': int(zero_aw.sum()),
                    'n_nonfinite': int((~finite & ~zero_aw).sum()),
                    'n_models': group.loc[valid, 'model'].nunique(),
                    'n_model_frames': len(group.loc[valid, KEY].drop_duplicates()),
                    'threshold': threshold,
                })
    return pd.DataFrame(rows), exclusions


def plot(summary, *, signed=False, output_dir=OUTPUT_DIR):
    """Plot the means and sample SDs already calculated by summarize()."""
    import matplotlib.pyplot as plt
    config = Config(figure_size=FIGURE_SIZE, dpi=DPI)
    setup_plots(config)
    difference = 'signed' if signed else 'absolute'
    suffix = '_signed' if signed else ''
    for metric, label in METRICS.items():
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        metric_summary = summary.loc[summary.metric_key.eq(metric)]
        for position, (atom_type, style) in enumerate(TYPE_STYLES.items()):
            type_summary = metric_summary.loc[metric_summary.atom_type.eq(atom_type)]
            by_exposure = type_summary.set_index('exposure').loc[CATEGORIES]
            bar_offset = (position - (len(TYPE_STYLES) - 1) / 2) * .18
            x_positions = np.arange(len(CATEGORIES)) + bar_offset
            means = by_exposure[f'mean_{difference}_percentage_difference'].to_numpy(dtype=float)
            sample_sd = by_exposure[f'sd_{difference}_percentage_difference'].to_numpy(dtype=float)
            ax.bar(x_positions, means, width=.175, color=CPK_COLORS[style['element']],
                   edgecolor='#333333', linewidth=.6, hatch=style['hatch'], label=atom_type)
            has_sd = np.isfinite(sample_sd)
            ax.errorbar(x_positions[has_sd], means[has_sd], yerr=sample_sd[has_sd], fmt='none',
                        ecolor='#929292' if style['element'] == 'C' else '#222222',
                        elinewidth=1, capsize=3)
            for x_position, value in zip(x_positions, means):
                if not np.isfinite(value):
                    ax.text(x_position, .03, 'n/a', ha='center', transform=ax.get_xaxis_transform())
        ax.set_xticks(range(len(CATEGORIES)), CATEGORIES)
        ax.set_xlabel(f'AW atomic water SA fraction (exposed > {100*EXPOSURE_THRESHOLD:g}%)')
        ax.set_ylabel(f'Mean {difference} percentage difference from AW (%)')
        if signed:
            ax.axhline(0, color='#555555', linewidth=.8)
        ax.set_title(f'PEG {label.lower()}: Power versus AW\nMean ± sample SD')
        handles, labels = ax.get_legend_handles_labels()
        fig.legend(handles, labels, loc='lower center', ncol=len(TYPE_STYLES), frameon=False)
        save_figure(fig, output_dir / f'F5A_{metric}{suffix}', config)


def volume_points(pairs, threshold=EXPOSURE_THRESHOLD, *, signed=False, by_water_contacts=False):
    """Build the exact point table used for each scatter plot and CSV."""
    aw_volume = pd.to_numeric(pairs.volume_aw, errors='coerce')
    power_volume = pd.to_numeric(pairs.volume_s, errors='coerce')
    valid_volume = np.isfinite(aw_volume) & aw_volume.ne(0) & np.isfinite(power_volume)
    classified_atom = pairs.atom_class_aw.isin(TYPE_STYLES)
    valid = pairs.pair_valid & valid_volume & classified_atom
    if by_water_contacts:
        # Counts are integers. Zero is a valid observation for this plot.
        aw_water_contacts = pd.to_numeric(pairs.water_coordination_aw, errors='coerce')
        valid &= (
            np.isfinite(aw_water_contacts)
            & aw_water_contacts.ge(0)
            & aw_water_contacts.eq(aw_water_contacts.round())
        )
    else:
        aw_water_fraction = pd.to_numeric(pairs.water_fraction_aw, errors='coerce')
        valid &= (
            np.isfinite(aw_water_fraction)
            & aw_water_fraction.between(0, 1)
            & aw_water_fraction.gt(threshold)
        )
    points = pairs.loc[valid, ATOM_KEY + ['atom_class_aw', 'volume_aw', 'volume_s']].copy()
    if by_water_contacts:
        points['aw_water_contacts'] = aw_water_contacts[valid].astype(int)
    else:
        points['solvent_exposure_percent'] = 100 * aw_water_fraction[valid]
    difference = 'signed' if signed else 'absolute'
    points[f'{difference}_volume_difference_percent'] = percentage_differences(
        aw_volume[valid], power_volume[valid], signed=signed,
    )
    return points


def plot_volume_scatter(points, *, signed=False, by_water_contacts=False, output_dir=OUTPUT_DIR):
    import matplotlib.pyplot as plt
    config = Config(figure_size=FIGURE_SIZE, dpi=DPI, font_size=12)
    setup_plots(config)
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    ax.tick_params(labelsize=config.font_size)
    difference = 'signed' if signed else 'absolute'
    suffix = '_signed' if signed else ''
    x_column = 'aw_water_contacts' if by_water_contacts else 'solvent_exposure_percent'
    y_column = f'{difference}_volume_difference_percent'
    # Draw terminal oxygens last so they stay visible over denser groups.
    for atom_type in ['Hydrogen', 'Carbon', 'Ether oxygen', 'Terminal oxygen']:
        atom_points = points.loc[points.atom_class_aw.eq(atom_type)]
        if atom_points.empty:
            continue
        style = TYPE_STYLES[atom_type]
        terminal = atom_type == 'Terminal oxygen'
        ax.scatter(atom_points[x_column], atom_points[y_column],
                   s=1.15 * (20 if terminal else 10), alpha=.75 if terminal else .5,
                   linewidths=.3, edgecolors='#666666',
                   color=CPK_COLORS[style['element']], marker=style['marker'], label=atom_type)
        if by_water_contacts:
            contact_stats = atom_points.groupby('aw_water_contacts')[y_column].agg(['mean', 'count'])
            supported_means = contact_stats.loc[contact_stats['count'].ge(10), 'mean']
            ax.scatter(supported_means.index, supported_means,
                       s=90 if terminal else 67.5, linewidths=1.4,
                       edgecolors='white' if style['element'] == 'C' else '#222222',
                       color=CPK_COLORS[style['element']], marker=style['marker'], zorder=4)
    if by_water_contacts:
        from matplotlib.ticker import MaxNLocator

        ax.set_xlabel('Solvent contacts in AW (water molecules)')
        ax.set_title('PEG atoms: Power versus AW\nLarge markers show means for n ≥ 10')
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    else:
        ax.set_xlabel('AW atomic solvent exposure: water SA fraction (%)')
        ax.set_title(f'Solvent-exposed PEG atoms: Power versus AW\nAW water SA fraction > {100*EXPOSURE_THRESHOLD:g}%')
        ax.set_xlim(0, 100)
    ax.set_ylabel(f'{difference.title()} volume difference from AW (%)')
    if signed:
        ax.axhline(0, color='#555555', linewidth=.8)
    else:
        ax.set_ylim(bottom=0)
    legend = ax.legend(title='PEG atom type', frameon=False, markerscale=2,
                       loc='upper left', bbox_to_anchor=(1.01, 1))
    for handle in legend.legend_handles:
        handle.set_alpha(1)
    x_name = 'water_contacts' if by_water_contacts else 'exposure'
    save_figure(fig, output_dir / f'F5A_volume_vs_{x_name}{suffix}', config)


def _self_test():
    """Check threshold equality, pooling, zero denominators, and missing data."""
    pairs = pd.DataFrame({
        'polymer': ['PEG'] * 5, 'model': ['test'] * 5, 'frame': ['frame_0001'] * 5,
        'atom_index': range(5), 'atom_class_aw': ['Carbon'] * 5,
        'pair_valid': [True, True, True, False, True],
        'water_fraction_aw': [.10, .100001, np.nan, .5, .05],
        'water_coordination_aw': [0., 2., 1., 4., np.nan],
        'volume_aw': [10., 20., 30., 10., 10.], 'volume_s': [12., 10., 30., 999., 14.],
        'surface_area_aw': [4., 0., 4., 4., 4.], 'surface_area_s': [6., 3., 4., 4., np.nan],
        'contacts_aw': [0., 2., 1., 1., 2.], 'contacts_s': [4., 3., 1., 1., 2.],
    })
    summary, excluded = summarize(pairs, .10)
    rows = summary[summary.atom_type.eq('Carbon')].set_index(['metric_key', 'exposure'])
    assert np.isclose(rows.loc[('volume', 'Not Exposed'), 'mean_absolute_percentage_difference'], 30.)
    assert rows.loc[('volume', 'Not Exposed'), 'n_valid'] == 2
    assert np.isclose(rows.loc[('volume', 'Not Exposed'), 'sd_absolute_percentage_difference'], np.sqrt(200))
    assert rows.loc[('volume', 'Exposed'), 'mean_absolute_percentage_difference'] == 50.
    assert np.isnan(rows.loc[('volume', 'Exposed'), 'sd_absolute_percentage_difference'])
    assert rows.loc[('surface_area', 'Exposed'), 'n_aw_zero'] == 1
    assert np.isnan(rows.loc[('surface_area', 'Exposed'), 'mean_absolute_percentage_difference'])
    assert rows.loc[('surface_area', 'Not Exposed'), 'n_nonfinite'] == 1
    assert rows.loc[('contacts', 'Not Exposed'), 'n_aw_zero'] == 1
    assert rows.loc[('contacts', 'Exposed'), 'mean_absolute_percentage_difference'] == 50.
    assert excluded == {'invalid_or_unmatched_pairs': 1, 'invalid_aw_exposure': 1, 'outside_plot_types': 0}
    assert (summary.n_candidates == summary.n_valid + summary.n_aw_zero + summary.n_nonfinite).all()
    points = volume_points(pairs, threshold=.10)
    assert points.atom_index.tolist() == [1]
    assert np.allclose(points.solvent_exposure_percent, [10.0001])
    assert np.allclose(points.absolute_volume_difference_percent, [50.])
    signed_summary, _ = summarize(pairs, .10, signed=True)
    signed_rows = signed_summary[signed_summary.atom_type.eq('Carbon')].set_index(['metric_key', 'exposure'])
    assert signed_rows.loc[('volume', 'Exposed'), 'mean_signed_percentage_difference'] == -50.
    assert np.allclose(volume_points(pairs, .10, signed=True).signed_volume_difference_percent, [-50.])
    contact_points = volume_points(pairs, .10, signed=True, by_water_contacts=True)
    assert contact_points.atom_index.tolist() == [0, 1, 2]
    assert contact_points.aw_water_contacts.tolist() == [0, 2, 1]
    assert np.allclose(contact_points.signed_volume_difference_percent, [20., -50., 0.])
    pairs.loc[1, 'volume_aw'] = 0
    assert volume_points(pairs, threshold=.10).empty
    geometry = pd.DataFrame({
        'polymer': ['PEG'] * 10, 'model': ['EG'] * 10, 'frame': ['frame_0001'] * 10,
        'residue_id': ['r'] * 10, 'atom_index': range(10), 'repeat_count': [1] * 10,
        'order_valid': [True] * 10, 'element': ['O', 'C', 'C', 'O'] + ['H'] * 6,
        'atom_name': ['O1', 'C2', 'C3', 'O4'] + [f'H{i}' for i in range(5, 11)],
        'x': [0, 1.4, 2.8, 4.2, 1.4, 1.4, 2.8, 2.8, -.95, 5.15],
        'y': [0, 0, 0, 0, 1, -1, 1, -1, 0, 0], 'z': [0.] * 10,
    })
    roles = assign_repeat_types(geometry.sample(frac=1, random_state=0)).sort_values('atom_index')
    assert roles.atom_class.tolist() == ['Terminal oxygen', 'Carbon', 'Carbon', 'Terminal oxygen',
        'Hydrogen', 'Hydrogen', 'Hydrogen', 'Hydrogen', None, None]
    ether = geometry.loc[:3].copy()
    ether = pd.concat([ether, ether.iloc[1:].assign(atom_index=[4, 5, 6],
        atom_name=['C5', 'C6', 'O7'], x=[5.6, 7., 8.4])], ignore_index=True).assign(repeat_count=2)
    assert assign_repeat_types(ether).loc[3, 'atom_class'] == 'Ether oxygen'
    assert roles.loc[4, 'attached_heavy_atom_index'] == 1
    assert roles.type_exclusion.str.startswith('O-H').sum() == 2
    geometry.loc[4, ['x', 'y']] = [2.1, 0.]
    assert assign_repeat_types(geometry).loc[4, 'type_exclusion'] == 'Unresolved hydrogen attachment'


def main(config=None):
    """Run the checks, then export each comparison from the same matched table."""
    _self_test()
    input_table = config.output_root / 'summaries/atoms.csv' if config else INPUT_TABLE
    output_dir = config.output_root / 'F5A' if config else OUTPUT_DIR
    pairs = load_pairs(input_table, config)
    summary, excluded = summarize(pairs)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(pairs[ATOM_KEY + ['atom_name_aw', 'atom_class_aw', 'attached_heavy_atom_index', 'type_exclusion']],
              output_dir / 'F5A_atom_types.csv')
    write_csv(summary, output_dir / 'F5A_summary.csv')
    plot(summary, output_dir=output_dir)
    absolute_exposure_points = volume_points(pairs)
    write_csv(absolute_exposure_points, output_dir / 'F5A_volume_vs_exposure_atoms.csv')
    plot_volume_scatter(absolute_exposure_points, output_dir=output_dir)
    signed_summary, _ = summarize(pairs, signed=True)
    write_csv(signed_summary, output_dir / 'F5A_summary_signed.csv')
    plot(signed_summary, signed=True, output_dir=output_dir)
    signed_exposure_points = volume_points(pairs, signed=True)
    write_csv(signed_exposure_points, output_dir / 'F5A_volume_vs_exposure_atoms_signed.csv')
    plot_volume_scatter(signed_exposure_points, signed=True, output_dir=output_dir)
    contact_points = volume_points(pairs, signed=True, by_water_contacts=True)
    write_csv(contact_points, output_dir / 'F5A_volume_vs_water_contacts_atoms_signed.csv')
    plot_volume_scatter(contact_points, signed=True, by_water_contacts=True, output_dir=output_dir)
    print(f'PEG exposure scatter: {len(absolute_exposure_points):,} atom/frame points, '
          f'{len(pairs) - len(absolute_exposure_points):,} excluded.')
    print(f'PEG water-contact scatter: {len(contact_points):,} atom/frame points, '
          f'{len(pairs) - len(contact_points):,} excluded.')
    print(f'AW atomic exposure threshold: {100*EXPOSURE_THRESHOLD:g}% (strictly greater is exposed)')
    print('Pooled atom/frame observations; existing F5 primary-valid pairs only.')
    print('Pair exclusions:', excluded)
    print('Atom-type exclusions:', pairs.type_exclusion.replace('', np.nan).dropna().value_counts().to_dict())
    print(summary[['metric', 'exposure', 'atom_type', 'mean_absolute_percentage_difference',
                   'sd_absolute_percentage_difference', 'n_valid', 'n_aw_zero',
                   'n_nonfinite']].to_string(index=False, float_format=lambda v: f'{v:.3f}'))
    print(f'Saved PNG/SVG charts and summary CSV to {output_dir.resolve()}')


if __name__ == '__main__':
    main()
