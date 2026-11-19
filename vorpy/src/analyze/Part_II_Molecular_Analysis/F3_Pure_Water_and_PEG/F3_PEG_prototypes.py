"""Generate Figure 3 PEG prototypes from the screened 69-frame cohort.

The signed atomic difference is 100 * (Power - AW) / AW.  Length and
exposure plots use its magnitude; the chemistry plot keeps the sign.  Atom
observations are reduced to frame means and then model means so long PEG
chains do not dominate merely because they contain more atoms.  Frames from
one trajectory are correlated; the seven molecular models are the displayed
system-level replicates.

This reader reuses the validated F5 PEG loaders.  It does not calculate new
Voronoi geometry or overwrite the validated F5 outputs.
"""
from pathlib import Path

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

from .F5A_PEG_chemistry_vs_hydration import MODELS, TYPE_STYLES, load_frame, plot_data
from vorpy.src.analyze.Part_II_Molecular_Analysis.Shared_Utilities.simple_models.F5_common import save_figure, setup_plots
from vorpy.src.analyze.Part_II_Molecular_Analysis.Shared_Utilities.simple_models.F5_config import ATOM_KEY, Config
try:
    from .support.F5_peg_dataset import screened_runs
except ModuleNotFoundError:
    from .F5_peg_dataset import screened_runs


REPO_ROOT = next(parent for parent in Path(__file__).resolve().parents
                 if (parent / 'pyproject.toml').is_file())
OUTPUT_DIR = REPO_ROOT / 'output' / 'figure3_prototypes'
METRICS = ('volume', 'surface_area')
METRIC_LABELS = {'volume': 'Atomic volume', 'surface_area': 'Atomic surface area'}
ATOM_CLASSES = tuple(TYPE_STYLES)
NEUTRAL = '#59636e'
EXPOSURE_EDGES = np.arange(0, 101, 10)
MIN_SYSTEMS_FOR_EXPOSURE_PLOT = 3


def infer_repeat_units(matched_atom_count):
    """Infer n for HO-(CH2CH2O)n-H from the complete matched atom count."""
    units, remainder = divmod(int(matched_atom_count) - 3, 7)
    if units < 1 or remainder:
        raise ValueError(f'Atom count {matched_atom_count} is not a complete PEG formula (7n + 3).')
    return units


def collect_observations(config=None):
    """Load all screened frames and return metric-long atom observations."""
    config = config or Config()
    runs = screened_runs()
    observations = []
    coverage = []
    for model in MODELS:
        frames = sorted(runs.loc[runs.molecule.eq(model), 'frame'].unique())
        for frame in frames:
            pairs, reason, invalid = load_frame(model, frame, config)
            if pairs is None:
                raise ValueError(f'Screened frame unavailable: {model}/{frame}: {reason}')
            matched_atoms = int(pairs.pair_valid.sum())
            repeat_units = infer_repeat_units(matched_atoms)
            exposure = pairs[ATOM_KEY + ['water_fraction_aw']]
            frame_counts = {}
            invalid_exposure = 0
            for metric in METRICS:
                points, _ = plot_data(pairs, metric)
                points = points.merge(exposure, on=ATOM_KEY, how='left', validate='one_to_one')
                points['metric'] = metric
                points['delta_pct'] = points.pop(f'delta_{metric}_pct')
                points['abs_delta_pct'] = points.delta_pct.abs()
                points['repeat_units'] = repeat_units
                valid_exposure = np.isfinite(points.water_fraction_aw) & points.water_fraction_aw.between(0, 1)
                invalid_exposure += int((~valid_exposure).sum())
                points['exposure_valid'] = valid_exposure
                observations.append(points)
                frame_counts[metric] = len(points)
            coverage.append({
                'model': model, 'frame': frame, 'repeat_units': repeat_units,
                'matched_polymer_atoms': matched_atoms,
                'layer1_water_molecules': invalid['layer1_waters'],
                'volume_atom_observations': frame_counts['volume'],
                'surface_area_atom_observations': frame_counts['surface_area'],
                'invalid_exposure_observations_across_metrics': invalid_exposure,
            })
    coverage = pd.DataFrame(coverage)
    if len(coverage) != 69:
        raise ValueError(f'Expected 69 screened model/frame sets, found {len(coverage)}.')
    return pd.concat(observations, ignore_index=True), coverage


def length_tables(observations):
    frame = (observations.groupby(['model', 'frame', 'repeat_units', 'metric'], observed=True)
             .agg(mean_absolute_difference_pct=('abs_delta_pct', 'mean'),
                  mean_signed_difference_pct=('delta_pct', 'mean'),
                  n_atom_observations=('atom_index', 'size')).reset_index())
    model = (frame.groupby(['model', 'repeat_units', 'metric'], observed=True)
             .agg(mean_absolute_difference_pct=('mean_absolute_difference_pct', 'mean'),
                  sd_frame_means_pct=('mean_absolute_difference_pct', 'std'),
                  n_frames=('frame', 'nunique'),
                  n_atom_frame_observations=('n_atom_observations', 'sum')).reset_index())
    return frame, model


def chemistry_tables(observations):
    frame = (observations.groupby(['model', 'frame', 'metric', 'atom_class_aw'], observed=True)
             .agg(mean_signed_difference_pct=('delta_pct', 'mean'),
                  n_atom_observations=('atom_index', 'size')).reset_index())
    model = (frame.groupby(['model', 'metric', 'atom_class_aw'], observed=True)
             .agg(mean_signed_difference_pct=('mean_signed_difference_pct', 'mean'),
                  sd_frame_means_pct=('mean_signed_difference_pct', 'std'),
                  n_frames=('frame', 'nunique'),
                  n_atom_frame_observations=('n_atom_observations', 'sum')).reset_index())
    summary = (model.groupby(['metric', 'atom_class_aw'], observed=True)
               .agg(equal_system_mean_pct=('mean_signed_difference_pct', 'mean'),
                    sd_across_systems_pct=('mean_signed_difference_pct', 'std'),
                    n_systems=('model', 'nunique'),
                    n_atom_frame_observations=('n_atom_frame_observations', 'sum')).reset_index())
    return frame, model, summary


def exposure_tables(observations):
    eligible = observations.loc[observations.exposure_valid].copy()
    eligible['aw_water_exposed_area_pct'] = 100 * eligible.water_fraction_aw
    eligible['exposure_bin_center_pct'] = pd.cut(
        eligible.aw_water_exposed_area_pct, EXPOSURE_EDGES,
        labels=(EXPOSURE_EDGES[:-1] + EXPOSURE_EDGES[1:]) / 2,
        include_lowest=True,
    ).astype(float)
    frame = (eligible.groupby(['model', 'frame', 'metric', 'exposure_bin_center_pct'], observed=True)
             .agg(mean_absolute_difference_pct=('abs_delta_pct', 'mean'),
                  n_atom_observations=('atom_index', 'size')).reset_index())
    model = (frame.groupby(['model', 'metric', 'exposure_bin_center_pct'], observed=True)
             .agg(mean_absolute_difference_pct=('mean_absolute_difference_pct', 'mean'),
                  n_frames=('frame', 'nunique'),
                  n_atom_frame_observations=('n_atom_observations', 'sum')).reset_index())
    summary = (model.groupby(['metric', 'exposure_bin_center_pct'], observed=True)
               .agg(equal_system_mean_pct=('mean_absolute_difference_pct', 'mean'),
                    sd_across_systems_pct=('mean_absolute_difference_pct', 'std'),
                    n_systems=('model', 'nunique'),
                    n_atom_frame_observations=('n_atom_frame_observations', 'sum')).reset_index())
    return frame, model, summary


def _plot_config():
    return Config(figure_size=(10, 4.2), dpi=240, font_size=10, export_formats=('png', 'svg'))


def plot_length(frame, model, output_dir):
    config = _plot_config()
    setup_plots(config)
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=config.figure_size, sharey=False)
    for ax, metric in zip(axes, METRICS):
        points = frame.loc[frame.metric.eq(metric)]
        means = model.loc[model.metric.eq(metric)].sort_values('repeat_units')
        ax.scatter(points.repeat_units, points.mean_absolute_difference_pct,
                   color=NEUTRAL, alpha=.18, s=16, linewidths=0)
        ax.errorbar(means.repeat_units, means.mean_absolute_difference_pct,
                    yerr=means.sd_frame_means_pct, color=NEUTRAL, marker='o', capsize=3)
        ax.set_xscale('log', base=2)
        ax.set_xticks(means.repeat_units, means.model)
        ax.tick_params(axis='x', rotation=35)
        ax.set(title=METRIC_LABELS[metric], xlabel='PEG repeat units, n',
               ylabel='Mean |Power − AW| / AW (%)')
    fig.suptitle('PEG length and atomic Voronoi-scheme differences')
    save_figure(fig, output_dir / 'F3_prototype_PEG_length', config)


def plot_chemistry(model, summary, output_dir):
    config = _plot_config()
    setup_plots(config)
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=config.figure_size)
    system_offsets = dict(zip(MODELS, np.linspace(-.24, .24, len(MODELS))))
    for ax, metric in zip(axes, METRICS):
        values = model.loc[model.metric.eq(metric)]
        totals = summary.loc[summary.metric.eq(metric)].set_index('atom_class_aw').reindex(ATOM_CLASSES)
        for position, atom_class in enumerate(ATOM_CLASSES):
            group = values.loc[values.atom_class_aw.eq(atom_class)]
            x = [position + system_offsets[name] for name in group.model]
            ax.scatter(x, group.mean_signed_difference_pct, color=NEUTRAL,
                       alpha=.6, s=22, linewidths=0)
        x = np.arange(len(ATOM_CLASSES))
        ax.errorbar(x, totals.equal_system_mean_pct, yerr=totals.sd_across_systems_pct,
                    color='#222222', marker='D', linestyle='none', capsize=3, zorder=3)
        ax.axhline(0, color='#777777', linewidth=.8)
        ax.set_xticks(x, ATOM_CLASSES, rotation=25, ha='right')
        ax.set(title=METRIC_LABELS[metric], ylabel='100 × (Power − AW) / AW (%)')
    fig.suptitle('Chemical environment and signed atomic differences')
    save_figure(fig, output_dir / 'F3_prototype_chemical_environment', config)


def plot_exposure(model, summary, output_dir):
    config = _plot_config()
    setup_plots(config)
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=config.figure_size)
    for ax, metric in zip(axes, METRICS):
        values = model.loc[model.metric.eq(metric)]
        totals = summary.loc[
            summary.metric.eq(metric) & summary.n_systems.ge(MIN_SYSTEMS_FOR_EXPOSURE_PLOT)
        ].sort_values('exposure_bin_center_pct')
        supported_bins = set(totals.exposure_bin_center_pct)
        for _, group in values.groupby('model', sort=False):
            group = group.loc[group.exposure_bin_center_pct.isin(supported_bins)].sort_values(
                'exposure_bin_center_pct')
            ax.plot(group.exposure_bin_center_pct, group.mean_absolute_difference_pct,
                    color=NEUTRAL, alpha=.20, linewidth=.8)
        ax.errorbar(totals.exposure_bin_center_pct, totals.equal_system_mean_pct,
                    yerr=totals.sd_across_systems_pct, color=NEUTRAL,
                    marker='o', capsize=3)
        for threshold in (10, 40):
            ax.axvline(threshold, color='#999999', linestyle='--', linewidth=.8)
        ax.set(xlim=(0, 100), title=METRIC_LABELS[metric],
               xlabel='AW water-exposed atomic surface area (%)',
               ylabel='Mean |Power − AW| / AW (%)')
    fig.suptitle('Solvent exposure and atomic Voronoi-scheme differences')
    save_figure(fig, output_dir / 'F3_prototype_solvent_exposure', config)


def main(output_dir=OUTPUT_DIR):
    output_dir = Path(output_dir)
    observations, coverage = collect_observations()
    length_frame, length_model = length_tables(observations)
    chemistry_frame, chemistry_model, chemistry_summary = chemistry_tables(observations)
    exposure_frame, exposure_model, exposure_summary = exposure_tables(observations)

    output_dir.mkdir(parents=True, exist_ok=True)
    tables = {
        'F3_prototype_coverage.csv': coverage,
        'F3_prototype_length_frame_means.csv': length_frame,
        'F3_prototype_length_model_summary.csv': length_model,
        'F3_prototype_chemistry_frame_means.csv': chemistry_frame,
        'F3_prototype_chemistry_model_means.csv': chemistry_model,
        'F3_prototype_chemistry_summary.csv': chemistry_summary,
        'F3_prototype_exposure_frame_bin_means.csv': exposure_frame,
        'F3_prototype_exposure_model_bin_means.csv': exposure_model,
        'F3_prototype_exposure_summary.csv': exposure_summary,
    }
    for name, table in tables.items():
        table.to_csv(output_dir / name, index=False)
    plot_length(length_frame, length_model, output_dir)
    plot_chemistry(chemistry_model, chemistry_summary, output_dir)
    plot_exposure(exposure_model, exposure_summary, output_dir)
    print(f'Generated three Figure 3 prototypes from {len(coverage)} screened frames in {output_dir}')
    print(f'Atom/frame observations: {len(observations) // len(METRICS)} per metric; '
          f'{coverage.matched_polymer_atoms.sum()} matched atoms before chemistry exclusions')


if __name__ == '__main__':
    main()
