"""Round-2 descriptive validation for Figure 3 PEG prototypes.

No p-values or replicate-level confidence intervals are reported: atoms and
frames within a molecular system are correlated.  Slopes and rank correlations
are sensitivity descriptions only.
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

from .F3_PEG_prototypes import (
    ATOM_CLASSES, METRICS, METRIC_LABELS, MODELS, NEUTRAL, REPO_ROOT,
    collect_observations, exposure_tables, length_tables,
)
from vorpy.src.analyze.Part_II_Molecular_Analysis.Shared_Utilities.simple_models.F5_common import save_figure, setup_plots
from vorpy.src.analyze.Part_II_Molecular_Analysis.Shared_Utilities.simple_models.F5_config import Config


OUTPUT_DIR = REPO_ROOT / 'output' / 'figure3_round2_validation'
PEG_MODELS = tuple(model for model in MODELS if model != 'EG')


def descriptive_trend(x, y):
    """Return a descriptive linear slope and Spearman rho without inference."""
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x)[valid], np.asarray(y)[valid]
    if len(x) < 2 or np.unique(x).size < 2:
        return np.nan, np.nan, len(x)
    slope = np.polyfit(x, y, 1)[0]
    rho = pd.Series(x).corr(pd.Series(y), method='spearman')
    return slope, rho, len(x)


def length_weighting_audit(observations):
    frame, model = length_tables(observations)
    rows = []
    for metric in METRICS:
        for weighting, table, value in [
            ('equal_frame', frame, 'mean_absolute_difference_pct'),
            ('equal_system', model, 'mean_absolute_difference_pct'),
        ]:
            selected = table.loc[table.metric.eq(metric)]
            slope, rho, count = descriptive_trend(
                np.log2(selected.repeat_units), selected[value])
            rows.append({
                'metric': metric, 'weighting': weighting,
                'slope_pct_points_per_doubling': slope,
                'spearman_rho': rho, 'n_units': count,
                'unit_definition': 'correlated frame mean' if weighting == 'equal_frame'
                                   else 'molecular-system mean',
            })
    return frame, model, pd.DataFrame(rows)


def within_environment_audit(observations):
    chemical_frame = (observations.groupby(
        ['model', 'frame', 'repeat_units', 'metric', 'atom_class_aw'], observed=True)
        .agg(mean_absolute_difference_pct=('abs_delta_pct', 'mean'),
             n_atom_observations=('atom_index', 'size')).reset_index())
    chemical_model = (chemical_frame.groupby(
        ['model', 'repeat_units', 'metric', 'atom_class_aw'], observed=True)
        .agg(mean_absolute_difference_pct=('mean_absolute_difference_pct', 'mean'),
             n_frames=('frame', 'nunique'),
             n_atom_frame_observations=('n_atom_observations', 'sum')).reset_index())

    _, exposure_model, _ = exposure_tables(observations)
    repeat_units = observations[['model', 'repeat_units']].drop_duplicates()
    exposure_model = exposure_model.merge(repeat_units, on='model', validate='many_to_one')

    rows = []
    for metric in METRICS:
        for atom_class in ATOM_CLASSES:
            group = chemical_model[
                chemical_model.metric.eq(metric) & chemical_model.atom_class_aw.eq(atom_class)]
            slope, rho, count = descriptive_trend(
                np.log2(group.repeat_units), group.mean_absolute_difference_pct)
            rows.append(dict(metric=metric, grouping='chemical_environment',
                             environment=atom_class, slope_pct_points_per_doubling=slope,
                             spearman_rho=rho, n_systems=count))
        for center, group in exposure_model.loc[
                exposure_model.metric.eq(metric)].groupby('exposure_bin_center_pct'):
            slope, rho, count = descriptive_trend(
                np.log2(group.repeat_units), group.mean_absolute_difference_pct)
            rows.append(dict(metric=metric, grouping='aw_exposure_bin',
                             environment=f'{center:g}%', slope_pct_points_per_doubling=slope,
                             spearman_rho=rho, n_systems=count))
    return chemical_frame, chemical_model, exposure_model, pd.DataFrame(rows)


def standardized_length_audit(observations, length_model, chemical_model, exposure_model):
    """Standardize PEG-only models to common chemistry and exposure mixtures."""
    rows = []
    peg = observations[observations.model.isin(PEG_MODELS)]
    class_counts = (peg.groupby(['model', 'metric', 'atom_class_aw'], observed=True)
                    .size().rename('n').reset_index())
    class_counts['weight'] = class_counts.n / class_counts.groupby(
        ['model', 'metric'], observed=True).n.transform('sum')
    reference_class_weights = (class_counts.groupby(['metric', 'atom_class_aw'], observed=True)
                               .weight.mean().rename('reference_weight').reset_index())

    for metric in METRICS:
        weights = reference_class_weights[
            reference_class_weights.metric.eq(metric)].set_index('atom_class_aw').reference_weight
        if set(weights.index) != set(ATOM_CLASSES):
            raise ValueError('PEG-only chemistry standardization requires all four atom classes.')
        for model in PEG_MODELS:
            values = chemical_model[
                chemical_model.metric.eq(metric) & chemical_model.model.eq(model)
            ].set_index('atom_class_aw').mean_absolute_difference_pct
            rows.append(dict(metric=metric, model=model, adjustment='chemistry_standardized',
                             mean_absolute_difference_pct=(values * weights).sum(),
                             support='four PEG atom classes; mean PEG-system class weights'))

        metric_bins = exposure_model[
            exposure_model.metric.eq(metric) & exposure_model.model.isin(PEG_MODELS)]
        common_bins = set.intersection(*(
            set(metric_bins.loc[metric_bins.model.eq(model), 'exposure_bin_center_pct'])
            for model in PEG_MODELS))
        common = metric_bins[metric_bins.exposure_bin_center_pct.isin(common_bins)].copy()
        common['model_weight'] = common.n_atom_frame_observations / common.groupby(
            'model', observed=True).n_atom_frame_observations.transform('sum')
        reference = common.groupby('exposure_bin_center_pct', observed=True).model_weight.mean()
        for model in PEG_MODELS:
            values = common.loc[common.model.eq(model)].set_index(
                'exposure_bin_center_pct').mean_absolute_difference_pct
            rows.append(dict(metric=metric, model=model, adjustment='exposure_standardized',
                             mean_absolute_difference_pct=(values * reference).sum(),
                             support='common exposure bins; mean PEG-system bin weights'))

    adjusted = pd.DataFrame(rows).merge(
        observations[['model', 'repeat_units']].drop_duplicates(), on='model', validate='many_to_one')
    observed = length_model[length_model.model.isin(PEG_MODELS)][
        ['metric', 'model', 'repeat_units', 'mean_absolute_difference_pct']].copy()
    observed['adjustment'] = 'observed'
    observed['support'] = 'all eligible atoms'
    combined = pd.concat([observed, adjusted], ignore_index=True)

    trend_rows = []
    for (metric, adjustment), group in combined.groupby(['metric', 'adjustment'], observed=True):
        slope, rho, count = descriptive_trend(
            np.log2(group.repeat_units), group.mean_absolute_difference_pct)
        trend_rows.append(dict(metric=metric, adjustment=adjustment,
                               slope_pct_points_per_doubling=slope,
                               spearman_rho=rho, n_systems=count))
    return combined, pd.DataFrame(trend_rows), reference_class_weights


def plot_standardized_trends(combined, output_dir):
    config = Config(figure_size=(10, 4.2), dpi=240, font_size=10, export_formats=('png', 'svg'))
    setup_plots(config)
    import matplotlib.pyplot as plt
    styles = {
        'observed': dict(color='#222222', marker='o', linestyle='-'),
        'chemistry_standardized': dict(color='#7b5ea7', marker='s', linestyle='--'),
        'exposure_standardized': dict(color='#2b8c7e', marker='^', linestyle=':'),
    }
    fig, axes = plt.subplots(1, 2, figsize=config.figure_size)
    for ax, metric in zip(axes, METRICS):
        for adjustment, style in styles.items():
            group = combined[
                combined.metric.eq(metric) & combined.adjustment.eq(adjustment)
            ].sort_values('repeat_units')
            ax.plot(group.repeat_units, group.mean_absolute_difference_pct,
                    label=adjustment.replace('_', ' '), **style)
        ax.set_xscale('log', base=2)
        ticks = combined.loc[combined.metric.eq(metric), ['repeat_units', 'model']].drop_duplicates()
        ax.set_xticks(ticks.repeat_units, ticks.model, rotation=35)
        ax.set(title=METRIC_LABELS[metric], xlabel='PEG repeat units, n',
               ylabel='Mean |Power − AW| / AW (%)')
    axes[1].legend(frameon=False)
    fig.suptitle('PEG length trends after descriptive composition standardization')
    save_figure(fig, output_dir / 'F3_round2_standardized_length_trends', config)


def plot_within_environment_slopes(audit, output_dir):
    config = Config(figure_size=(10, 4.2), dpi=240, font_size=10, export_formats=('png', 'svg'))
    setup_plots(config)
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=config.figure_size)
    for ax, metric in zip(axes, METRICS):
        data = audit[
            audit.metric.eq(metric) & audit.grouping.eq('chemical_environment')
        ].set_index('environment').reindex(ATOM_CLASSES)
        ax.bar(data.index, data.slope_pct_points_per_doubling, color=NEUTRAL)
        ax.axhline(0, color='#777777', linewidth=.8)
        ax.tick_params(axis='x', rotation=25)
        ax.set(title=METRIC_LABELS[metric],
               ylabel='Descriptive slope (% points per doubling)')
    fig.suptitle('Within-chemical-environment PEG length trends')
    save_figure(fig, output_dir / 'F3_round2_within_chemistry_slopes', config)


def main(output_dir=OUTPUT_DIR):
    output_dir = Path(output_dir)
    observations, coverage = collect_observations()
    length_frame, length_model, weighting = length_weighting_audit(observations)
    chemical_frame, chemical_model, exposure_model, within = within_environment_audit(observations)
    standardized, standardized_trends, reference_weights = standardized_length_audit(
        observations, length_model, chemical_model, exposure_model)
    output_dir.mkdir(parents=True, exist_ok=True)
    tables = {
        'F3_round2_coverage.csv': coverage,
        'F3_round2_length_weighting_audit.csv': weighting,
        'F3_round2_within_environment_audit.csv': within,
        'F3_round2_standardized_model_means.csv': standardized,
        'F3_round2_standardized_trends.csv': standardized_trends,
        'F3_round2_reference_chemistry_weights.csv': reference_weights,
        'F3_round2_chemical_frame_means.csv': chemical_frame,
    }
    for name, table in tables.items():
        table.to_csv(output_dir / name, index=False)
    plot_standardized_trends(standardized, output_dir)
    plot_within_environment_slopes(within, output_dir)
    print(f'Validated Figure 3 weighting and environment trends across {len(coverage)} frames.')


if __name__ == '__main__':
    main()
