"""Paired scheme deviations against AW, in native units or percentage points."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if not __package__:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[5]))
    __package__ = 'vorpy.src.analyze.Part_II_Molecular_Analysis.F5'

from .F5_config import KEY, METRICS, UNITS, SCHEME_LABELS
from .F5_common import write_csv, save_figure, empty_panel, distribution_panels


def run(data,config):
    pairs=data['atom_pairs'];unit_pairs=data['unit_pairs']
    environment_rows=[];bin_rows=[]
    for level,table in [('atom',pairs),('repeat',unit_pairs),('residue',data.get('residue_pairs',pd.DataFrame()))]:
        if table.empty:continue
        valid=table[table.pair_valid].copy()
        valid['exposure_bin']=pd.cut(valid.water_fraction_aw,np.linspace(0,1,11),include_lowest=True,
            labels=(np.arange(10)+.5)/10).astype(float)
        valid['environment']=valid.environment_aw
        for metric in METRICS:
            col='delta_'+metric
            for metadata,group in valid.groupby(KEY+['scheme','environment'],observed=True,dropna=False):
                values=group[col].dropna()
                environment_rows.append(dict(zip(KEY+['scheme','environment'],metadata),level=level,metric=metric,
                    n_matched=len(group),n_valid=len(values),mean_signed_deviation=values.mean(),
                    mean_absolute_deviation=values.abs().mean(),median_absolute_deviation=values.abs().median()))
            for metadata,group in valid.groupby(KEY+['scheme','exposure_bin'],observed=True):
                values=group[col].dropna()
                bin_rows.append(dict(zip(KEY+['scheme','exposure_bin'],metadata),level=level,metric=metric,
                    n_matched=len(group),n_valid=len(values),mean_signed_deviation=values.mean(),
                    mean_absolute_deviation=values.abs().mean(),median_absolute_deviation=values.abs().median()))
    environments=pd.DataFrame(environment_rows,columns=KEY+['scheme','environment','level','metric','n_matched','n_valid',
        'mean_signed_deviation','mean_absolute_deviation','median_absolute_deviation'])
    bins=pd.DataFrame(bin_rows,columns=KEY+['scheme','exposure_bin','level','metric','n_matched','n_valid',
        'mean_signed_deviation','mean_absolute_deviation','median_absolute_deviation'])
    write_csv(environments,config.output_root/'summaries/F5E_environment_deviations.csv')
    write_csv(bins,config.output_root/'summaries/F5E_exposure_binned_deviations.csv')
    model_rows=[]
    for metadata,group in environments.groupby(['polymer','model','scheme','level','metric','environment'],dropna=False):
        model_rows.append(dict(zip(['polymer','model','scheme','level','metric','environment'],metadata),
            n_frames=group.frame.nunique(),n_matched=int(group.n_matched.sum()),n_valid=int(group.n_valid.sum()),
            mean_signed_deviation=group.mean_signed_deviation.mean(),
            frame_sd_signed_deviation=group.mean_signed_deviation.std(),
            mean_absolute_deviation=group.mean_absolute_deviation.mean(),
            mean_frame_median_absolute_deviation=group.median_absolute_deviation.mean()))
    write_csv(pd.DataFrame(model_rows),config.output_root/'summaries/F5E_model_deviation_summary.csv')
    for family in config.polymers:
        for level in ['atom','repeat','residue']:
            fig,axes=plt.subplots(2,3,figsize=config.figure_size)
            for ax,metric in zip(axes.flat,METRICS):
                found=False
                for scheme in ['pow','prm']:
                    subset=bins[bins.polymer.eq(family)&bins.level.eq(level)&bins.scheme.eq(scheme)&bins.metric.eq(metric)]
                    if subset.empty:continue
                    # Average frames within a model first; give each model equal weight.
                    models=subset.groupby(['model','exposure_bin'],observed=True).mean_signed_deviation.mean().reset_index()
                    summary=models.groupby('exposure_bin').mean_signed_deviation.agg(['mean','std']).reset_index()
                    if summary['mean'].notna().any():found=True
                    x=100*summary.exposure_bin.to_numpy(dtype=float);y=summary['mean'].to_numpy();sd=summary['std'].to_numpy()
                    ax.plot(x,y,color=dict(config.colors)[scheme],marker='o',ms=3,label=SCHEME_LABELS[scheme])
                    ax.fill_between(x,y-sd,y+sd,color=dict(config.colors)[scheme],alpha=.12)
                ax.axhline(0,color='gray',lw=.7)
                ax.set_xlabel('AW water surface-area fraction (%)')
                label='Δ water SA fraction (percentage points)' if metric=='water_fraction' else 'Δ '+UNITS[metric]
                if level=='repeat' and metric=='volume':label='Δ Repeat volume (Å³)'
                if level=='repeat' and metric=='surface_area':label='Δ Sum of atomic SA / repeat (Å²)'
                if level=='residue' and metric=='volume':label='Δ PDB residue volume (Å³)'
                if level=='residue' and metric=='surface_area':label='Δ Sum of atomic SA / residue (Å²)'
                ax.set_ylabel(label)
                if not found:empty_panel(ax)
            handles,labels=axes.flat[0].get_legend_handles_labels()
            if handles:fig.legend(handles,labels,loc='lower center',ncol=2)
            fig.suptitle(f'{family}: paired {level} deviations from AW\nEqual model weighting; bands: between-model SD after averaging frames')
            save_figure(fig,config.output_root/family/'F5E'/f'{level}_deviations_vs_exposure',config)
        # Distribution of matched atom frame means, using AW repeat environments.
        subset=environments[environments.polymer.eq(family)&environments.level.eq('atom')]
        wide=subset.pivot_table(index=KEY+['scheme','environment'],columns='metric',values='mean_signed_deviation',observed=True).reset_index()
        if wide.empty:
            wide=pd.DataFrame(columns=KEY+['scheme','environment']+METRICS)
        for metric in METRICS:
            if metric not in wide:wide[metric]=np.nan
        # Rename to communicate native difference units in the shared distribution helper.
        from .F5_config import UNITS as labels
        delta_names={m:'delta_'+m for m in METRICS}
        for metric,new in delta_names.items():
            labels[new]='Δ water fraction (percentage points)' if metric=='water_fraction' else 'Δ '+labels[metric]
        distribution_panels(wide.rename(columns=delta_names),list(delta_names.values()),
            f'{family}: paired atom deviations by AW reference environment',
            config.output_root/family/'F5E/deviations_by_environment',config)


if __name__=='__main__':
    from .F5_run_all import main
    main(default_figures=['E'])
