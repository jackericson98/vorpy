"""Terminal profiles of validated repeat units (no convergence-model fitting)."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if not __package__:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[5]))
    __package__ = 'vorpy.src.analyze.Part_II_Molecular_Analysis.F5'

from .F5_config import KEY, UNITS, SCHEME_LABELS
from .F5_common import frame_summary, temporal_summary, write_csv, save_figure, empty_panel


def run(data, config):
    metrics=['water_fraction','mean_atomic_volume','water_coordination','polymer_coordination']
    units=data['units']
    frames=frame_summary(units[units.order_valid],KEY+['scheme','end_distance','repeat_count'],metrics)
    write_csv(frames,config.output_root/'summaries/F5A_frame_profiles.csv')
    summary=temporal_summary(frames,['polymer','model','scheme','end_distance','repeat_count'],metrics)
    write_csv(summary,config.output_root/'summaries/F5A_terminal_profiles.csv')
    contrasts=[]
    for metadata,group in units[units.primary_valid & units.order_valid].groupby(KEY+['scheme','repeat_count']):
        n=metadata[-1]
        if n<4:continue
        terminal=group[group.end_distance.eq(0)]
        interior=group[group.end_distance.ge(n/4)]
        for metric in metrics:
            contrasts.append(dict(zip(KEY+['scheme','repeat_count'],metadata),metric=metric,
                terminal_mean=terminal[metric].mean(),interior_mean=interior[metric].mean(),
                terminal_minus_interior=terminal[metric].mean()-interior[metric].mean(),
                n_terminal_units=len(terminal),n_interior_units=len(interior)))
    write_csv(pd.DataFrame(contrasts),config.output_root/'summaries/F5A_terminal_interior_frame_contrasts.csv')
    markers=['o','s','^','v','D','x','+']
    for family in config.polymers:
        fig,axes=plt.subplots(4,len(config.schemes),figsize=config.figure_size,squeeze=False,sharex='col')
        family_data=summary[summary.polymer.eq(family)] if not summary.empty else summary
        models=family_data.sort_values('repeat_count').model.unique() if not family_data.empty else []
        handles=[]
        for k,scheme in enumerate(config.schemes):
            color=dict(config.colors)[scheme]
            for row,metric in enumerate(metrics):
                ax=axes[row,k]
                for j,model in enumerate(models):
                    group=family_data[family_data.model.eq(model)&family_data.scheme.eq(scheme)].sort_values('end_distance')
                    if group.empty:continue
                    x=group.end_distance.to_numpy(dtype=float)
                    y=group[f'{metric}_mean'].to_numpy(dtype=float)
                    sd=group[f'{metric}_std'].to_numpy(dtype=float)
                    line,=ax.plot(x,y,color=color,marker=markers[j%len(markers)],markersize=3,
                        linestyle=['-','--',':','-.'][j%4],label=model)
                    ax.fill_between(x,y-sd,y+sd,color=color,alpha=.05)
                    if row==0 and k==0:handles.append(line)
                if len(models)==0:empty_panel(ax,'No validated repeat ordering')
                if row==0:ax.set_title(SCHEME_LABELS[scheme],color=color)
                if k==0:ax.set_ylabel({'water_fraction':'Water SA\nfraction',
                    'mean_atomic_volume':'Mean atomic\nvolume (Å³)',
                    'water_coordination':'Water molecules\nper repeat',
                    'polymer_coordination':'Polymer atom\nneighbors per repeat'}[metric])
                if row==3:ax.set_xlabel('Repeats from nearest chain end')
        if handles:fig.legend(handles,[h.get_label() for h in handles],loc='lower center',ncol=min(7,len(handles)))
        fig.suptitle(f'{family}: terminal profiles (bands: frame SD, descriptive)')
        save_figure(fig,config.output_root/family/'F5A/terminal_effects',config)


if __name__=='__main__':
    from .F5_run_all import main
    main(default_figures=['A'])
