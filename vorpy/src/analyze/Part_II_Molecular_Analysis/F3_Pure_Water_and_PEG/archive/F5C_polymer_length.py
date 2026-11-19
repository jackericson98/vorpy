"""Length dependence; error bars are temporal SD, not independent-replicate CIs."""
import pandas as pd
import matplotlib.pyplot as plt

if not __package__:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[5]))
    __package__ = 'vorpy.src.analyze.Part_II_Molecular_Analysis.F5'

from .F5_config import KEY, UNITS, SCHEME_LABELS
from .F5_common import exposure, temporal_summary, write_csv, save_figure, empty_panel


def run(data,config):
    units=data['units']; rows=[]
    metrics=['water_fraction','buried_fraction','mean_atomic_volume','water_coordination','polymer_coordination']
    for metadata,group in units[units.primary_valid & units.order_valid].groupby(KEY+['scheme','repeat_count']):
        row=dict(zip(KEY+['scheme','repeat_count'],metadata))
        row.update(n_units=len(group),n_atoms=int(group.atom_count.sum()),
            water_fraction=group.water_area.sum()/group.surface_area.sum(),
            buried_fraction=exposure(group.water_fraction,config).eq('Buried').mean(),
            mean_atomic_volume=group.volume.sum(min_count=len(group))/group.atom_count.sum(),
            water_coordination=group.water_coordination.mean(),polymer_coordination=group.polymer_coordination.mean())
        rows.append(row)
    frames=pd.DataFrame(rows,columns=KEY+['scheme','repeat_count','n_units','n_atoms']+metrics)
    write_csv(frames,config.output_root/'summaries/F5C_length_frame_means.csv')
    summary=temporal_summary(frames,['polymer','model','scheme','repeat_count'],metrics)
    write_csv(summary,config.output_root/'summaries/F5C_length_summary.csv')
    for family in config.polymers:
        fig,axes=plt.subplots(2,3,figsize=config.figure_size)
        for ax,metric in zip(axes.flat,metrics):
            found=False
            for scheme in config.schemes:
                group=summary[summary.polymer.eq(family)&summary.scheme.eq(scheme)].sort_values('repeat_count') if not summary.empty else summary
                if group.empty:continue
                found=True
                ax.plot(group.repeat_count,group[f'{metric}_mean'],
                    color=dict(config.colors)[scheme],marker='o',ms=4,label=SCHEME_LABELS[scheme])
                measured=group[group[f'{metric}_std'].notna()]
                ax.errorbar(measured.repeat_count,measured[f'{metric}_mean'],yerr=measured[f'{metric}_std'],
                    color=dict(config.colors)[scheme],fmt='none',capsize=2)
                if scheme==config.schemes[0]:
                    for row in group.itertuples():
                        ax.annotate(row.model,(row.repeat_count,getattr(row,f'{metric}_mean')),fontsize=6,xytext=(2,5),textcoords='offset points')
            ax.set_xlabel('Repeat-unit count, N')
            ax.set_ylabel(UNITS[metric])
            if not found:empty_panel(ax)
        axes.flat[-1].set_visible(False)
        handles,labels=axes.flat[0].get_legend_handles_labels()
        if handles:fig.legend(handles,labels,loc='lower center',ncol=3)
        fig.suptitle(f'{family}: length dependence (error bars: frame SD)\nBuried fraction uses each scheme’s own provisional exposure bins')
        save_figure(fig,config.output_root/family/'F5C/polymer_length',config)


if __name__=='__main__':
    from .F5_run_all import main
    main(default_figures=['C'])
