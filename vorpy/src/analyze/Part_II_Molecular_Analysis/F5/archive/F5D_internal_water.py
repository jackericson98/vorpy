"""Available boundary diagnostics only; internal occupancy is explicitly unavailable."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if not __package__:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[5]))
    __package__ = 'vorpy.src.analyze.Part_II_Molecular_Analysis.F5'

from .F5_config import KEY, SCHEME_LABELS
from .F5_common import write_csv, save_figure, empty_panel


def run(data,config):
    topology=data['topology'].copy()
    topology['conditional_2pi_chi']=np.where(topology.gauss_bonnet_applicable,2*np.pi*pd.to_numeric(topology['Euler Characteristic'],errors='coerce'),np.nan)
    topology['frame_number']=topology.frame.str.extract(r'(\d+)$').astype(int)
    write_csv(topology,config.output_root/'summaries/F5D_topology_diagnostics.csv')
    write_csv(topology[KEY+['scheme','internal_occupancy_status']],config.output_root/'summaries/F5D_internal_occupancy_unavailable.csv')
    markers=['o','s','^','v','D','x','+']
    for family in config.polymers:
        subset=topology[topology.polymer.eq(family)]
        fig,axes=plt.subplots(3,len(config.schemes),figsize=config.figure_size,squeeze=False)
        handles=[]
        for k,scheme in enumerate(config.schemes):
            for row,(metric,label) in enumerate([
                ('Euler Characteristic','Boundary Euler\ncharacteristic (logged)'),
                ('Integrated Gaussian Curvature','Integrated Gaussian\ncurvature\n(dimensionless)'),
                ('Gauss-Bonnet Error','Logged Gauss–Bonnet\nerror (dimensionless)')]):
                ax=axes[row,k];found=False
                for j,(model,group) in enumerate(subset[subset.scheme.eq(scheme)].groupby('model')):
                    group=group.sort_values('frame_number');y=pd.to_numeric(group[metric],errors='coerce')
                    if not y.notna().any():continue
                    found=True
                    line,=ax.plot(group.frame_number,y,color=dict(config.colors)[scheme],
                        marker=markers[j%7],markersize=3,linestyle=['-','--',':','-.'][j%4],label=model)
                    if row==0 and k==0:handles.append(line)
                if not found:empty_panel(ax,'Not available / assumptions unverified')
                if row==0:ax.set_title(SCHEME_LABELS[scheme],color=dict(config.colors)[scheme])
                if k==0:ax.set_ylabel(label)
                if row==2:ax.set_xlabel('Frame (descriptive snapshots)')
        if handles:fig.legend(handles,[h.get_label() for h in handles],loc='lower center',ncol=min(7,len(handles)))
        fig.suptitle(f'{family}: logged topology diagnostics\nInternal water occupancy unavailable; Euler characteristic is not a cavity count')
        save_figure(fig,config.output_root/family/'F5D/topology_diagnostics',config)
        flags=['Boundary Complete','Boundary Closed','Boundary Manifold','Boundary Orientable']
        flag_rows=[]
        for (model,scheme),group in subset.groupby(['model','scheme']):
            for flag in flags:
                values=group[flag]
                flag_rows.append(dict(model=model,scheme=scheme,flag=flag,n_frames=group.frame.nunique(),
                    n_available=values.notna().sum(),fraction_true=values.dropna().astype(bool).mean() if values.notna().any() else np.nan))
        flags_table=pd.DataFrame(flag_rows)
        write_csv(flags_table,config.output_root/family/'F5D/boundary_flags.csv')
        fig,ax=plt.subplots(figsize=config.figure_size)
        if flags_table.empty:empty_panel(ax)
        else:
            matrix=flags_table.pivot(index=['model','scheme'],columns='flag',values='fraction_true').reindex(columns=flags)
            cmap=plt.get_cmap('viridis').copy();cmap.set_bad('#dddddd')
            im=ax.imshow(matrix.to_numpy(dtype=float),aspect='auto',vmin=0,vmax=1,cmap=cmap)
            ax.set_xticks(range(4),[s.removeprefix('Boundary ') for s in flags])
            ax.set_yticks(range(len(matrix)),[f'{model} / {scheme.upper()}' for model,scheme in matrix.index],fontsize=8)
            fig.colorbar(im,ax=ax,label='Fraction of available frames flagged true')
        fig.suptitle(f'{family}: boundary diagnostic flags (gray = unavailable)')
        save_figure(fig,config.output_root/family/'F5D/boundary_flags',config)


if __name__=='__main__':
    from .F5_run_all import main
    main(default_figures=['D'])
