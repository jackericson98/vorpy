"""Binned hydration geometry by chemical class; no hydrogen-bond inference."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if not __package__:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[5]))
    __package__ = 'vorpy.src.analyze.Part_II_Molecular_Analysis.F5'

from .F5_config import KEY, SCHEME_LABELS
from .F5_common import frame_summary, write_csv, save_figure, empty_panel


CLASS_COLORS={'carbon':'#555555','PEG ether oxygen':'#e69f00','PEG terminal oxygen':'#56b4e9',
              'oxygen':'#aa4499','sulfur':'#ccbb44','nitrogen':'#332288','hydrogen':'#999999','other':'#009988'}


def run(data,config):
    atoms=data['atoms'].copy()
    atoms['fraction_bin']=pd.cut(atoms.water_fraction,np.linspace(0,1,11),include_lowest=True,
        labels=(np.arange(10)+.5)/10).astype(float)
    metrics=['water_coordination','area_per_water','water_area']
    frames=frame_summary(atoms,KEY+['scheme','atom_class','fraction_bin'],metrics)
    write_csv(frames,config.output_root/'summaries/F5F_coordination_frame_bins.csv')
    for family in config.polymers:
        fig,axes=plt.subplots(2,len(config.schemes),figsize=config.figure_size,squeeze=False)
        handles=[]
        for k,scheme in enumerate(config.schemes):
            for row,metric in enumerate(['water_coordination','area_per_water']):
                ax=axes[row,k];found=False
                subset=frames[frames.polymer.eq(family)&frames.scheme.eq(scheme)]
                for atom_class,group in subset.groupby('atom_class'):
                    group=group.dropna(subset=['fraction_bin',metric])
                    if group.empty:continue
                    found=True
                    models=group.groupby(['model','fraction_bin'])[metric].mean().reset_index()
                    summary=models.groupby('fraction_bin')[metric].mean()
                    line,=ax.plot(100*summary.index.to_numpy(),summary,color=CLASS_COLORS.get(atom_class,'black'),
                        marker='o',ms=3,label=atom_class)
                    if row==0 and k==0:handles.append(line)
                if not found:empty_panel(ax)
                ax.set_xlabel('Water surface-area fraction (%)')
                if k==0:ax.set_ylabel('Unique water molecules' if row==0 else 'Water area / water molecule (Å²)')
                if row==0:ax.set_title(SCHEME_LABELS[scheme],color=dict(config.colors)[scheme])
        if handles:fig.legend(handles,[h.get_label() for h in handles],loc='lower center',ncol=min(4,len(handles)))
        fig.suptitle(f'{family}: coordination and interface geometry\nBinned model/frame means, with equal weighting of models')
        save_figure(fig,config.output_root/family/'F5F/water_coordination',config)
        selected=atoms[atoms.polymer.eq(family)]
        maximum=selected.water_area.max()
        if pd.isna(maximum) or maximum<=0:continue
        selected=selected.copy()
        edges=np.linspace(0,maximum,13)
        selected['area_bin']=pd.cut(selected.water_area,edges,include_lowest=True,labels=(edges[:-1]+edges[1:])/2).astype(float)
        area_frames=frame_summary(selected,KEY+['scheme','atom_class','area_bin'],['water_coordination'])
        write_csv(area_frames,config.output_root/family/'F5F/water_area_coordination_frame_bins.csv')
        fig,axes=plt.subplots(1,len(config.schemes),figsize=config.figure_size,squeeze=False)
        for k,scheme in enumerate(config.schemes):
            ax=axes[0,k];found=False
            for atom_class,group in area_frames[area_frames.scheme.eq(scheme)].groupby('atom_class'):
                models=group.groupby(['model','area_bin']).water_coordination.mean().reset_index()
                values=models.groupby('area_bin').water_coordination.mean()
                if values.empty:continue
                found=True
                ax.plot(values.index,values,color=CLASS_COLORS.get(atom_class,'black'),marker='o',ms=3,label=atom_class)
            if not found:empty_panel(ax)
            ax.set_xlabel('Water interface area (Å²)');ax.set_ylabel('Unique water molecules')
            ax.set_title(SCHEME_LABELS[scheme],color=dict(config.colors)[scheme])
        if handles:fig.legend(handles,[h.get_label() for h in handles],loc='lower center',ncol=min(4,len(handles)))
        fig.suptitle(f'{family}: water area versus coordination (binned model/frame means)')
        save_figure(fig,config.output_root/family/'F5F/coordination_vs_water_area',config)


if __name__=='__main__':
    from .F5_run_all import main
    main(default_figures=['F'])
