#!/usr/bin/env python3
"""VorPy bulk-water AW/Power/Primitive comparison, log-only, Tkinter selectors.

Requires: numpy, pandas, matplotlib (Tkinter normally ships with Python).
Run: python bulk_water_compare.py
Optional: python bulk_water_compare.py --aw AW.csv --power POW.csv --primitive PRM.csv --out results

Edit CONFIG to customize plots and filtering. No property-based outlier rejection.
"""
import argparse
import ast
import csv
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

CONFIG = {
    'primary_buffer_A': 5.0,
    'buffers_A': [0, 2.5, 5, 7.5, 10],
    'filter_mode': 'cell_envelope',  # 'cell_envelope' or 'oxygen'
    'boundary_reference': 'all_atom_bbox', # coordinates in log, NOT true simulation box
    'require_complete': True,
    'require_resolved_contacts': True,
    'water_residues': {'SOL','HOH','WAT','TIP3','TIP3P'},
    'oxygen_names': {'OW','O','OH2'},
    'hydrogen_names': {'HW1','HW2','H1','H2'},
    'figure_size': (9, 6), 'dpi': 200, 'font_size': 12,
    'marker_size': 12, 'alpha': 0.45, 'hist_bins': 45,
    'distance_bin_width_A': 1.5, 'min_bin_n': 3,
    'show_sd': True, 'show_grid': True,
    'formats': ['png','svg'],
    'colors': {'AW':'#2874b2','Power':'#c44747','Primitive':'#348a51'},
    'metrics': {'volume': 'Volume (Å³)', 'surface_area': 'Summed atomic surface area (Å²)',
                'contacts': 'Unique neighboring molecules'},
    'x_limits': {}, 'y_limits': {},  # plot-name -> (min,max), e.g. 'volume_distribution':(0,70)
    'percent_differences': True,
}

def select_file(title):
    import tkinter as tk
    from tkinter import filedialog
    root=tk.Tk();root.withdraw()
    try: return filedialog.askopenfilename(title=title,filetypes=[('VorPy CSV','*.csv'),('All files','*.*')])
    finally:root.destroy()

def select_directory():
    import tkinter as tk
    from tkinter import filedialog
    root=tk.Tk();root.withdraw()
    try:return filedialog.askdirectory(title='Choose output folder')
    finally:root.destroy()

def parse_neighbors(s):
    try:
        x=ast.literal_eval(s)
        return [int(v) for v in x] if isinstance(x,list) else []
    except (ValueError,TypeError,SyntaxError):return []

def read_atoms(path):
    with open(path,encoding='utf-8-sig',newline='') as f:
        reader=csv.reader(f)
        for row in reader:
            if row and row[0].strip()=='Atoms':
                header=next(reader);break
        else:raise ValueError(f'No Atoms section in {path}')
        data=[]
        for row in reader:
            if row and row[0].strip().isdigit():data.append((row+['']*len(header))[:len(header)])
            elif data:break
    df=pd.DataFrame(data,columns=header)
    required={'Index','Name','Chain','Residue','Residue Sequence','X','Y','Z','Volume',
              'Surface Area','Complete Cell?','Maximum Point Distance','neighbors'}
    missing=required-set(df.columns)
    if missing:raise ValueError(f'{path}: missing columns: {sorted(missing)}')
    for c in ['Index','Residue Sequence','X','Y','Z','Volume','Surface Area','Maximum Point Distance']:
        df[c]=pd.to_numeric(df[c],errors='coerce')
    df['Complete Cell?']=df['Complete Cell?'].astype(str).str.strip().str.lower().eq('true')
    for c in ['Name','Chain','Residue']:df[c]=df[c].fillna('').astype(str).str.strip()
    df['neighbors_parsed']=df['neighbors'].map(parse_neighbors)
    # Stable identity: chain + residue type + residue sequence, not atom index.
    df['molecule']=list(zip(df['Chain'],df['Residue'],df['Residue Sequence']))
    return df

def molecular_table(atoms, scheme, cfg):
    # Domain bounds use ALL logged atoms to avoid shrinking the boundary after filtering.
    xyz=atoms[['X','Y','Z']].to_numpy(float)
    if not np.isfinite(xyz).all():raise ValueError(f'{scheme}: missing atomic coordinates')
    lo=xyz.min(axis=0);hi=xyz.max(axis=0)
    waters=atoms[atoms['Residue'].isin(cfg['water_residues'])].copy()
    idx_to_mol={int(r.Index):r.molecule for r in atoms.itertuples(index=False) if pd.notna(r.Index)}
    result=[]
    for mol,g in waters.groupby('molecule',sort=False):
        oxy=g[g['Name'].isin(cfg['oxygen_names'])]
        hyd=g[g['Name'].isin(cfg['hydrogen_names'])]
        if len(g)!=3 or len(oxy)!=1 or len(hyd)!=2:continue
        if cfg['require_complete'] and not g['Complete Cell?'].all():continue
        if not np.isfinite(g[['Volume','Surface Area','Maximum Point Distance']].to_numpy(float)).all():continue
        neighbors=set();unresolved=set()
        for a in g.itertuples(index=False):
            for neighbor_index in a.neighbors_parsed:
                other=idx_to_mol.get(neighbor_index)
                if other is None:unresolved.add(neighbor_index)
                elif other!=mol:neighbors.add(other)
        if cfg['require_resolved_contacts'] and unresolved:continue
        center=oxy[['X','Y','Z']].iloc[0].to_numpy(float)
        coords=g[['X','Y','Z']].to_numpy(float)
        radii=g['Maximum Point Distance'].to_numpy(float)
        oxygen_margin=float(np.min(np.r_[center-lo,hi-center]))
        # Conservative radial envelope; logged point-distance interpretation should be verified.
        cell_margin=float(np.min(np.column_stack([coords-lo,hi-coords])-radii[:,None]))
        result.append({'molecule':repr(mol),'scheme':scheme,'x':center[0],'y':center[1],'z':center[2],
                       'volume':float(g['Volume'].sum()),'surface_area':float(g['Surface Area'].sum()),
                       'contacts':len(neighbors),'unresolved_neighbor_indices':len(unresolved),
                       'oxygen_boundary_distance':oxygen_margin,'cell_envelope_margin':cell_margin,
                       'complete':bool(g['Complete Cell?'].all())})
    return pd.DataFrame(result),lo,hi

def save(fig,out,name,cfg):
    fig.tight_layout()
    for fmt in cfg['formats']:fig.savefig(out/f'{name}.{fmt}',dpi=cfg['dpi'],bbox_inches='tight')
    plt.close(fig)

def decorate(ax,name,cfg):
    if name in cfg['x_limits']:ax.set_xlim(*cfg['x_limits'][name])
    if name in cfg['y_limits']:ax.set_ylim(*cfg['y_limits'][name])
    if cfg['show_grid']:ax.grid(alpha=.2)

def plot_all(data,matched,out,cfg):
    plt.rcParams.update({'font.size':cfg['font_size']})
    metric_names=cfg['metrics']; schemes=list(data)
    for metric,label in metric_names.items():
        name=f'{metric}_distribution';fig,ax=plt.subplots(figsize=cfg['figure_size'])
        for s in schemes:
            d=data[s];d=d[d['selected']]
            if len(d):ax.hist(d[metric],bins=cfg['hist_bins'],density=True,histtype='step',linewidth=1.8,
                              label=f'{s} (n={len(d)})',color=cfg['colors'][s])
        ax.set(xlabel=label,ylabel='Probability density',title=f'Interior-water {label.lower()}');ax.legend();decorate(ax,name,cfg);save(fig,out,name,cfg)
        name=f'{metric}_boundary';fig,ax=plt.subplots(figsize=cfg['figure_size'])
        for s in schemes:
            d=data[s];ax.scatter(d['oxygen_boundary_distance'],d[metric],s=cfg['marker_size'],alpha=cfg['alpha'],color=cfg['colors'][s],label=s)
        ax.set(xlabel='Oxygen distance to logged-coordinate bounding box (Å)',ylabel=label,title='Molecular properties versus boundary distance')
        ax.legend();decorate(ax,name,cfg);save(fig,out,name,cfg)
        records=[]
        for s in schemes:
            d=data[s].copy();width=cfg['distance_bin_width_A'];d['bin']=np.floor(d['oxygen_boundary_distance']/width)*width
            agg=d.groupby('bin')[metric].agg(['mean','std','count']).reset_index()
            agg['scheme']=s;records.append(agg)
        pd.concat(records,ignore_index=True).to_csv(out/f'{metric}_distance_bins.csv',index=False)
        name=f'{metric}_binned';fig,ax=plt.subplots(figsize=cfg['figure_size'])
        for s,agg in zip(schemes,records):
            a=agg[agg['count']>=cfg['min_bin_n']];x=a['bin']+cfg['distance_bin_width_A']/2;y=a['mean']
            ax.plot(x,y,'o-',label=s,color=cfg['colors'][s])
            if cfg['show_sd']:ax.fill_between(x,y-a['std'],y+a['std'],alpha=.12,color=cfg['colors'][s])
        ax.set(xlabel='Oxygen distance to logged-coordinate bounding box (Å)',ylabel=f'Mean {label}',title='Boundary-distance binned means ± SD')
        ax.legend();decorate(ax,name,cfg);save(fig,out,name,cfg)
        sensitivity=[]
        for s in schemes:
            for buffer in cfg['buffers_A']:
                d=data[s];a=d[d['interior_margin']>=buffer][metric]
                sensitivity.append({'scheme':s,'buffer_A':buffer,'n':len(a),'mean':a.mean(),'sd':a.std(ddof=1)})
        sens=pd.DataFrame(sensitivity);sens.to_csv(out/f'{metric}_buffer_sensitivity.csv',index=False)
        name=f'{metric}_buffer_sensitivity';fig,ax=plt.subplots(figsize=cfg['figure_size'])
        for s in schemes:
            a=sens[sens['scheme']==s];ax.plot(a['buffer_A'],a['mean'],'o-',label=s,color=cfg['colors'][s])
            if cfg['show_sd']:ax.fill_between(a['buffer_A'],a['mean']-a['sd'],a['mean']+a['sd'],alpha=.12,color=cfg['colors'][s])
        ax.set(xlabel='Minimum geometric interior margin (Å)',ylabel=f'Mean {label}',title='Buffer sensitivity (no outlier filtering)')
        ax.legend();decorate(ax,name,cfg);save(fig,out,name,cfg)
        if 'AW' not in schemes:continue
        for other in ['Power','Primitive']:
            if other not in schemes:continue
            a=matched[matched['selected_common']].copy()
            if a.empty:continue
            for suffix,scale in [('absolute',False),('percent',True)]:
                if scale and not cfg['percent_differences']:continue
                base=a[f'{metric}_AW'];val=a[f'{metric}_{other}'];valid=base.ne(0) if scale else base.notna()
                d=(val-base)[valid]
                if scale:d=100*d/base[valid]
                name=f'{metric}_{other.lower()}_minus_aw_{suffix}'
                fig,ax=plt.subplots(figsize=cfg['figure_size'])
                ax.scatter(base[valid],d,s=cfg['marker_size'],alpha=cfg['alpha'],color=cfg['colors'][other])
                ax.axhline(0,color='black',ls='--',lw=1)
                ax.set(xlabel=f'AW {label}',ylabel=f'{other} − AW '+('(%)' if scale else f'({"Å³" if metric=="volume" else "Å²" if metric=="surface_area" else "contacts"})'),
                       title=f'{other} vs AW: matched interior waters (n={len(d)})')
                decorate(ax,name,cfg);save(fig,out,name,cfg)
    # A long-form per-water table is also exported separately.

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for flag in ['aw','power','primitive','out']:p.add_argument('--'+flag)
    p.add_argument('--buffer',type=float)
    args=p.parse_args();cfg=CONFIG.copy()
    if args.buffer is not None:cfg['primary_buffer_A']=args.buffer
    paths={'AW':args.aw,'Power':args.power,'Primitive':args.primitive}
    if not any(paths.values()):
        paths['AW']=select_file('Select AW VorPy logs (Cancel to skip)')
        paths['Power']=select_file('Select Power VorPy logs (Cancel to skip)')
        paths['Primitive']=select_file('Select Primitive VorPy logs (Cancel to skip)')
    paths={k:v for k,v in paths.items() if v}
    if not paths:raise SystemExit('No logs selected.')
    out=Path(args.out or select_directory() or 'bulk_water_comparison_output');out.mkdir(parents=True,exist_ok=True)
    tables={};bounds={}
    for scheme,path in paths.items():
        tab,lo,hi=molecular_table(read_atoms(path),scheme,cfg)
        if tab.empty:raise ValueError(f'No eligible water molecules in {scheme} logs')
        tables[scheme]=tab;bounds[scheme]=(lo,hi)
    # Use identical spatial reference across schemes, not separately fitted bboxes.
    # Coordinate extrema from AW (or first available) are the reference, and a shared
    # molecule's margin is the minimum across schemes (conservative common selection).
    for s,t in tables.items():
        t['interior_margin']=t['cell_envelope_margin'] if cfg['filter_mode']=='cell_envelope' else t['oxygen_boundary_distance']
        t['selected']=t['interior_margin']>=cfg['primary_buffer_A']
    combined=pd.concat(tables.values(),ignore_index=True);combined.to_csv(out/'water_molecules_by_scheme.csv',index=False)
    merged=None
    for s,t in tables.items():
        cols=['molecule','volume','surface_area','contacts','interior_margin','selected']
        q=t[cols].rename(columns={c:f'{c}_{s}' for c in cols if c!='molecule'})
        merged=q if merged is None else merged.merge(q,on='molecule',how='outer',validate='one_to_one')
    selected_cols=[f'selected_{s}' for s in tables]
    merged['selected_common']=merged[selected_cols].fillna(False).all(axis=1)
    merged.to_csv(out/'matched_water_molecules.csv',index=False)
    plot_all(tables,merged,out,cfg)
    with open(out/'README.txt','w') as f:
        f.write('Water volume = sum of three atomic Voronoi cell volumes.\n'
                'Surface area = sum of three atomic cell surface areas, NOT exterior molecular area.\n'
                'Contacts = distinct neighboring molecules from atomic neighbor indices; intramolecular contacts excluded.\n'
                'Only complete three-atom waters with resolvable neighbor indices are used by default.\n'
                'Boundary proxy = bounding box of logged coordinates; not a physical periodic box.\n'
                'The cell envelope uses logged maximum point distance as a conservative radial bound.\n'
                'Cross-scheme differences use the intersection of waters passing each scheme filter.\n'
                'No filtering based on volume, surface area, contact count or their statistical extremes.\n'
                'WARNING: For a validated bulk baseline, reconstruct periodic geometry or establish a trustworthy physical boundary.\n')
    for s,t in tables.items():print(f'{s}: eligible {len(t)}, interior {int(t.selected.sum())}, logged coordinate bounds {bounds[s]}')
    print('Matched interior waters across all supplied schemes:',int(merged.selected_common.sum()))
    print('Results:',out.resolve())

if __name__=='__main__':main()
