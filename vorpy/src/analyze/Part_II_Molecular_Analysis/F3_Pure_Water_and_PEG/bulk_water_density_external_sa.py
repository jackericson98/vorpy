#!/usr/bin/env python3
"""VorPy bulk-water KDE: Tkinter-select AW, Power, Primitive logs; no PDB.

Requires numpy, pandas, matplotlib, scipy; tkinter for interactive selection.
Run: python bulk_water_density.py
Or: python bulk_water_density.py --aw aw.csv --power power.csv --primitive primitive.csv --out figures
Edit CONFIG to customize all plot settings. Only three density plots are produced. Surface area uses inter-molecular interfaces from the Surfaces section, excluding all intra-water interfaces.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import ast
import csv

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

def read_external_surfaces(path, atoms, cfg):
    """Sum Surface Area once per interface for each water it bounds.

    Intramolecular interfaces are excluded. A water-water interface contributes
    once to EACH water's external area, not twice to one water.
    """
    idx_to_mol = {int(r.Index): r.molecule for r in atoms.itertuples(index=False)
                  if pd.notna(r.Index)}
    water_mols = set(atoms.loc[atoms['Residue'].isin(cfg['water_residues']), 'molecule'])
    area_by_mol = {m: 0.0 for m in water_mols}
    face_counts = {m: 0 for m in water_mols}
    seen = set()
    with open(path, encoding='utf-8-sig', newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            if row and row[0].strip() == 'Surfaces':
                header = next(reader)
                break
        else:
            raise ValueError(f'{path}: missing Surfaces section; cannot calculate molecular exterior area')
        required = {'Ball 1', 'Ball 2', 'Surface Area'}
        if not required.issubset(header):
            raise ValueError(f'{path}: Surfaces missing {sorted(required-set(header))}')
        columns = {name: header.index(name) for name in required}
        n_records = 0
        for row in reader:
            try:
                a = int(row[columns['Ball 1']]); b = int(row[columns['Ball 2']])
                area = float(row[columns['Surface Area']])
            except (IndexError, ValueError, TypeError):
                if n_records: break
                continue
            n_records += 1
            if not np.isfinite(area) or area < 0:
                raise ValueError(f'{path}: invalid surface area for {a}, {b}')
            key = (min(a, b), max(a, b))
            if key in seen:
                raise ValueError(f'{path}: duplicate surface pair {key}; cannot safely sum')
            seen.add(key)
            ma = idx_to_mol.get(a); mb = idx_to_mol.get(b)
            # Logs can omit atomic rows for some cells; a missing endpoint
            # is external to any fully represented three-atom water molecule.
            if ma is not None and ma == mb:
                continue
            for mol in (ma, mb):
                if mol in water_mols:
                    area_by_mol[mol] += area
                    face_counts[mol] += 1
    if not n_records:
        raise ValueError(f'{path}: Surfaces section contains no usable interfaces')
    return area_by_mol, face_counts

def molecular_table(atoms, scheme, cfg, area_by_mol, face_counts):
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
                       'volume':float(g['Volume'].sum()),'surface_area':float(area_by_mol.get(mol, 0.0)),
                       'external_face_count':face_counts.get(mol, 0),
                       'contacts':len(neighbors),'unresolved_neighbor_indices':len(unresolved),
                       'oxygen_boundary_distance':oxygen_margin,'cell_envelope_margin':cell_margin,
                       'complete':bool(g['Complete Cell?'].all())})
    return pd.DataFrame(result),lo,hi



CONFIG = {
    'primary_buffer_A': 5.0,
    'filter_mode': 'cell_envelope', # 'cell_envelope', 'oxygen', or 'none'
    'require_complete': True,
    'require_resolved_contacts': True,
    'water_residues': {'SOL','HOH','WAT','TIP3','TIP3P'},
    'oxygen_names': {'OW','O','OH2'},
    'hydrogen_names': {'HW1','HW2','H1','H2'},
    'metrics': {'volume':'Molecular volume (Å³)',
                'surface_area':'External molecular surface area (Å²)',
                'contacts':'Unique neighboring molecules'},
    'colors': {'AW':'#2676bb','Power':'#cf4848','Primitive':'#32965a'},
    'figure_size': (8.2,5.4), 'dpi':250, 'font_size':12,
    'line_width':2.6, 'mean_line_width':1.5, 'mean_line_style':'--',
    'fill_alpha':0.12, 'show_fill':True, 'show_grid':True,
    'show_title':False, 'show_mean_lines':True,
    'show_rug':False, 'rug_alpha':0.12, 'rug_height_fraction':0.035,
    'show_histogram':False, 'hist_bins':'fd', 'hist_alpha':0.08,
    'kde_bandwidth':1.2, # multiplier on Scott's bandwidth: >1 smoother, <1 more detail
    'kde_points':600, 'kde_padding_bandwidths':3.5,
    'x_limits':{}, 'y_limits':{}, # e.g. {'volume':(15,55)}
    'formats':['png','svg'],
    'legend_font_size':10, 'legend_location':'best',
    'font_family':'DejaVu Sans',
}

def density_plot(tables,metric,label,out,cfg):
    arrays={s:t.loc[t.selected,metric].dropna().to_numpy(float) for s,t in tables.items()}
    arrays={s:a for s,a in arrays.items() if len(a)}
    if not arrays:return
    all_values=np.concatenate(list(arrays.values()))
    fig,ax=plt.subplots(figsize=cfg['figure_size'])
    span=max(float(np.ptp(all_values)),1.0)
    left=float(all_values.min()-0.08*span);right=float(all_values.max()+0.08*span)
    for scheme,a in arrays.items():
        color=cfg['colors'][scheme]
        mu=float(np.mean(a));sd=float(np.std(a,ddof=1)) if len(a)>1 else float('nan')
        legend=f'{scheme}: {mu:.2f} ± {sd:.2f} (n={len(a)})'
        if len(a)>1 and np.std(a)>0:
            kde=gaussian_kde(a,bw_method=lambda obj:obj.scotts_factor()*cfg['kde_bandwidth'])
            bandwidth=float(np.sqrt(kde.covariance[0,0]))
            x=np.linspace(min(left,float(a.min()-cfg['kde_padding_bandwidths']*bandwidth)),
                          max(right,float(a.max()+cfg['kde_padding_bandwidths']*bandwidth)),cfg['kde_points'])
            y=kde(x)
            ax.plot(x,y,color=color,lw=cfg['line_width'],label=legend)
            if cfg['show_fill']:ax.fill_between(x,y,color=color,alpha=cfg['fill_alpha'])
        else:
            ax.axvline(mu,color=color,lw=cfg['line_width'],label=legend)
        if cfg['show_mean_lines'] and len(a)>1:
            ax.axvline(mu,color=color,lw=cfg['mean_line_width'],ls=cfg['mean_line_style'])
        if cfg['show_histogram']:
            ax.hist(a,bins=cfg['hist_bins'],density=True,color=color,alpha=cfg['hist_alpha'])
        if cfg['show_rug']:
            ax.plot(a,np.zeros_like(a), '|',color=color,alpha=cfg['rug_alpha'],markersize=5)
    ax.set_xlabel(label);ax.set_ylabel('Probability density')
    if cfg['show_title']:ax.set_title(f'Bulk water: {label}')
    if cfg['show_grid']:ax.grid(axis='y',alpha=.18)
    ax.spines[['top','right']].set_visible(False)
    ax.set_ylim(bottom=0)
    if metric in cfg['x_limits']:ax.set_xlim(*cfg['x_limits'][metric])
    if metric in cfg['y_limits']:ax.set_ylim(*cfg['y_limits'][metric])
    ax.legend(loc=cfg['legend_location'],fontsize=cfg['legend_font_size'],frameon=False)
    fig.tight_layout()
    for ext in cfg['formats']:fig.savefig(out/f'{metric}_density.{ext}',dpi=cfg['dpi'],bbox_inches='tight')
    plt.close(fig)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('aw','power','primitive','out'):p.add_argument('--'+name)
    p.add_argument('--buffer',type=float)
    p.add_argument('--bandwidth',type=float)
    args=p.parse_args()
    cfg=CONFIG.copy()
    if args.buffer is not None:cfg['primary_buffer_A']=args.buffer
    if args.bandwidth is not None:cfg['kde_bandwidth']=args.bandwidth
    if cfg['kde_bandwidth']<=0:raise ValueError('KDE bandwidth must be positive')
    paths={'AW':args.aw,'Power':args.power,'Primitive':args.primitive}
    if not any(paths.values()):
        paths={s:select_file(f'Select {s} logs (Cancel to skip)') for s in paths}
    paths={s:path for s,path in paths.items() if path}
    if not paths:raise SystemExit('No logs selected')
    out=Path(args.out or select_directory() or 'bulk_water_density_output');out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':cfg['font_size'],'font.family':cfg['font_family']})
    tables={};records=[]
    for scheme,path in paths.items():
        atoms=read_atoms(path)
        areas,faces=read_external_surfaces(path,atoms,cfg)
        tab,_,_=molecular_table(atoms,scheme,cfg,areas,faces)
        if tab.empty:raise ValueError(f'No eligible water molecules for {scheme}')
        if (tab['external_face_count']==0).any():
            raise ValueError(f'{scheme}: some eligible waters have no external surface records')
        if cfg['filter_mode']=='cell_envelope':margin=tab['cell_envelope_margin']
        elif cfg['filter_mode']=='oxygen':margin=tab['oxygen_boundary_distance']
        elif cfg['filter_mode']=='none':margin=pd.Series(np.inf,index=tab.index)
        else:raise ValueError('Invalid filter_mode')
        tab['selected']=margin>=cfg['primary_buffer_A']
        tables[scheme]=tab
        for metric in cfg['metrics']:
            values=tab.loc[tab.selected,metric].dropna()
            records.append({'scheme':scheme,'metric':metric,'n':len(values),
                            'mean':values.mean(),'sd':values.std(ddof=1)})
        print(f'{scheme}: {len(tab)} eligible, {int(tab.selected.sum())} selected')
    for metric,label in cfg['metrics'].items():density_plot(tables,metric,label,out,cfg)
    pd.DataFrame(records).to_csv(out/'density_statistics.csv',index=False)
    print('Plots and statistics:',out.resolve())
    print('NOTE: density curves are KDE estimates; reported mean and SD use unsmoothed molecular values.')
    print('NOTE: geometric bounding box is a proxy, not a validated periodic boundary.')

if __name__=='__main__':main()
