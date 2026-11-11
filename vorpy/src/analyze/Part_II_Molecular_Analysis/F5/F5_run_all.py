"""Run F5 analysis directly or as a module (see README)."""
import argparse
from dataclasses import asdict
import importlib
import json
from pathlib import Path
import pickle
import time
import traceback

import numpy as np
import pandas as pd

if not __package__:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[5]))
    __package__ = 'vorpy.src.analyze.Part_II_Molecular_Analysis.F5'

from .F5_config import Config, KEY
from .F5_common import (read_pdb, load_scheme, source_signature, independent_area_check, write_csv,
    attach_aw_exposure, match_schemes, setup_plots, exposure)


MODEL_NAMES={'PEG':('EG','PEG300','PEG1000','PEG2000','PEG3000','PEG4000','PEG6000'),
             'Heparin':('dp2','dp12','dp18','dp24','dp30')}
FIGURES={'A':'F5A_exposed_vs_buried','B':'F5B_water_burial'}


def add_residue_pairs(data,config):
    """Compare original PDB residues in addition to the inferred PEG repeat units."""
    residues=data['residues'].copy()
    keys=KEY+['residue_id']
    ref=residues[residues.scheme.eq('aw')].copy()
    ref['aw_water_fraction']=ref.water_fraction.where(ref.primary_valid)
    ref['environment']=exposure(ref.aw_water_fraction,config)
    residues=residues.merge(ref[keys+['aw_water_fraction','environment']],on=keys,how='left',validate='many_to_one')
    data['residues']=residues
    # The shared group matcher accepts unit_id as its generic aggregation key.
    pairs,report=match_schemes(residues.assign(unit_id=residues.residue_id),config,level='unit')
    pairs=pairs.rename(columns={'unit_id':'residue_id'})
    pairs['level']='residue'
    report['level']='residue'
    data['residue_pairs']=pairs
    return report


def frame_selection(tokens):
    frames=[]
    for token in tokens:
        for value in token.split(','):
            value=value.removeprefix('frame_')
            bounds=value.split('-')
            if len(bounds)==2:
                low,high=map(int,bounds)
                if low>high:raise ValueError('Reversed frame range')
                frames.extend(range(low,high+1))
            else:frames.append(int(value))
    if any(n<1 for n in frames):raise ValueError('Frame numbers must be positive')
    return tuple(f'frame_{n:04d}' for n in sorted(set(frames)))


def discover(config):
    """Include expected missing frames in the inventory, alongside discovered frames."""
    records=[]
    known={m for values in MODEL_NAMES.values() for m in values}
    if config.models and not set(config.models).issubset(known):
        raise ValueError('Unknown requested models: '+', '.join(set(config.models)-known))
    for family in config.polymers:
        for model in MODEL_NAMES[family]:
            if config.models and model not in config.models:continue
            folder=config.dataset_root/family/model
            discovered={p.name for p in folder.glob('frame_*') if p.is_dir()}
            frames=config.frames or tuple(sorted(discovered | {f'frame_{i:04d}' for i in range(1,12)}))
            for frame in frames:
                path=folder/frame
                for scheme in config.schemes:
                    log,pdb=path/f'{scheme}_logs.csv',path/'frame.pdb'
                    missing=[str(p) for p in [log,pdb] if not p.is_file() or p.stat().st_size==0]
                    records.append(dict(polymer=family,model=model,frame=frame,scheme=scheme,
                        log=log,pdb=pdb,status='missing' if missing else 'available',missing_files='; '.join(missing)))
    return pd.DataFrame(records)


def run(config,figures):
    started=time.perf_counter()
    config.dataset_root=Path(config.dataset_root).resolve()
    config.output_root=Path(config.output_root).resolve()
    if config.output_root==config.dataset_root or config.dataset_root in config.output_root.parents:
        raise ValueError('Output root must be outside the scientific source dataset')
    if not 0<=config.buried_threshold<config.exposed_threshold<=1:
        raise ValueError('Exposure thresholds must satisfy 0 <= buried < exposed <= 1')
    if config.area_atol<0 or config.area_rtol<0 or config.minimum_interface_area<0:
        raise ValueError('Area tolerances and minimum interface area must be nonnegative')
    out=config.output_root
    (out/'validation').mkdir(parents=True,exist_ok=True)
    (out/'summaries').mkdir(parents=True,exist_ok=True)
    inventory=discover(config)
    if inventory.empty:raise ValueError('No model/family combinations selected')
    write_csv(inventory,out/'validation/file_inventory.csv')
    (out/'validation/run_config.json').write_text(json.dumps(asdict(config),indent=2,default=str),encoding='utf-8')
    tables={name:[] for name in ['atoms','residues','units','topology','validation']}
    failures=[];checks=[];shards=[];interface_totals=[];pdb_reports=[]
    for metadata,rows in inventory.groupby(KEY,sort=False):
        family,model,frame=metadata
        pdb=None
        print(f'Loading {family}/{model}/{frame}',flush=True)
        for row in rows.itertuples():
            md=dict(polymer=family,model=model,frame=frame,scheme=row.scheme)
            if row.status=='missing':
                failures.append(dict(md,status='missing',reason=row.missing_files));continue
            try:
                signature=source_signature([row.log,row.pdb],config)
                cache=out/'cache'/family/model/frame/f'{row.scheme}.pkl'
                bundle=None
                if config.cache and cache.is_file():
                    try:
                        with cache.open('rb') as stream:stored=pickle.load(stream)
                        if stored['signature']==signature:bundle=stored['tables']
                    except (OSError,EOFError,pickle.UnpicklingError,KeyError):
                        bundle=None
                if bundle is None:
                    if pdb is None:
                        pdb=read_pdb(row.pdb,family,config)
                    bundle=load_scheme(row.log,pdb,config,md)
                    if model in {'PEG3000','dp24'} and frame=='frame_0001':
                        check=independent_area_check(row.log,bundle['atoms'],pdb,config,md)
                        if not check.passed.all():raise ValueError('Independent raw water-area check failed')
                        bundle['area_checks']=check
                    if config.cache:
                        cache.parent.mkdir(parents=True,exist_ok=True)
                        temporary=cache.with_suffix('.tmp')
                        with temporary.open('wb') as stream:pickle.dump({'signature':signature,'tables':bundle},stream,pickle.HIGHEST_PROTOCOL)
                        temporary.replace(cache)
                if 'area_checks' in bundle:checks.append(bundle['area_checks'])
                pdb_reports.append(bundle['pdb_inventory'])
                for name in tables:tables[name].append(bundle[name])
                surface_path=out/'summaries/interfaces'/family/model/frame/f'{row.scheme}.csv'
                write_csv(bundle['interfaces'],surface_path)
                shards.append(dict(md,path=str(surface_path.relative_to(out)),n_interfaces=len(bundle['interfaces']),signature=signature))
                for kind,group in bundle['interfaces'].groupby('interface_type'):
                    interface_totals.append(dict(md,interface_type=kind,n_interfaces=len(group),
                        surface_area=group.surface_area.sum(min_count=len(group)),contact_area=group.contact_area.sum(min_count=len(group))))
                print(f'  {row.scheme}: {len(bundle["atoms"])} atoms, {len(bundle["interfaces"])} interfaces, '
                    f'{int(bundle["atoms"].primary_valid.sum())} valid atoms',flush=True)
            except Exception as error:
                failures.append(dict(md,status='failed',reason=f'{type(error).__name__}: {error}'))
                print(f'  {row.scheme}: FAILED {error}',flush=True)
        # Release the full PDB before loading the next frame (up to ~200,000 atoms).
        del pdb
    data={name:pd.concat(values,ignore_index=True) for name,values in tables.items() if values}
    failure_table=pd.DataFrame(failures,columns=KEY+['scheme','status','reason'])
    write_csv(failure_table,out/'validation/failures.csv')
    write_csv(pd.concat(checks,ignore_index=True) if checks else pd.DataFrame(columns=KEY+['scheme','atom_index','passed']),out/'validation/independent_area_checks.csv')
    write_csv(pd.concat(pdb_reports,ignore_index=True).drop_duplicates() if pdb_reports else pd.DataFrame(),out/'validation/pdb_residue_inventory.csv')
    write_csv(pd.DataFrame(shards),out/'summaries/interface_shards.csv')
    write_csv(pd.DataFrame(interface_totals),out/'summaries/interface_totals.csv')
    if not tables['atoms']:
        result={'status':'no_valid_logs','requested_logs':len(inventory),'failed_or_missing':len(failures)}
        (out/'validation/run_report.json').write_text(json.dumps(result,indent=2))
        return result
    data['atoms'],data['units']=attach_aw_exposure(data['atoms'],data['units'],config)
    # Scheme comparisons use raw identity-verified pairs; invalid rows remain in exported tables.
    data['atom_pairs'],atom_report=match_schemes(data['atoms'],config)
    data['unit_pairs'],unit_report=match_schemes(data['units'],config,level='unit')
    residue_report=add_residue_pairs(data,config)
    write_csv(pd.concat([atom_report,unit_report,residue_report],ignore_index=True),out/'validation/scheme_matching.csv')
    for name,table in data.items():
        write_csv(table,out/('validation' if name=='validation' else 'summaries')/f'{name}.csv')
    setup_plots(config)
    figure_failures=[]
    for letter in figures:
        try:
            print(f'Plotting F5{letter}',flush=True)
            importlib.import_module('.'+FIGURES[letter],__package__).main(config)
        except Exception as error:
            figure_failures.append({'figure':letter,'reason':f'{type(error).__name__}: {error}','traceback':traceback.format_exc()})
            print(f'F5{letter} FAILED: {error}',flush=True)
    write_csv(pd.DataFrame(figure_failures,columns=['figure','reason','traceback']),out/'validation/figure_failures.csv')
    result={'status':'complete' if not failures and not figure_failures else 'completed_with_exclusions',
        'requested_logs':len(inventory),'loaded_logs':len(data['validation']),'failed_or_missing':len(failures),
        'models_loaded':data['atoms'][['polymer','model']].drop_duplicates().shape[0],
        'model_frames_loaded':data['atoms'][KEY].drop_duplicates().shape[0],
        'atom_observations':len(data['atoms']),'invalid_primary_atoms':int((~data['atoms'].primary_valid).sum()),
        'repeat_observations':len(data['units']),'invalid_primary_units':int((~data['units'].primary_valid).sum()),
        'independent_area_checks':sum(len(check) for check in checks),
        'independent_area_checks_passed':all(check.passed.all() for check in checks) if checks else None,
        'matched_atom_pairs':int(atom_report.n_valid.sum()) if not atom_report.empty else 0,
        'matched_residue_pairs':int(residue_report.n_valid.sum()) if not residue_report.empty else 0,
        'identity_mismatches':int(atom_report.n_identity_mismatch.sum()) if not atom_report.empty else 0,
        'figures_requested':figures,'figure_failures':figure_failures,'elapsed_seconds':time.perf_counter()-started,
        'internal_water_occupancy':'unavailable','statistics':'Descriptive frame variability; no independent molecular-replicate claims'}
    (out/'validation/run_report.json').write_text(json.dumps(result,indent=2,default=str),encoding='utf-8')
    report=['# F5 validation report','',json.dumps(result,indent=2,default=str),'',
        'Validated: raw/reader row counts; zero-based endpoint-to-PDB mapping; names, residue IDs and coordinates; '
        'water-site identifiers; undirected interface deduplication; atomic area partitions; cross-scheme identities. '
        'Independent raw area checks are sampled, not a revalidation of the VorPy solver.',
        '', 'Exploratory: atomic exposure, water polymer-contact categories, and coordination/interface associations. '
        'Frames are correlated snapshots. Error bars/bands describe variability, not confidence intervals.',
        '', 'Unavailable: internal/external water occupancy, cavity counts, solvent channels, penetration or diffusion. '
        'No validated envelope/segmentation is supplied. Topology flags are retained; non-orientable or incomplete '
        'boundaries do not support Gauss–Bonnet/cavity interpretations.',
        '', 'See failures.csv, validation.csv, scheme_matching.csv, independent_area_checks.csv, and F5B_coverage.csv '
        'for model/frame exclusions, missing fields, and boundary limitations.']
    (out/'validation/validation_report.md').write_text('\n'.join(report),encoding='utf-8')
    print(json.dumps(result,indent=2,default=str),flush=True)
    return result



def main(default_figures=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset-root',type=Path)
    parser.add_argument('--output-root',type=Path)
    parser.add_argument('--polymer',nargs='+',choices=['PEG','Heparin'])
    parser.add_argument('--models',nargs='+')
    parser.add_argument('--frames',nargs='+',help='1 2 3, 1-11, or frame_0001')
    parser.add_argument('--schemes',nargs='+',choices=['aw','pow','prm'])
    parser.add_argument('--figures',nargs='*',choices=list(FIGURES),default=default_figures or list(FIGURES),
                        help='Panels A and B; pass --figures alone to prepare shared tables only')
    parser.add_argument('--include-incomplete',action='store_true')
    parser.add_argument('--exclude-hydrogens',action='store_true')
    parser.add_argument('--no-cache',action='store_true')
    parser.add_argument('--show',action='store_true')
    parser.add_argument('--dpi',type=int)
    parser.add_argument('--width',type=float)
    parser.add_argument('--height',type=float)
    parser.add_argument('--area-rtol',type=float)
    parser.add_argument('--area-atol',type=float)
    args=parser.parse_args()
    config=Config()
    for argument,field in [('dataset_root','dataset_root'),('output_root','output_root'),('polymer','polymers'),
                           ('models','models'),('schemes','schemes'),('dpi','dpi'),('area_rtol','area_rtol'),('area_atol','area_atol')]:
        value=getattr(args,argument)
        if value is not None:setattr(config,field,tuple(value) if isinstance(value,list) else value)
    if args.frames:
        try:config.frames=frame_selection(args.frames)
        except ValueError as error:parser.error(str(error))
    config.figure_size=(args.width or config.figure_size[0],args.height or config.figure_size[1])
    config.include_incomplete=args.include_incomplete
    config.include_hydrogens=not args.exclude_hydrogens
    config.cache=not args.no_cache
    config.show=args.show
    result=run(config,args.figures)
    if result.get('figure_failures') or result['status']=='no_valid_logs':raise SystemExit(1)


if __name__=='__main__':
    main()
