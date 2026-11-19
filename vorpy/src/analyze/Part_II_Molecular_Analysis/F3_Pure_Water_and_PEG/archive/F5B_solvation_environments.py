"""AW repeat exposure defines provisional environments for every scheme."""
if not __package__:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[5]))
    __package__ = 'vorpy.src.analyze.Part_II_Molecular_Analysis.F5'

from .F5_config import KEY
from .F5_common import frame_summary, write_csv, distribution_panels


def run(data,config):
    atoms=data['atoms']
    metrics=['volume','surface_area','water_fraction','polymer_coordination','water_coordination','area_per_water']
    groups=KEY+['scheme','environment','atom_class']
    frames=frame_summary(atoms,groups,metrics)
    write_csv(frames,config.output_root/'summaries/F5B_chemical_environment_frame_means.csv')
    all_atoms=frame_summary(atoms,KEY+['scheme','environment'],metrics)
    write_csv(all_atoms,config.output_root/'summaries/F5B_environment_frame_means.csv')
    zero=data['units'][KEY+['scheme','unit_id','water_area','zero_water','aw_zero_water','primary_valid','environment']]
    write_csv(zero,config.output_root/'summaries/F5B_zero_water_units.csv')
    for family in config.polymers:
        distribution_panels(all_atoms[all_atoms.polymer.eq(family)],metrics,
            f'{family}: AW reference exposure bins',config.output_root/family/'F5B/solvation_environments',config)
        for atom_class in ['PEG ether oxygen','carbon']:
            subset=frames[frames.polymer.eq(family)&frames.atom_class.eq(atom_class)]
            if subset.n_observations.sum()<5:continue
            distribution_panels(subset,metrics,f'{family}: {atom_class}, AW reference exposure bins',
                config.output_root/family/'F5B'/atom_class.replace(' ','_'),config)


if __name__=='__main__':
    from .F5_run_all import main
    main(default_figures=['B'])
