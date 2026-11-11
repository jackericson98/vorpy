"""Configuration shared by the F5 loader, runner, and figures."""
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Config:
    dataset_root: Path = Path(r'C:\Users\Ryzen-3600\Versiti Dropbox\John Ericson\Versiti lab\Jack\Ericson Voronoi DNA\OrganizedModels')
    output_root: Path = Path(__file__).parent / 'output'
    polymers: tuple = ('PEG', 'Heparin')
    models: tuple = ()
    frames: tuple = ()
    schemes: tuple = ('aw', 'pow', 'prm')
    water_residues: tuple = ('SOL', 'HOH', 'WAT', 'H2O', 'TIP3', 'TIP3P', 'SPC', 'SPCE')
    water_oxygens: tuple = ('OW', 'O', 'OH2')
    water_hydrogens: tuple = ('HW1', 'HW2', 'H1', 'H2', 'H01', 'H02', '1H', '2H')
    ion_residues: tuple = ('SOD', 'NA', 'NA+', 'CLA', 'CL', 'CL-', 'K', 'POT', 'MG', 'CA', 'CAL')
    other_solvents: tuple = ('DMS', 'DMSO', 'ETH', 'MEO', 'GOL')
    peg_residues: tuple = ('PEG', 'EG_', 'EG')
    heparin_residues: tuple = ('AGLC', 'AIDO', 'AGL', 'AID')
    buried_threshold: float = 0.10
    exposed_threshold: float = 0.40
    area_rtol: float = 1e-4
    area_atol: float = 2e-5
    minimum_interface_area: float = 0.0
    coordinate_tolerance: float = 0.002
    include_hydrogens: bool = True
    include_incomplete: bool = False
    figure_size: tuple = (12., 7.)
    dpi: int = 200
    font_size: int = 10
    colors: tuple = (('aw', '#2166ac'), ('pow', '#d73027'), ('prm', '#1a9850'))
    export_formats: tuple = ('png', 'svg')
    show: bool = False
    cache: bool = True


SCHEME_LABELS = {'aw': 'AW', 'pow': 'Power', 'prm': 'Primitive'}
KEY = ['polymer', 'model', 'frame']
ATOM_KEY = KEY + ['atom_index']
UNIT_KEY = KEY + ['unit_id']
METRICS = ['water_fraction', 'water_area', 'volume', 'surface_area',
           'water_coordination', 'polymer_coordination']
UNITS = {'water_fraction': 'Water SA fraction', 'water_area': 'Water interface area (Å²)',
         'volume': 'Atomic volume (Å³)', 'mean_atomic_volume': 'Mean atomic volume (Å³)',
         'surface_area': 'Atomic surface area (Å²)',
         'water_coordination': 'Unique water molecules',
         'polymer_coordination': 'Polymer atom neighbors',
         'area_per_water': 'Water interface area / water molecule (Å²)',
         'buried_fraction': 'Fraction of repeat units buried'}
