from Data.Analyze.read_logs import read_logs
from System.sys_funcs.input.pdb import read_pdb


"""
The point of this file is to be able to take in multiple logs from different solved models and generate compariso data

Things to compare:
    1. Different solutions - Additively Weighted, Power, Primitive
    2. Different coarse graining schemes
    3. Different versions of the same molecule - frames, hybrid, bound, pressurized
    4. Different groupings - Atom, element type, backbone, side-chains, residues, at
"""


def analyze_atoms():
    """

    """
    pass


def analyze_elements(elements):
    pass


def analyze_residues(residues):
    pass


def compare_by_chain(chains):

    pass


def compare_interfaces():
    pass


def compare_models(log_files, pdb_files):
    """
    Takes in logs and pdbs and returns a comparison breakdowns
    e.g. edta vor and edta pow, grouping - residue. Each residue will have their volume, surface area
    """
    pass


def compare_models_network_type(models, grouping):
    """
    Compares vor, pow, del models to one another and produces a grouping level % diff in vol sa and max curvature
    """
    # Check to make sure that the models are the same
    same = True
    model_name = models[0]['data']['name']
    for model in models:

