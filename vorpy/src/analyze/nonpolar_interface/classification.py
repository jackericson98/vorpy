"""Explicit, versioned chemical and interface labels for analysis tables."""

from vorpy.src.chemistry.chemistry_interpreter import (
    ion_residue_names,
    water_residue_names,
)


NONPOLAR_RULE_VERSION = "strict_element_CH_v1"
EXCLUDED_RESIDUES = frozenset(
    str(name).upper() for name in water_residue_names | ion_residue_names
)


def classify_atom(atom, rule=NONPOLAR_RULE_VERSION):
    """Return ``nonpolar``, ``other``, or ``excluded`` for one atom record.

    The strict baseline deliberately uses C/H plus residue exclusions. It is
    a transparent screening rule, not a solvent-accessibility or force-field
    hydrophobicity assignment.
    """
    if rule != NONPOLAR_RULE_VERSION:
        raise ValueError(f"Unknown nonpolar classification rule: {rule!r}")
    residue = str(atom.get("res_name", atom.get("residue", ""))).strip().upper()
    element = str(atom.get("element", "")).strip().upper()
    if residue in EXCLUDED_RESIDUES or element in {"M", "MW"}:
        return "excluded"
    if element in {"C", "H"}:
        return "nonpolar"
    return "other"


def classify_surface_pair(atom_1, atom_2, rule=NONPOLAR_RULE_VERSION):
    """Classify a surface pair while retaining mixed/unknown chemistry."""
    labels = (classify_atom(atom_1, rule), classify_atom(atom_2, rule))
    if "excluded" in labels:
        return "excluded"
    if labels == ("nonpolar", "nonpolar"):
        return "nonpolar"
    if "nonpolar" in labels:
        return "mixed"
    return "other"
