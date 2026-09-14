"""Project atomic-property overrides shared by the viewer and solver adapter."""
from dataclasses import replace
from math import isfinite

from vorpy.src.chemistry import element_radii, my_masses

PROPERTIES = ("radius", "mass", "charge")
# Match the viewer's ion classification without importing its rendering stack.
ION_RESIDUES = {"LI", "NA", "K", "RB", "CS", "MG", "CA", "SR", "BA", "ZN", "CD",
                "FE", "FE2", "FE3", "MN", "CU", "CU1", "CU2", "CO", "NI", "AL", "F",
                "CL", "BR", "I", "IOD", "SO4", "PO4", "NH4"}


def validate_defaults(values):
    """Return a detached, normalized override mapping; reject invalid numbers."""
    result = {}
    for scope, properties in values.items():
        normalized = {}
        for name, value in properties.items():
            if name not in PROPERTIES:
                raise ValueError(f"Unknown atomic property: {name}")
            number = float(value)
            if not isfinite(number) or (name != "charge" and number <= 0):
                raise ValueError(f"{scope}: {name} must be {'finite' if name == 'charge' else 'positive and finite'}")
            normalized[name] = number
        if normalized:
            result[str(scope).upper()] = normalized
    return result


def element_defaults(element):
    return {"radius": element_radii.get(element.upper()),
            "mass": my_masses.get(element.lower()), "charge": None}


def effective_overrides(defaults, element, residue, name):
    """Resolve each property separately: residue atom > ion > element."""
    values = dict(defaults.get(element.upper(), {}))
    if residue.upper() in ION_RESIDUES:
        values.update(defaults.get(f"ION:{element.upper()}", {}))
    values.update(defaults.get(f"RES:{residue}:{name}".upper(), {}))
    return values


def apply_atom_defaults(atom, defaults):
    if not defaults and not atom.source_properties:
        return atom
    source = atom.source_properties or {
        "radius": atom.radius, "mass": atom.mass, "charge": atom.charge,
    }
    values = dict(source)
    values.update(effective_overrides(defaults, atom.element, atom.residue_name, atom.name))
    if all(getattr(atom, name) == value for name, value in values.items()):
        return atom
    return replace(atom, **values, source_properties=dict(source))


def apply_system_defaults(system, defaults):
    """Apply edits before Group creation or spatial indexing, leaving algorithms intact."""
    defaults = validate_defaults(defaults)
    if not defaults:
        return
    # PDB charge columns may use pandas string dtype; preserve source strings
    # while permitting explicit numeric charges. Numeric properties may be integer-typed.
    for name, column in (("radius", "rad"), ("mass", "mass"), ("charge", "charge")):
        if any(name in values for values in defaults.values()) and column in system.balls:
            system.balls[column] = system.balls[column].astype(object if name == "charge" else float)
    radius_changed = False
    for index, row in system.balls.iterrows():
        values = effective_overrides(defaults, str(row['element']),
                                     str(row.get('res_name', '')), str(row['name']))
        for name, value in values.items():
            column = 'rad' if name == 'radius' else name
            system.balls.at[index, column] = value
            radius_changed |= name == 'radius'
    if radius_changed:
        system.max_atom_rad = float(system.balls['rad'].max())
        system.invalidate_spatial_index()
