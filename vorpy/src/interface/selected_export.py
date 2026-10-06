"""Export exactly the physical representation cached by Phase-1 analysis."""
import json
from pathlib import Path

from vorpy.src.geometry.visualization.export import (
    _write_interface_layers, _write_csv, pymol_off_loader_lines,
)
from vorpy.src.output.draw import DEFAULT_EDGE_RADIUS


def export_selected_interface(iface):
    analysis = getattr(iface, 'geometry_analysis', None)
    rep = getattr(iface, '_geometry_representation', None)
    if analysis is None or rep is None:
        raise ValueError('Selected geometry requires cached InterfaceGeometryAnalysis')
    pairs = {tuple(pair) for pair in analysis.alpha_selection['selected_generator_pairs']}
    mappings = {(int(r['generator_i']), int(r['generator_j'])) for r in rep.selection_mappings}
    physical = {tuple(sorted(s.generator_ids)) for s in rep.selected_surfaces}
    mapped = {(int(r['generator_i']), int(r['generator_j'])) for r in rep.selection_mappings if r['selected']}
    if pairs != mappings or physical != mapped:
        raise ValueError('Cached analysis, selected surfaces, and mapping disagree')
    for surface in rep.selected_surfaces:
        row = iface.net.surfs.loc[surface.feature_id]
        if tuple(sorted(map(int, row['balls']))) != tuple(sorted(surface.generator_ids)):
            raise ValueError('Selected surface does not match its generator pair')
    root = Path(iface.dir).resolve() / 'alpha_interface'
    coordinates = {int(i): row['loc'] for i, row in iface.net.balls.iterrows()}
    layers = _write_interface_layers(iface.net, rep, root, 'selected', coordinates, DEFAULT_EDGE_RADIUS)
    _write_csv(root / 'mapping.csv', rep.selection_mappings)
    metadata = {'scheme': rep.scheme, 'alpha': rep.alpha,
                'selected_generator_pairs': sorted(pairs), 'physical_generator_pairs': sorted(physical),
                'selected_surface_ids': [s.feature_id for s in rep.selected_surfaces],
                'additional_network_solves': 0, 'area_A2': rep.area_A2,
                'full_physical_interface': '../surfs.off', 'layers': layers}
    (root / 'selection.json').write_text(json.dumps(metadata, indent=2) + '\n', encoding='utf-8')
    layer = layers['selected']
    lines = ['import os', 'from pymol import cmd',
             'from pymol.cgo import BEGIN, END, TRIANGLES, COLOR, VERTEX',
             *pymol_off_loader_lines()]
    for key, name, color in [('surfaces', 'alpha_interface', (1., .15, .25)),
                             ('physical_edges', 'alpha_edges', (.3, .3, .3)),
                             ('physical_vertices', 'alpha_vertices', (.9, .6, .1))]:
        if layer[key]:
            lines.append(f'_load_off_layer({str(root / layer[key])!r}, {name!r}, {color!r})')
    lines.append("cmd.zoom('alpha_interface')")
    (root / 'selected_interface.py').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    (root / 'selected_interface.pml').write_text(f'run "{(root / "selected_interface.py").as_posix()}"\n', encoding='utf-8')
    return metadata
