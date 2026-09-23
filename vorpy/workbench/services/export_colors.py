"""Map viewer colors onto backend geometry without changing the viewer."""
from dataclasses import replace

import numpy as np
from matplotlib import colormaps
from matplotlib.colors import to_rgb

from vorpy.workbench.services.display_colors import DisplayColors


def make_color_provider(result, overrides):
    from vorpy.workbench.services.vorpy_backend import _layers_from_network
    group = result.export_group
    network = group.net
    if all(getattr(network, key, None) is None for key in ('surfs', 'edges', 'verts')):
        return None
    if group.layer_surfs is None and all(getattr(network, key, None) is not None for key in ('surfs', 'edges', 'verts')):
        group.get_layers(max_layers=1)
    shell_indices = {
        'surfaces': list(group.layer_surfs[0]) if group.layer_surfs else [],
        'edges': list(group.layer_edges[0]) if group.layer_edges else [],
        'vertices': list(group.layer_verts[0]) if group.layer_verts else [],
    }
    # The result can omit geometry hidden at solve time. Generate fallback layers
    # for those exports, but use the current viewer layers wherever available.
    fallback = _layers_from_network(network, shell_indices['surfaces'],
                                   shell_indices['edges'], shell_indices['vertices'])
    current = {layer.name: layer for layer in result.layers}
    layers = [current.get(layer.name, layer) for layer in fallback]
    normalizer = DisplayColors(result.layers or layers)
    tables = {'surfaces': network.surfs, 'edges': network.edges, 'vertices': network.verts}
    colors = {}
    for layer in layers:
        kind = layer.kind
        mode = 'boundary' if 'shell' in layer.name.lower() else 'magnitude'
        indices = shell_indices[kind] if mode == 'boundary' else list(range(len(tables[kind])))
        if kind in overrides:
            layer = replace(layer, color_scheme='solid', color=overrides[kind])
        raw = layer.cell_scalars.get(layer.color_scheme)
        if layer.color_scheme != 'solid' and raw is not None:
            values, (low, high) = normalizer._display_scalars(layer, raw)
            rgb = colormaps[layer.color_map](np.clip((values - low) / (high - low), 0, 1))[:, :3]
        else:
            rgb = None
        offset = 0
        for index in indices:
            row = tables[kind].iloc[index]
            count = len(row['tris']) if kind == 'surfaces' else max(len(row['points']) - 1, 0) if kind == 'edges' else 1
            if rgb is None or count == 0:
                value = np.asarray([to_rgb(layer.color)])
            else:
                value = rgb[offset:offset + count] if kind == 'surfaces' else rgb[offset:offset + 1]
            colors[kind, mode, index] = value
            offset += count

    def provider(kind, index, mode, count):
        # Cell files show the same component colors as the all-geometry viewer
        # layer; shell files use the independently styled shell viewer layer.
        display_mode = 'boundary' if mode == 'boundary' else 'magnitude'
        value = colors.get((kind, display_mode, index))
        if value is None:
            value = colors[kind, 'magnitude', index]
        if len(value) == count:
            return value
        return np.tile(value[0], (count, 1))
    return provider
