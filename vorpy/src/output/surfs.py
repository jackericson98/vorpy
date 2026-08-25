import os
import time
import numpy as np
from vorpy.src.output.color_tris import color_tris


def write_surfs(net, surfs, file_name, color=False, directory=None, concave_colors=False, ref_surfs=None,
                universal_max=True, chunk_size=10000):
    """Export selected network surfaces to an OFF file."""

    if directory is not None:
        os.chdir(directory)

    if surfs is None or len(surfs) == 0:
        return

    surfs = list(surfs)

    if color is False:
        color = (1, 0, 0)

    if ref_surfs is None:
        ref_surfs = []

    if 'tri_colors' not in net.surfs:
        net.surfs['tri_colors'] = [[] for _ in range(len(net.surfs))]

    # ------------------------------------------------------------------
    # Gather surfaces
    # ------------------------------------------------------------------

    surf_rows = [net.surfs.iloc[ndx] for ndx in surfs]

    # ------------------------------------------------------------------
    # Count geometry
    # ------------------------------------------------------------------
    num_points = sum(len(surf['points']) for surf in surf_rows)
    num_tris = sum(len(surf['tris']) for surf in surf_rows)

    # ------------------------------------------------------------------
    # Determine curvature scale
    # ------------------------------------------------------------------

    surf_scheme = net.settings['surf_scheme'].lower()

    # Map surface coloring schemes to their DataFrame columns
    scheme_columns = {
        'mean': 'mean_curv',
        'mean_curv': 'mean_curv',

        'gauss': 'gauss_curv',
        'gauss_curv': 'gauss_curv',

        'avg_mean': 'avg_mean_curv',
        'avg_gauss': 'avg_gauss_curv',

        'max_mean': 'mean_curv',
        'max_gauss': 'gauss_curv',

        'int_mean_curv': 'int_mean_curv',
        'int_mean_curv_sq': 'int_mean_curv_sq',
        'int_gauss_curv': 'int_gauss_curv',

        # Representative surface energy
        'surf_energy': 'surf_energy',
    }

    value_column = scheme_columns.get(
        surf_scheme,
        'mean_curv'
    )

    if universal_max:
        values = np.asarray(
            net.surfs[value_column],
            dtype=float
        )
    else:
        values = np.asarray(
            [surf[value_column] for surf in surf_rows],
            dtype=float
        )

    # Remove NaN / inf
    values = values[np.isfinite(values)]

    if len(values) == 0:
        min_val = 0.0
        max_val = 1.0

    else:
        min_val = float(np.min(values))
        max_val = float(np.max(values))

        # Prevent divide-by-zero during normalization
        if max_val == min_val:
            max_val = min_val + 1.0

    # ------------------------------------------------------------------
    # Generate triangle colors
    # ------------------------------------------------------------------
    tri_colors = []

    if net.settings['net_type'] == 'aw':
        ref_set = set(ref_surfs)

        for surf in surf_rows:
            if concave_colors:
                ref_ball = next(ball for ball in surf['balls'] if ball in ref_set)
                non_ref_ball = next(ball for ball in surf['balls'] if ball not in ref_set)
                inverse = net.balls.iloc[ref_ball]['rad'] <= net.balls.iloc[non_ref_ball]['rad']

                tri_colors.append(color_tris(surf=surf, color_map=net.settings['surf_col'],
                                             color_scheme=net.settings['surf_scheme'],
                                             color_factor=net.settings['scheme_factor'], max_val=max_val,
                                             min_val=min_val))
            else:
                tri_colors.append(color_tris(
                    surf=surf, color_map=net.settings['surf_col'], color_scheme=net.settings['surf_scheme'],
                    color_factor=net.settings['scheme_factor'], max_val=max_val
                ))
    else:
        tri_colors = [[color] * len(surf['tris']) for surf in surf_rows]

    # ------------------------------------------------------------------
    # Write OFF file
    # ------------------------------------------------------------------

    with open(file_name + ".off", 'w', buffering=1024 * 1024) as file:
        file.write(f"OFF\n{num_points} {num_tris} 0\n\n\n")

        # --------------------------------------------------------------
        # Write points in chunks
        # --------------------------------------------------------------

        buffer = []

        for surf in surf_rows:
            for point in surf['points']:
                buffer.append(
                    f"{round(float(point[0]), 4)} "
                    f"{round(float(point[1]), 4)} "
                    f"{round(float(point[2]), 4)}\n"
                )

                if len(buffer) >= chunk_size:
                    file.write(''.join(buffer))
                    buffer.clear()

        if buffer:
            file.write(''.join(buffer))
            buffer.clear()

        # --------------------------------------------------------------
        # Write faces in chunks
        # --------------------------------------------------------------

        vertex_offset = 0
        buffer = []

        for surf, colors in zip(surf_rows, tri_colors):
            for tri, tri_color in zip(surf['tris'], colors):
                buffer.append(
                    f"3 {tri[0] + vertex_offset} "
                    f"{tri[1] + vertex_offset} "
                    f"{tri[2] + vertex_offset} "
                    f"{tri_color[0]} {tri_color[1]} {tri_color[2]}\n"
                )

                if len(buffer) >= chunk_size:
                    file.write(''.join(buffer))
                    buffer.clear()

            vertex_offset += len(surf['points'])

        if buffer:
            file.write(''.join(buffer))


def write_surfs1(surfs, file_name, settings, color=False, directory=None, chunk_size=10000):
    """Legacy DataFrame-based OFF surface writer."""

    if directory is not None:
        os.chdir(directory)

    if surfs is None or len(surfs) == 0:
        return

    if color is False:
        color = np.random.rand(3)

    surf_rows = [surf for _, surf in surfs.iterrows()]
    num_points = sum(len(surf['points']) for surf in surf_rows)
    num_tris = sum(len(surf['tris']) for surf in surf_rows)

    if settings['net_type'] == 'aw':
        max_val = max(surfs['curv'])
        tri_colors = [
            color_tris(
                surf=surf, color_map=settings['surf_col'], color_scheme=settings['surf_scheme'],
                color_factor=settings['scheme_factor'], max_val=max_val
            )
            for surf in surf_rows
        ]
    else:
        tri_colors = [[color] * len(surf['tris']) for surf in surf_rows]

    with open(file_name + ".off", 'w', buffering=1024 * 1024) as file:
        file.write(f"OFF\n{num_points} {num_tris} 0\n\n\n")

        buffer = []

        for surf in surf_rows:
            for point in surf['points']:
                buffer.append(
                    f"{round(float(point[0]), 4)} "
                    f"{round(float(point[1]), 4)} "
                    f"{round(float(point[2]), 4)}\n"
                )

                if len(buffer) >= chunk_size:
                    file.write(''.join(buffer))
                    buffer.clear()

        if buffer:
            file.write(''.join(buffer))
            buffer.clear()

        vertex_offset = 0

        for surf, colors in zip(surf_rows, tri_colors):
            for tri, tri_color in zip(surf['tris'], colors):
                buffer.append(
                    f"3 {tri[0] + vertex_offset} "
                    f"{tri[1] + vertex_offset} "
                    f"{tri[2] + vertex_offset} "
                    f"{tri_color[0]} {tri_color[1]} {tri_color[2]}\n"
                )

                if len(buffer) >= chunk_size:
                    file.write(''.join(buffer))
                    buffer.clear()

            vertex_offset += len(surf['points'])

        if buffer:
            file.write(''.join(buffer))