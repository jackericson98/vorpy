import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..', '..'))
sys.path.append(vorpy_root)

from vorpy.src.system.system import System
from vorpy.src.visualize.mpl_visualize import plot_balls, plot_verts, plot_edges, plot_surfs


"""
Complete Cell Plotting

Creates a random ensemble of balls, places them into a VorPy System,
creates one group containing the entire ensemble, solves that group's
network, selects a complete interior cell, and plots it.
"""


# =============================================================================
# Settings
# =============================================================================

# Random ensemble
seed = 5
num_balls = 60
ensemble_size = 12.0
min_rad = 0.7
max_rad = 1.9

# Network
net_type = 'aw'
surf_res = 0.2
surf_color = 'purple'
box_size = 1.0
max_vert = 40

# Plot
show_cell_atom = True
show_neighbor_atoms = True
show_all_atoms = False
show_surfs = True
show_edges = True
show_verts = True
show_vert_balls = False

atom_alpha = 0.05
all_atom_alpha = 0.05
surf_alpha = 0.35
edge_thickness = 2
vert_alpha = 0.15

plot_title = False
show_axes = False


# =============================================================================
# Generate Random Ensemble
# =============================================================================

rng = np.random.default_rng(seed)

locs = [np.array(_) for _ in rng.uniform(-ensemble_size / 2, ensemble_size / 2, size=(num_balls, 3))]
rads = rng.uniform(min_rad, max_rad, size=num_balls)

balls = pd.DataFrame({
    'loc': locs,
    'rad': rads,
    'num': list(range(num_balls)),
    'name': [str(i) for i in range(num_balls)],
    'mass': [None] * num_balls
})


# =============================================================================
# Create the VorPy System
# =============================================================================

my_sys = System(
    spheres=balls,
    simple=True,
    name='Random Ensemble'
)


# =============================================================================
# Create One Group Containing the Entire Ensemble
# =============================================================================

all_balls = list(range(num_balls))

my_sys.create_group(
    atoms=all_balls,
    make_net=True
)

group = my_sys.groups[0]
group.name = 'Random Ensemble'


# =============================================================================
# Build the Group Network
# =============================================================================

net = group.net

net.default_settings(
    surf_res=surf_res,
    box_size=box_size,
    max_vert=max_vert,
    net_type=net_type,
    print_metrics=False
)

net.build()


# =============================================================================
# Find a Complete Cell
# =============================================================================

complete_cells = [i for i in net.balls.index if bool(net.balls.at[i, 'complete'])]

if len(complete_cells) == 0:
    raise RuntimeError(
        'No complete cells were found in this random ensemble. '
        'Try increasing num_balls or changing seed.'
    )

# Pick the complete cell closest to the center of the ensemble
ensemble_center = np.mean(locs, axis=0)

cell_ndx = min(
    complete_cells,
    key=lambda i: np.linalg.norm(np.asarray(net.balls.at[i, 'loc']) - ensemble_center)
)

print(f'\nSelected cell: {cell_ndx}')
print(f'Complete cells found: {len(complete_cells)}')


# =============================================================================
# Get the Network Elements Belonging to the Cell
# =============================================================================

vert_ndxs = list(net.balls.at[cell_ndx, 'verts'])
edge_ndxs = list(net.balls.at[cell_ndx, 'edges'])
surf_ndxs = list(net.balls.at[cell_ndx, 'surfs'])

print(f'Vertices: {len(vert_ndxs)}')
print(f'Edges:    {len(edge_ndxs)}')
print(f'Surfaces: {len(surf_ndxs)}')


# =============================================================================
# Find the Balls that Actually Define the Cell
# =============================================================================

neighbor_ndxs = set()

for surf_ndx in surf_ndxs:
    for ball_ndx in net.surfs.at[surf_ndx, 'balls']:
        if ball_ndx != cell_ndx:
            neighbor_ndxs.add(ball_ndx)

neighbor_ndxs = sorted(neighbor_ndxs)

cell_ball_ndxs = [cell_ndx] + neighbor_ndxs

print(f'Neighbor balls: {neighbor_ndxs}')
print(f'Balls needed for cell: {cell_ball_ndxs}')


# =============================================================================
# Create Plot
# =============================================================================

fig = plt.figure()
ax = fig.add_subplot(projection='3d')


# =============================================================================
# Plot Full Original Ensemble
# =============================================================================

if show_all_atoms:
    plot_balls(
        alocs=[net.balls.at[i, 'loc'] for i in net.balls.index],
        arads=[net.balls.at[i, 'rad'] for i in net.balls.index],
        fig=fig,
        ax=ax,
        alpha=all_atom_alpha,
        res=10
    )


# =============================================================================
# Plot Cell Surfaces
# =============================================================================

if show_surfs:
    surf_points = [net.surfs.at[i, 'points'] for i in surf_ndxs]
    surf_tris = [net.surfs.at[i, 'tris'] for i in surf_ndxs]

    plot_surfs(
        surf_points,
        surf_tris,
        fig=fig,
        ax=ax,
        alpha=surf_alpha,
        colors=[surf_color] * len(surf_ndxs)
    )


# =============================================================================
# Plot Cell Edges
# =============================================================================

if show_edges:
    edge_points = [net.edges.at[i, 'points'] for i in edge_ndxs]

    plot_edges(
        edge_points,
        fig=fig,
        ax=ax,
        colors=['blue'] * len(edge_ndxs),
        thickness=edge_thickness
    )


# =============================================================================
# Plot Cell Vertices
# =============================================================================

if show_verts:
    vert_locs = [net.verts.at[i, 'loc'] for i in vert_ndxs]
    vert_rads = [abs(net.verts.at[i, 'rad']) for i in vert_ndxs]

    plot_verts(
        vert_locs,
        vert_rads,
        fig=fig,
        ax=ax,
        colors=['red'] * len(vert_ndxs),
        spheres=show_vert_balls,
        alpha=vert_alpha,
        res=10
    )


# =============================================================================
# Plot Central Cell Atom
# =============================================================================

if show_cell_atom:
    plot_balls(
        alocs=[net.balls.at[cell_ndx, 'loc']],
        arads=[net.balls.at[cell_ndx, 'rad']],
        fig=fig,
        ax=ax,
        colors=['black'],
        alpha=atom_alpha,
        res=10
    )


# =============================================================================
# Plot Neighboring Atoms
# =============================================================================

if show_neighbor_atoms:
    plot_balls(
        alocs=[net.balls.at[i, 'loc'] for i in neighbor_ndxs],
        arads=[net.balls.at[i, 'rad'] for i in neighbor_ndxs],
        fig=fig,
        ax=ax,
        alpha=atom_alpha,
        res=10
    )


# =============================================================================
# Set Plot Limits Around the Selected Cell
# =============================================================================

plot_ball_ndxs = [cell_ndx] + neighbor_ndxs
plot_locs = np.asarray([net.balls.at[i, 'loc'] for i in plot_ball_ndxs])
plot_rads = np.asarray([net.balls.at[i, 'rad'] for i in plot_ball_ndxs])

mins = np.min(plot_locs - plot_rads[:, None], axis=0)
maxs = np.max(plot_locs + plot_rads[:, None], axis=0)

center = (mins + maxs) / 2
half_range = max(maxs - mins) / 2 * 1.15

ax.set_xlim(center[0] - half_range, center[0] + half_range)
ax.set_ylim(center[1] - half_range, center[1] + half_range)
ax.set_zlim(center[2] - half_range, center[2] + half_range)


# =============================================================================
# Axes / Title
# =============================================================================

if not show_axes:
    ax.set_axis_off()

if plot_title:
    ax.set_title(
        f'Random Complete Cell {cell_ndx}\n'
        f'{len(surf_ndxs)} Surfaces, {len(edge_ndxs)} Edges, {len(vert_ndxs)} Vertices'
    )


# =============================================================================
# Show
# =============================================================================

plt.show()