import numpy as np
import matplotlib.pyplot as plt
from vorpy.src.calculations import rotate_points
from vorpy.src.calculations import calc_circ
from vorpy.src.visualize import *


def plot_vertex_2d(calc_verts, ob, pv_loc, edge_ball_locs, edge_ball_rads, edge_normal, real_vert):
    my_edge = calc_circ(*edge_ball_locs, *edge_ball_rads)

    # Calculate the edge plane
    ep_norm = np.cross(pv_loc - my_edge[0], edge_normal)
    new_chosen_vert_point = rotate_points(ep_norm, np.array([chosen_vert['loc']]))
    new_cv_points = rotate_points(ep_norm, np.array([_['loc'] for _ in calc_verts]))
    new_pv_point = rotate_points(ep_norm, np.array([pv_loc]))
    new_oa_point = rotate_points(ep_norm, np.array([ob]))
    new_real_vert_loc = rotate_points(ep_norm, np.array([real_vert['loc']]))
    plt.scatter([_[0] for _ in new_cv_points], [_[1] for _ in new_cv_points])
    plt.scatter([new_pv_point[0][0]], [new_pv_point[0][1]], marker='o', s=20)
    plt.scatter([new_oa_point[0][0]], [new_oa_point[0][1]])
    plt.scatter([new_real_vert_loc[0][0]], [new_real_vert_loc[0][1]], marker='x', s=20)
    plt.scatter([new_chosen_vert_point[0][0]], [new_chosen_vert_point[0][1]], marker='.', s=40)

    plt.show()


def plot_vert_situation(edge_balls, my_vert, vn_1_loc, vn_1_rad, b_locs, b_rads, a0=None, a1=None, a2=None, v0=None):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    # actual_vert_other_ball = [_ for _ in actual_vert['balls'] if _ not in edge_ndxs][0]
    # edge balls
    plot_balls([b_locs[_] for _ in edge_balls], [b_rads[_] for _ in edge_balls], fig=fig, ax=ax, colors=['r', 'r', 'r'])
    # other ball
    if a0 is not None:
        plot_balls([b_locs[a0]], [b_rads[1779]], fig=fig, ax=ax, colors=['pink'])
    # interfering ball
    if a1 is not None:
        plot_balls([b_locs[3144]], [b_rads[3144]], fig=fig, ax=ax, colors=['orange'])
    # # actual vert other ball
    if a2 is not None:
        plot_balls([b_locs[7]], [b_rads[7]], fig=fig, ax=ax, colors=['purple'])
    # actual vert
    plot_verts([my_vert['loc2']], [my_vert['rad2']], fig=fig, ax=ax, spheres=True, colors=['b'])
    # closest vert
    if v0 is not None:
        plot_verts([v0['loc']], [v0['rad']], fig=fig, ax=ax, spheres=True, colors=['white'])
    # previous vert
    plot_verts([vn_1_loc], [vn_1_rad], fig=fig, ax=ax, spheres=True, colors=['green'], Show=True)


