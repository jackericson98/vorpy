import numpy as np
from System.sys_funcs.calcs.vert import calc_vert
from System.sys_funcs.calcs.edge import calc_circ
from System.sys_funcs.calcs.calcs import calc_dist


"""
Detects and calculates and loop edges
"""


def detect_loop_edge(ball_loc, ball_rad, sur_locs, sur_rads):
    """
    Takes in a ball and determines if any of the surrounding balls make a vertex
    """
    # Initial condition where you at least need three balls
    if len(sur_locs) < 2:
        return False
    # Closest balls variable instantiation
    c_locs, c_rads = sur_locs[:2], sur_rads[:2]
    c_dists = [calc_dist(ball_loc, sur_locs[i]) - (sur_rads[i] + ball_rad) for i in range(2)]
    close_balls = sorted(zip(c_locs, c_rads, c_dists), key=lambda x: x[2])
    # Find the closest balls to the ball
    for loc, rad in zip(sur_locs[2:], sur_rads[2:]):
        # Calculate the distance between the loc and the ball loc
        b_dist = calc_dist(loc, ball_loc) - rad - ball_rad
        if b_dist < close_balls[0][2]:
            # Reassign the close_balls
            close_balls = [(loc, rad, b_dist), close_balls[0]]
        elif b_dist < close_balls[1][2]:
            # Reassign the close_balls
            close_balls = [close_balls[0], (loc, rad, b_dist)]

    # Determine the edge loop if one exists
    my_circ = calc_circ(*[_[0] for _ in close_balls], ball_loc, *[_[1] for _ in close_balls], True)
    # Determine if there is not a loop edge
    