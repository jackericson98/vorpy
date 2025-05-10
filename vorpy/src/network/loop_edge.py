from vorpy.src.calculations import calc_circ
from vorpy.src.calculations import calc_dist


"""
Detects and calculates and loop edges
"""


def detect_loop_edge(ball_loc, ball_rad, sur_locs, sur_rads):
    """
    Detects if a ball and its surrounding balls form a loop edge by checking if they can create a valid vertex.

    Parameters
    ----------
    ball_loc : numpy.ndarray
        Location of the central ball in 3D space
    ball_rad : float
        Radius of the central ball
    sur_locs : list of numpy.ndarray
        List of locations of surrounding balls
    sur_rads : list of float
        List of radii of surrounding balls

    Returns
    -------
    bool
        True if a valid loop edge is detected, False otherwise

    Notes
    -----
    - A loop edge requires at least 3 balls to form a valid vertex
    - The function first finds the two closest balls to the central ball
    - The validity of the loop edge is determined by checking if these balls can form a circle
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
    