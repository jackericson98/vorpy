import time
import numpy as np
from vorpy.src.calculations import calc_dist
from vorpy.src.calculations import calc_com
from vorpy.src.calculations import calc_angle_jit
from vorpy.src.calculations import calc_circ
from vorpy.src.calculations import calc_edge_dir
from vorpy.src.network.edge_project import edge_project
from vorpy.src.calculations import calc_surf_func
from vorpy.src.visualize import plot_edges
from vorpy.src.visualize import plot_balls
from vorpy.src.visualize import plot_verts
import matplotlib.pyplot as plt


def build_straight_edge(locs, rads, vlocs, res):
    try:
        loc, rad = calc_circ(locs[0], locs[1], locs[2], rads[0], rads[1], rads[2])
    except TypeError:
        loc = calc_com([locs[0], locs[1], locs[2]])
        rad = calc_dist(loc, locs[0]) - rads[0]
    # Create the vals dictionary
    vals = {'loc': loc, 'rad': rad}
    # Determine the edge length
    edge_dist = calc_dist(vlocs[0], vlocs[1])
    # Divide the edge length by the resolution to find the number of points
    num_points = max(int(edge_dist / res) + 1, 3)
    # Create the new resolution to get even divisions of the edge
    new_res = edge_dist / num_points
    # Find the direction the edge heads
    edge_dir = vlocs[1] - vlocs[0]
    edge_dir_norm = np.linalg.norm(edge_dir)

    if edge_dir_norm == 0 or not np.isfinite(edge_dir_norm):
        return [vlocs[0], vlocs[1]], vals

    e_hat = edge_dir / edge_dir_norm
    e_points = [vlocs[0] + i * new_res * e_hat for i in range(num_points + 1)]
    # Return the edge
    return e_points, vals


def mid_edge_point(ep1, ep2, func, vmid, direction, new_direction=True):
    """
    Calculates a middle point on an edge based on provided edge points, function, and direction.

    Parameters:
        ep1, ep2 (tuple): Edge points between which to calculate the middle point.
        func (function): Function describing the surface on which the edge lies.
        vmid (np.ndarray): Midpoint used for reference in calculations.
        direction (tuple): Initial direction for edge calculation.
        new_direction (bool, optional): Flag to indicate if direction calculation is required; default True.

    Returns:
        ndarray: New edge point projected onto the surface defined by 'func'.
    """
    # If the point is the first point we just need to move in the direction of the direction vector
    if new_direction:
        # Get the direction between the edges
        edir = ep2 - ep1
        # Get the distance between the points
        edist = calc_dist(ep1, ep2)

        if edist == 0 or not np.isfinite(edist):
            return None

        proj_point = ep1 + 0.5 * edir
        rn = proj_point - vmid
        rn_norm = np.linalg.norm(rn)

        if rn_norm == 0 or not np.isfinite(rn_norm):
            return None

        direction = rn / rn_norm

    direction = np.array(direction, dtype=float)
    dnorm = np.linalg.norm(direction)

    if dnorm == 0 or not np.isfinite(dnorm):
        return None

    direction = direction / dnorm

    return edge_project(np.array(direction), np.array(vmid), np.array(func))


def choose_stable_dnorm(vlocs, func, vmid, dnorm):
    mid_a = mid_edge_point(
        vlocs[0],
        vlocs[1],
        func,
        vmid,
        dnorm,
        new_direction=False
    )

    mid_b = mid_edge_point(
        vlocs[0],
        vlocs[1],
        func,
        vmid,
        -dnorm,
        new_direction=False
    )

    if mid_a is None or not np.all(np.isfinite(mid_a)):
        return -dnorm

    if mid_b is None or not np.all(np.isfinite(mid_b)):
        return dnorm

    score_a = calc_dist(vlocs[0], mid_a) + calc_dist(mid_a, vlocs[1])
    score_b = calc_dist(vlocs[0], mid_b) + calc_dist(mid_b, vlocs[1])

    return -dnorm if score_b < score_a else dnorm


def build_edge(
        locs,
        rads,
        vlocs,
        res,
        blocs,
        brads,
        eballs,
        straight=False,
        vmid=None,
        dnorm=None,
        edub=False,
        edge_points1=None,
        edge_verts=None,
        redone_edge=False,
        edge_index=None,
        debug=False,
        timeout=5.0
):
    """
    Constructs an edge based on various parameters describing the geometry and properties of the network elements.

    Parameters:
        locs (list): Locations of interest points.
        rads (list): Radii corresponding to each location.
        vlocs (list): Vertex locations defining the bounds of the edge.
        blocs (list): The balls in the network's locations
        brads (list): The balls in the network's radd
        res (float): Resolution for determining the detail of the edge computation.
        blocs, brads (list): Additional network-specific parameters, locations, and radii used in complex edge calculations.
        eballs (list): Indices or identifiers for the balls involved in edge calculation.
        straight (bool): Whether the edge should be constructed as a straight line.
        vmid (tuple): Midpoint from which to project
        dnorm (tuple): Normal direction for dynamic calculation segments.
        edub (bool): Indicates whether to use a double precision or higher accuracy mode.
        edge_points1 (list): Previously calculated edge points for visualization or debugging.
        edge_verts (list): Vertices associated with the edge for visualization or debugging.
        redone_edge (bool): Flag indicating whether the edge is being recalculated.

    Returns:
        tuple: A tuple containing the list of computed edge points and additional values for further processing.
    """

    def dprint(msg):
        if debug:
            print(msg)

    vlocs = [np.array(vlocs[0], dtype=float), np.array(vlocs[1], dtype=float)]
    locs = [np.array(_, dtype=float) for _ in locs]
    rads = [float(_) for _ in rads]

    if straight or round(rads[0], 3) == round(rads[1], 3) == round(rads[2], 3):
        return build_straight_edge(locs, rads, vlocs, res)

    # Choose a curved surface to project onto. If the edge isn't straight at least 2 surfs are curved.
    if round(rads[0], 10) == round(rads[1], 10):
        func = calc_surf_func(locs[1], rads[1], locs[2], rads[2])
    else:
        func = calc_surf_func(locs[0], rads[0], locs[1], rads[1])

    # Get the edge direction
    edge_vals = None

    if vmid is None or dnorm is None:

        if calc_dist(vlocs[0], vlocs[1]) < 1e-8:
            return [vlocs[0], vlocs[1]], {
                "loc": np.array(vlocs[0], dtype=float),
                "rad": 0.0,
                "case": "zero_length"
            }

        edge_vals = calc_edge_dir(blocs, brads, eballs, vlocs, edub=edub)

        if edge_vals is None:
            return build_straight_edge(locs, rads, vlocs, res)

        vmid = edge_vals["vmid"]
        dnorm = edge_vals["dnorm"]
        dnorm = choose_stable_dnorm(vlocs, func, vmid, dnorm)

    else:
        vmid = np.array(vmid, dtype=float)
        dnorm = np.array(dnorm, dtype=float)
        dnorm_norm = np.linalg.norm(dnorm)

        if dnorm_norm == 0 or not np.isfinite(dnorm_norm):

            if calc_dist(vlocs[0], vlocs[1]) < 1e-8:
                return [vlocs[0], vlocs[1]], {
                    "loc": np.array(vlocs[0], dtype=float),
                    "rad": 0.0,
                    "case": "zero_length"
                }

            edge_vals = calc_edge_dir(blocs, brads, eballs, vlocs, edub=edub)
            vmid = edge_vals["vmid"]
            dnorm = edge_vals["dnorm"]
            dnorm = choose_stable_dnorm(vlocs, func, vmid, dnorm)
        else:
            dnorm = dnorm / dnorm_norm

    dprint(f"\nBUILD_EDGE edge={edge_index}, eballs={eballs}, redone={redone_edge}")
    dprint(f"vlocs distance: {calc_dist(vlocs[0], vlocs[1])}")
    dprint(f"rads: {rads}")
    dprint(f"vmid: {vmid}")
    dprint(f"dnorm: {dnorm}")

    if edge_points1 is not None and edge_vals is not None:
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        plot_edges([edge_points1], fig, ax)
        plot_balls([blocs[_] for _ in eballs], [brads[_] for _ in eballs], fig=fig, ax=ax)
        ax.plot(
            [edge_vals['vmid'][0], edge_vals['vmid'][0] + edge_vals['dnorm'][0]],
            [edge_vals['vmid'][1], edge_vals['vmid'][1] + edge_vals['dnorm'][1]],
            [edge_vals['vmid'][2], edge_vals['vmid'][2] + edge_vals['dnorm'][2]]
        )
        plot_balls([edge_vals['loc']], [edge_vals['rad']], fig=fig, ax=ax, colors=['red'])
        plot_verts(vlocs, [1, 1], fig=fig, ax=ax, colors=['g', 'g'])
        print(edge_vals)
        print(vlocs)
        print(edge_verts)
        plt.show()

    if edge_vals is not None and edge_vals.get("case") == 5:
        edge0 = build_edge(
            locs, rads, [vlocs[0], edge_vals['loc']], res, blocs, brads, eballs,
            straight=straight,
            vmid=edge_vals['vmid0'],
            dnorm=edge_vals['dnorm0'],
            edub=edub,
            edge_verts=edge_verts,
            edge_index=edge_index,
            debug=debug,
            timeout=timeout
        )

        edge = build_edge(
            locs, rads, [edge_vals['loc'], edge_vals['loc2']], res, blocs, brads, eballs,
            straight=straight,
            vmid=edge_vals['vmid'],
            dnorm=edge_vals['dnorm'],
            edub=edub,
            edge_verts=edge_verts,
            edge_index=edge_index,
            debug=debug,
            timeout=timeout
        )

        edge1 = build_edge(
            locs, rads, [edge_vals['loc2'], vlocs[1]], res, blocs, brads, eballs,
            straight=straight,
            vmid=edge_vals['vmid1'],
            dnorm=edge_vals['dnorm1'],
            edub=edub,
            edge_verts=edge_verts,
            edge_index=edge_index,
            debug=debug,
            timeout=timeout
        )

        return edge0[0] + edge[0] + edge1[0], edge_vals
    # Create a catch for the time
    start = time.perf_counter()
    # Instantiate the edge points list with the vertices
    e_points = [*vlocs]
    loop_n = 0

    while True:
        loop_n += 1
        loop_start = time.perf_counter()
        new_points_added = False
        insertions = 0
        max_seg = 0.0

        # Loop through the edge points
        i = 0
        while i < len(e_points) - 1:
            # Get the middle points
            ep1, ep2 = e_points[i], e_points[i + 1]
            seg_len = calc_dist(ep1, ep2)
            max_seg = max(max_seg, seg_len)

            if seg_len > res:
                mid_point = mid_edge_point(
                    ep1,
                    ep2,
                    func,
                    vmid,
                    dnorm,
                    new_direction=len(e_points) > 2
                )

                if mid_point is None or not np.all(np.isfinite(mid_point)):
                    # Degenerate midpoint fallback.
                    mid_point = 0.5 * (ep1 + ep2)

                e_points.insert(i + 1, mid_point)
                new_points_added = True
                insertions += 1
                i += 1

            i += 1

        if debug:
            dprint(
                f"edge={edge_index}, redone={redone_edge}, "
                f"loop={loop_n}, points={len(e_points)}, "
                f"insertions={insertions}, max_seg={max_seg:.6f}, "
                f"loop_time={time.perf_counter() - loop_start:.3f}s"
            )

        if not new_points_added:
            return e_points, edge_vals

        if time.perf_counter() - start > timeout:
            if not redone_edge:
                dprint(
                    f"\nEDGE TIMEOUT RETRY: edge={edge_index}, eballs={eballs}, "
                    f"points={len(e_points)}, max_seg={max_seg:.6f}"
                )

                return build_edge(
                    locs,
                    rads,
                    vlocs,
                    res,
                    blocs,
                    brads,
                    eballs,
                    straight=straight,
                    vmid=vmid,
                    dnorm=-dnorm,
                    edub=edub,
                    edge_points1=edge_points1,
                    edge_verts=edge_verts,
                    redone_edge=True,
                    edge_index=edge_index,
                    debug=debug,
                    timeout=timeout
                )

            raise RuntimeError(
                f"Unable to build edge after retry. "
                f"edge_index={edge_index}, eballs={eballs}, edge_verts={edge_verts}"
            )