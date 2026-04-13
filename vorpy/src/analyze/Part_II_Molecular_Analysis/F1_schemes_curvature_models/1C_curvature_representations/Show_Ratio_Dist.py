import numpy as np
import matplotlib.pyplot as plt


def plot_circle(rad, loc, color='k', linewidth=2, show_axes=False):
    x0, y0 = loc

    xs = np.linspace(x0 - rad, x0 + rad, 500)

    ys_upper = y0 + np.sqrt(rad**2 - (xs - x0)**2)
    ys_lower = y0 - np.sqrt(rad**2 - (xs - x0)**2)

    plt.plot(xs, ys_upper, c=color, linewidth=linewidth)
    plt.plot(xs, ys_lower, c=color, linewidth=linewidth)

    plt.gca().set_aspect('equal', adjustable='box')
    if not show_axes:
        plt.axis('off')


def plot_circle_rads(rad,
                     loc,
                     angle_deg=0,
                     color='k',
                     linewidth=2,
                     label=None,
                     label_offset=0.05,
                     fontsize=12):

    x0, y0 = loc
    theta = np.deg2rad(angle_deg)

    # Endpoint of radius
    x1 = x0 + rad * np.cos(theta)
    y1 = y0 + rad * np.sin(theta)

    # Draw radius line
    plt.plot([x0, x1], [y0, y1], c=color, linewidth=linewidth)

    # -------- LABEL --------
    if label is not None:
        # Midpoint of radius
        xm = x0 + 0.5 * rad * np.cos(theta)
        ym = y0 + 0.5 * rad * np.sin(theta)

        # Perpendicular offset for readability
        perp = np.array([-np.sin(theta), np.cos(theta)])
        xm += label_offset * perp[0]
        ym += label_offset * perp[1]

        plt.text(xm, ym,
                 label,
                 ha='center',
                 va='center',
                 fontsize=fontsize,
                 color=color)


def plot_line_between(circle1, circle2,
                      color='k',
                      linewidth=2,
                      tab_height=0.1,
                      vertical_offset=1,
                      fontsize=12):

    c1 = np.array(circle1[0], dtype=float)
    r1 = circle1[1]

    c2 = np.array(circle2[0], dtype=float)
    r2 = circle2[1]

    # Direction from circle1 → circle2
    direction = c2 - c1
    dist = np.linalg.norm(direction)

    if dist == 0:
        raise ValueError("Circle centers are identical; direction undefined.")

    drhat = direction / dist

    # Surface points
    start = c1 + r1 * drhat
    end = c2 - r2 * drhat

    # Apply vertical offset (to everything consistently)
    start_off = start.copy()
    end_off = end.copy()
    start_off[1] -= vertical_offset
    end_off[1] -= vertical_offset

    # Main line
    plt.plot([start_off[0], end_off[0]],
             [start_off[1], end_off[1]],
             c=color, linewidth=linewidth)

    # Perpendicular unit vector
    perp = np.array([-drhat[1], drhat[0]])
    h = tab_height / 2

    # Tabs
    s1 = start_off - h * perp
    s2 = start_off + h * perp
    e1 = end_off - h * perp
    e2 = end_off + h * perp

    plt.plot([s1[0], s2[0]], [s1[1], s2[1]],
             c=color, linewidth=linewidth)
    plt.plot([e1[0], e2[0]], [e1[1], e2[1]],
             c=color, linewidth=linewidth)

    # -------- DISTANCE LABEL --------
    mid = 0.5 * (start_off + end_off)
    d_surface = dist - r1 - r2

    plt.text(mid[0], mid[1] + 0.15,
             f"d = {d_surface:.2f}",
             ha='center',
             va='bottom',
             fontsize=fontsize,
             color=color)

    return start, end


def main():

    CIRCLE1 = [-1.0, 0], 0.5
    CIRCLE2 = [1.0, 0], 1

    CIRCLE1_COLOR = 'k'
    CIRCLE2_COLOR = 'k'

    LINEWIDTH = 3

    plot_circle(CIRCLE1[1], CIRCLE1[0], color=CIRCLE1_COLOR, linewidth=LINEWIDTH)
    plot_circle(CIRCLE2[1], CIRCLE2[0], color=CIRCLE2_COLOR, linewidth=LINEWIDTH)
    plot_circle_rads(CIRCLE1[1], CIRCLE1[0],
                     angle_deg=0,
                     label="r₁",
                     color=CIRCLE1_COLOR,
                     linewidth=LINEWIDTH)

    plot_circle_rads(CIRCLE2[1], CIRCLE2[0],
                     angle_deg=180,
                     label="r₂",
                     color=CIRCLE2_COLOR,
                     linewidth=LINEWIDTH)
    plot_line_between(CIRCLE1, CIRCLE2, color='k', linewidth=2, tab_height=0.1)
    plt.show()


if __name__ == '__main__':
    main()