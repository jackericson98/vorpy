import numpy as np
import matplotlib.pyplot as plt


# ----------------------------- #
# Global style / appearance
# ----------------------------- #

FRONT_COLOR = 'white'
BACK_COLOR = 'white'

LINE1_COLOR = '#8b0000'
LINE2_COLOR = '#003b8e'

SURFACE_ALPHA_FRONT = 0.92
SURFACE_ALPHA_BACK = 0.92

LINE_WIDTH = 4
POINT_SIZE = 90
POINT_COLOR = 'black'

SHELL_THICKNESS = 0.015
LINE_OFFSET = 0.03

CYLINDER_CAP_COLOR = 'white'
CYLINDER_CAP_ALPHA = 1.0

POINT_OFFSET = 0.06

EDGE_COLOR = 'black'
EDGE_WIDTH = 5


def offset_point(x, y, z, nx, ny, nz, eps=POINT_OFFSET):
    return x + eps * nx, y + eps * ny, z + eps * nz


def normalize_vectors(nx, ny, nz):
    mag = np.sqrt(nx ** 2 + ny ** 2 + nz ** 2)
    mag = np.where(mag == 0, 1.0, mag)

    return nx / mag, ny / mag, nz / mag


def set_axes_equal(ax):
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    y_range = abs(y_limits[1] - y_limits[0])
    z_range = abs(z_limits[1] - z_limits[0])

    x_middle = np.mean(x_limits)
    y_middle = np.mean(y_limits)
    z_middle = np.mean(z_limits)

    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

    try:
        ax.set_box_aspect((1, 1, 1))
    except Exception:
        pass


def style_axis(ax, elev=24, azim=-55):
    ax.axis('off')
    ax.view_init(elev=elev, azim=azim)
    set_axes_equal(ax)


def create_figure(figsize=(8, 8)):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')

    return fig, ax


def offset_curve(x, y, z, nx, ny, nz, eps=LINE_OFFSET):
    return x + eps * nx, y + eps * ny, z + eps * nz


def plot_double_sided_surface(ax, x, y, z, nx, ny, nz, thickness=SHELL_THICKNESS):
    x_front = x + thickness * nx
    y_front = y + thickness * ny
    z_front = z + thickness * nz

    x_back = x - thickness * nx
    y_back = y - thickness * ny
    z_back = z - thickness * nz

    ax.plot_surface(
        x_back,
        y_back,
        z_back,
        color=BACK_COLOR,
        alpha=SURFACE_ALPHA_BACK,
        linewidth=0,
        edgecolor='none',
        antialiased=True,
        shade=True
    )

    ax.plot_surface(
        x_front,
        y_front,
        z_front,
        color=FRONT_COLOR,
        alpha=SURFACE_ALPHA_FRONT,
        linewidth=0,
        edgecolor='none',
        antialiased=True,
        shade=True
    )


def plot_principal_line(ax, x, y, z, color, lw=LINE_WIDTH):
    ax.plot(
        x,
        y,
        z,
        color=color,
        lw=lw,
        zorder=100
    )


def plot_intersection_point(ax, x, y, z):
    ax.scatter(
        [x], [y], [z],
        s=160,
        color='white',
        edgecolors='black',
        linewidths=2,
        depthshade=False
    )


def plot_boundary_curve(ax, x, y, z, color=EDGE_COLOR, lw=EDGE_WIDTH):
    ax.plot(x, y, z, color=color, lw=lw)


def plot_flat_surface(size=3.6, n=120):
    x = np.linspace(-size, size, n)
    y = np.linspace(-size, size, n)
    x, y = np.meshgrid(x, y)
    z = np.zeros_like(x)

    nx = np.zeros_like(x)
    ny = np.zeros_like(y)
    nz = np.ones_like(z)

    fig, ax = create_figure()

    plot_double_sided_surface(ax, x, y, z, nx, ny, nz)

    # outer boundary of plane
    eps = SHELL_THICKNESS
    plot_boundary_curve(ax, np.linspace(-size, size, 300), np.full(300, -size), np.full(300, eps))
    plot_boundary_curve(ax, np.linspace(-size, size, 300), np.full(300, size), np.full(300, eps))
    plot_boundary_curve(ax, np.full(300, -size), np.linspace(-size, size, 300), np.full(300, eps))
    plot_boundary_curve(ax, np.full(300, size), np.linspace(-size, size, 300), np.full(300, eps))

    t = np.linspace(-size, size, 400)

    x_line1 = t
    y_line1 = np.zeros_like(t)
    z_line1 = np.zeros_like(t)

    nx_line1 = np.zeros_like(t)
    ny_line1 = np.zeros_like(t)
    nz_line1 = np.ones_like(t)

    x_line1, y_line1, z_line1 = offset_curve(
        x_line1, y_line1, z_line1,
        nx_line1, ny_line1, nz_line1
    )

    x_line2 = np.zeros_like(t)
    y_line2 = t
    z_line2 = np.zeros_like(t)

    nx_line2 = np.zeros_like(t)
    ny_line2 = np.zeros_like(t)
    nz_line2 = np.ones_like(t)

    x_line2, y_line2, z_line2 = offset_curve(
        x_line2, y_line2, z_line2,
        nx_line2, ny_line2, nz_line2
    )

    plot_principal_line(ax, x_line1, y_line1, z_line1, LINE1_COLOR)
    plot_principal_line(ax, x_line2, y_line2, z_line2, LINE2_COLOR)

    point_x, point_y, point_z = offset_point(
        np.array([0.0]),
        np.array([0.0]),
        np.array([0.0]),
        np.array([0.0]),
        np.array([0.0]),
        np.array([1.0])
    )
    plot_intersection_point(ax, point_x[0], point_y[0], point_z[0])

    style_axis(ax, elev=25, azim=-45)
    plt.show()


def plot_sphere(radius=3.0, n_phi=120, n_theta=180):
    phi = np.linspace(0, np.pi, n_phi)
    theta = np.linspace(0, 2 * np.pi, n_theta)
    phi, theta = np.meshgrid(phi, theta)

    x = radius * np.sin(phi) * np.cos(theta)
    y = radius * np.sin(phi) * np.sin(theta)
    z = radius * np.cos(phi)

    nx = x / radius
    ny = y / radius
    nz = z / radius

    fig, ax = create_figure()

    plot_double_sided_surface(ax, x, y, z, nx, ny, nz)

    theta_line = np.linspace(0, 2 * np.pi, 500)
    phi0 = np.pi / 2

    x_line1 = radius * np.sin(phi0) * np.cos(theta_line)
    y_line1 = radius * np.sin(phi0) * np.sin(theta_line)
    z_line1 = radius * np.cos(phi0) * np.ones_like(theta_line)

    nx_line1 = x_line1 / radius
    ny_line1 = y_line1 / radius
    nz_line1 = z_line1 / radius

    x_line1, y_line1, z_line1 = offset_curve(
        x_line1, y_line1, z_line1,
        nx_line1, ny_line1, nz_line1
    )

    phi_line = np.linspace(0.12, np.pi - 0.12, 400)
    theta0 = np.pi / 5 + np.pi

    x_line2 = radius * np.sin(phi_line) * np.cos(theta0)
    y_line2 = radius * np.sin(phi_line) * np.sin(theta0)
    z_line2 = radius * np.cos(phi_line)

    nx_line2 = x_line2 / radius
    ny_line2 = y_line2 / radius
    nz_line2 = z_line2 / radius

    x_line2, y_line2, z_line2 = offset_curve(
        x_line2, y_line2, z_line2,
        nx_line2, ny_line2, nz_line2
    )

    plot_principal_line(ax, x_line1, y_line1, z_line1, LINE1_COLOR)
    plot_principal_line(ax, x_line2, y_line2, z_line2, LINE2_COLOR)

    x0 = radius * np.cos(theta0)
    y0 = radius * np.sin(theta0)
    z0 = 0.0

    nx0 = x0 / radius
    ny0 = y0 / radius
    nz0 = 0.0

    point_x, point_y, point_z = offset_point(
        np.array([x0]),
        np.array([y0]),
        np.array([z0]),
        np.array([nx0]),
        np.array([ny0]),
        np.array([nz0])
    )
    plot_intersection_point(ax, point_x[0], point_y[0], point_z[0])

    style_axis(ax, elev=22, azim=-52)
    plt.show()


def plot_cylinder_caps(ax, radius, height, n_r=80, n_theta=220):
    r = np.linspace(0, radius, n_r)
    theta = np.linspace(0, 2 * np.pi, n_theta)
    r, theta = np.meshgrid(r, theta)

    x = r * np.cos(theta)
    y = r * np.sin(theta)

    z_top = np.full_like(x, height / 2)
    z_bottom = np.full_like(x, -height / 2)

    # filled top disk
    ax.plot_surface(
        x, y, z_top,
        color='white',
        alpha=1.0,
        linewidth=0,
        edgecolor='none',
        shade=False
    )

    # filled bottom disk
    ax.plot_surface(
        x, y, z_bottom,
        color='white',
        alpha=1.0,
        linewidth=0,
        edgecolor='none',
        shade=False
    )

    # top and bottom rim outlines
    theta_rim = np.linspace(0, 2 * np.pi, 500)
    x_rim = radius * np.cos(theta_rim)
    y_rim = radius * np.sin(theta_rim)

    ax.plot(x_rim, y_rim, np.full_like(theta_rim, height / 2), color='black', lw=2.5)
    ax.plot(x_rim, y_rim, np.full_like(theta_rim, -height / 2), color='black', lw=2.5)


def plot_cylinder(radius=2.0, height=6.0, n_theta=180, n_z=120):
    theta = np.linspace(0, 2 * np.pi, n_theta)
    z = np.linspace(-height / 2, height / 2, n_z)
    theta, z = np.meshgrid(theta, z)

    x = radius * np.cos(theta)
    y = radius * np.sin(theta)

    nx = np.cos(theta)
    ny = np.sin(theta)
    nz = np.zeros_like(z)

    fig, ax = create_figure()

    plot_double_sided_surface(ax, x, y, z, nx, ny, nz)
    plot_cylinder_caps(ax, radius=radius, height=height)

    z_line1 = np.linspace(-height / 2, height / 2, 400)
    theta0 = np.pi / 5 + np.pi

    x_line1 = radius * np.cos(theta0) * np.ones_like(z_line1)
    y_line1 = radius * np.sin(theta0) * np.ones_like(z_line1)

    nx_line1 = np.cos(theta0) * np.ones_like(z_line1)
    ny_line1 = np.sin(theta0) * np.ones_like(z_line1)
    nz_line1 = np.zeros_like(z_line1)

    x_line1, y_line1, z_line1 = offset_curve(
        x_line1, y_line1, z_line1,
        nx_line1, ny_line1, nz_line1
    )

    theta_line2 = np.linspace(0, 2 * np.pi, 500)
    z0 = np.zeros_like(theta_line2)

    x_line2 = radius * np.cos(theta_line2)
    y_line2 = radius * np.sin(theta_line2)

    nx_line2 = np.cos(theta_line2)
    ny_line2 = np.sin(theta_line2)
    nz_line2 = np.zeros_like(theta_line2)

    x_line2, y_line2, z_line2 = offset_curve(
        x_line2, y_line2, z0,
        nx_line2, ny_line2, nz_line2
    )

    plot_principal_line(ax, x_line1, y_line1, z_line1, LINE1_COLOR)
    plot_principal_line(ax, x_line2, y_line2, z_line2, LINE2_COLOR)

    x0 = radius * np.cos(theta0)
    y0 = radius * np.sin(theta0)
    z0 = 0.0

    nx0 = np.cos(theta0)
    ny0 = np.sin(theta0)
    nz0 = 0.0

    point_x, point_y, point_z = offset_point(
        np.array([x0]),
        np.array([y0]),
        np.array([z0]),
        np.array([nx0]),
        np.array([ny0]),
        np.array([nz0])
    )
    plot_intersection_point(ax, point_x[0], point_y[0], point_z[0])

    style_axis(ax, elev=18, azim=-55)
    plt.show()


def plot_saddle_surface(size=2.0, a=0.7, n=160):
    x = np.linspace(-size, size, n)
    y = np.linspace(-size, size, n)
    x, y = np.meshgrid(x, y)

    z = a * (x ** 2 - y ** 2)

    fx = 2 * a * x
    fy = -2 * a * y

    nx = -fx
    ny = -fy
    nz = np.ones_like(z)

    nx, ny, nz = normalize_vectors(nx, ny, nz)

    fig, ax = create_figure()

    plot_double_sided_surface(ax, x, y, z, nx, ny, nz)

    # saddle outer boundary
    t = np.linspace(-size, size, 400)

    xb = t
    yb = np.full_like(t, -size)
    zb = a * (xb ** 2 - yb ** 2)
    plot_boundary_curve(ax, xb, yb, zb)

    xb = t
    yb = np.full_like(t, size)
    zb = a * (xb ** 2 - yb ** 2)
    plot_boundary_curve(ax, xb, yb, zb)

    yb = t
    xb = np.full_like(t, -size)
    zb = a * (xb ** 2 - yb ** 2)
    plot_boundary_curve(ax, xb, yb, zb)

    yb = t
    xb = np.full_like(t, size)
    zb = a * (xb ** 2 - yb ** 2)
    plot_boundary_curve(ax, xb, yb, zb)

    x_line1 = t
    y_line1 = np.zeros_like(t)
    z_line1 = a * (t ** 2)

    fx_line1 = 2 * a * x_line1
    fy_line1 = np.zeros_like(t)

    nx_line1 = -fx_line1
    ny_line1 = -fy_line1
    nz_line1 = np.ones_like(t)

    nx_line1, ny_line1, nz_line1 = normalize_vectors(nx_line1, ny_line1, nz_line1)

    x_line1, y_line1, z_line1 = offset_curve(
        x_line1, y_line1, z_line1,
        nx_line1, ny_line1, nz_line1
    )

    x_line2 = np.zeros_like(t)
    y_line2 = t
    z_line2 = -a * (t ** 2)

    fx_line2 = np.zeros_like(t)
    fy_line2 = -2 * a * y_line2

    nx_line2 = -fx_line2
    ny_line2 = -fy_line2
    nz_line2 = np.ones_like(t)

    nx_line2, ny_line2, nz_line2 = normalize_vectors(nx_line2, ny_line2, nz_line2)

    x_line2, y_line2, z_line2 = offset_curve(
        x_line2, y_line2, z_line2,
        nx_line2, ny_line2, nz_line2
    )

    plot_principal_line(ax, x_line1, y_line1, z_line1, LINE1_COLOR)
    plot_principal_line(ax, x_line2, y_line2, z_line2, LINE2_COLOR)

    nx0, ny0, nz0 = normalize_vectors(
        np.array([0.0]),
        np.array([0.0]),
        np.array([1.0])
    )
    point_x, point_y, point_z = offset_point(
        np.array([0.0]),
        np.array([0.0]),
        np.array([0.0]),
        nx0,
        ny0,
        nz0
    )
    plot_intersection_point(ax, point_x[0], point_y[0], point_z[0])

    style_axis(ax, elev=28, azim=130)
    plt.show()


def main():
    plot_sphere()
    plot_flat_surface()
    plot_cylinder()
    plot_saddle_surface()


if __name__ == "__main__":
    main()