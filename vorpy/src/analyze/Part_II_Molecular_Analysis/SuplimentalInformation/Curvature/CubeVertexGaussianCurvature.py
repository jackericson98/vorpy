"""
CubeVertexGaussianCurvature.py

Generate three separate figures illustrating Gaussian curvature
as angle defect at a cube vertex.

Outputs:
    cube.png
    cube.svg

    cube_vertex.png
    cube_vertex.svg

    flat_vertex.png
    flat_vertex.svg

Figure 1:
    Cube with a highlighted and labeled selected vertex.

Figure 2:
    3D cube vertex with three face-interior angles of pi/2.

Figure 3:
    Flat vertex with three equal face angles of 2*pi/3.
"""

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.patches import Arc, Polygon, Circle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


# ============================================================
# Settings
# ============================================================

DPI = 300

plt.rcParams.update({
    "font.size": 14,
    "mathtext.fontset": "stix",
    "font.family": "DejaVu Sans",
})


# ============================================================
# Helper functions
# ============================================================

def save_figure(fig, basename):
    """Save a figure as both PNG and SVG."""

    fig.savefig(
        f"{basename}.png",
        dpi=DPI,
        bbox_inches="tight",
        transparent=True,
        pad_inches=0.15,
    )

    fig.savefig(
        f"{basename}.svg",
        bbox_inches="tight",
        transparent=True,
        pad_inches=0.15,
    )


def draw_arc(
    ax,
    center,
    radius,
    theta1,
    theta2,
    label,
    label_radius=None,
    fontsize=18,
):
    """Draw a 2D angle arc and label it."""

    cx, cy = center

    arc = Arc(
        (cx, cy),
        2 * radius,
        2 * radius,
        angle=0,
        theta1=theta1,
        theta2=theta2,
        linewidth=2,
    )

    ax.add_patch(arc)

    if label_radius is None:
        label_radius = radius * 1.35

    theta_mid = np.deg2rad(
        (theta1 + theta2) / 2
    )

    x = cx + label_radius * np.cos(theta_mid)
    y = cy + label_radius * np.sin(theta_mid)

    ax.text(
        x,
        y,
        label,
        ha="center",
        va="center",
        fontsize=fontsize,
    )


# ============================================================
# Figure 1 — Cube
# ============================================================

def draw_cube(ax):
    """
    Draw a projected cube with one highlighted vertex.
    """

    # Front square
    A = np.array([0.0, 0.0])
    B = np.array([2.2, 0.0])
    C = np.array([2.2, 2.2])
    D = np.array([0.0, 2.2])

    # Projection offset for back square
    offset = np.array([0.85, 0.65])

    A2 = A + offset
    B2 = B + offset
    C2 = C + offset
    D2 = D + offset

    # --------------------------------------------------------
    # Visible faces
    # --------------------------------------------------------

    front = Polygon(
        [A, B, C, D],
        closed=True,
        alpha=0.10,
    )

    top = Polygon(
        [D, C, C2, D2],
        closed=True,
        alpha=0.16,
    )

    side = Polygon(
        [B, B2, C2, C],
        closed=True,
        alpha=0.22,
    )

    ax.add_patch(front)
    ax.add_patch(top)
    ax.add_patch(side)

    # --------------------------------------------------------
    # Cube edges
    # --------------------------------------------------------

    edges = [
        (A, B),
        (B, C),
        (C, D),
        (D, A),

        (A2, B2),
        (B2, C2),
        (C2, D2),
        (D2, A2),

        (A, A2),
        (B, B2),
        (C, C2),
        (D, D2),
    ]

    for p1, p2 in edges:
        ax.plot(
            [p1[0], p2[0]],
            [p1[1], p2[1]],
            linewidth=2,
        )

    # --------------------------------------------------------
    # Selected vertex
    # --------------------------------------------------------

    ax.scatter(
        [D[0]],
        [D[1]],
        s=120,
        zorder=10,
    )

    ax.annotate(
        "selected vertex",
        xy=D,
        xytext=(-0.75, 3.25),
        arrowprops=dict(
            arrowstyle="->",
            linewidth=1.5,
        ),
        fontsize=14,
        ha="left",
    )

    # --------------------------------------------------------
    # Appearance
    # --------------------------------------------------------

    ax.set_title(
        "Cube",
        fontsize=22,
        weight="bold",
        pad=12,
    )

    ax.set_xlim(-1.0, 3.7)
    ax.set_ylim(-0.5, 3.7)

    ax.set_aspect("equal")
    ax.axis("off")


# ============================================================
# Figure 2 — 3D Cube Vertex
# ============================================================

def draw_cube_vertex(ax):
    """
    Draw a true 3D cube corner.

    Three mutually perpendicular square faces meet at the origin.

    Each face contributes an interior angle of pi/2 at the vertex.

    IMPORTANT:
    The labeled angles are face-interior angles, not dihedral angles.
    """

    O = np.array([0.0, 0.0, 0.0])

    L = 2.0

    # --------------------------------------------------------
    # Three mutually perpendicular faces
    # --------------------------------------------------------

    # xy face: z = 0
    face_xy = np.array([
        O,
        [L, 0, 0],
        [L, L, 0],
        [0, L, 0],
    ])

    # xz face: y = 0
    face_xz = np.array([
        O,
        [L, 0, 0],
        [L, 0, L],
        [0, 0, L],
    ])

    # yz face: x = 0
    face_yz = np.array([
        O,
        [0, L, 0],
        [0, L, L],
        [0, 0, L],
    ])

    faces = [
        face_xy,
        face_xz,
        face_yz,
    ]

    poly = Poly3DCollection(
        faces,
        alpha=0.15,
        edgecolor="black",
        linewidth=1.8,
    )

    ax.add_collection3d(poly)

    # --------------------------------------------------------
    # Three edges incident to the vertex
    # --------------------------------------------------------

    edge_vectors = [
        np.array([L, 0, 0]),
        np.array([0, L, 0]),
        np.array([0, 0, L]),
    ]

    for endpoint in edge_vectors:
        ax.plot(
            [O[0], endpoint[0]],
            [O[1], endpoint[1]],
            [O[2], endpoint[2]],
            linewidth=3,
        )

    # --------------------------------------------------------
    # Highlight the vertex
    # --------------------------------------------------------

    ax.scatter(
        [0],
        [0],
        [0],
        s=110,
        zorder=20,
    )

    # --------------------------------------------------------
    # Draw pi/2 angle arcs
    # --------------------------------------------------------

    r = 0.65

    theta = np.linspace(
        0,
        np.pi / 2,
        100,
    )

    # xy face
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    z = np.zeros_like(theta)

    ax.plot(
        x,
        y,
        z,
        linewidth=2.5,
    )

    # xz face
    x = r * np.cos(theta)
    y = np.zeros_like(theta)
    z = r * np.sin(theta)

    ax.plot(
        x,
        y,
        z,
        linewidth=2.5,
    )

    # yz face
    x = np.zeros_like(theta)
    y = r * np.cos(theta)
    z = r * np.sin(theta)

    ax.plot(
        x,
        y,
        z,
        linewidth=2.5,
    )

    # --------------------------------------------------------
    # pi/2 labels
    # --------------------------------------------------------

    label_r = 0.95

    q = np.sqrt(2) / 2

    # xy face
    ax.text(
        label_r * q,
        label_r * q,
        0,
        r"$\frac{\pi}{2}$",
        fontsize=20,
        ha="center",
        va="center",
    )

    # xz face
    ax.text(
        label_r * q,
        0,
        label_r * q,
        r"$\frac{\pi}{2}$",
        fontsize=20,
        ha="center",
        va="center",
    )

    # yz face
    ax.text(
        0,
        label_r * q,
        label_r * q,
        r"$\frac{\pi}{2}$",
        fontsize=20,
        ha="center",
        va="center",
    )

    # --------------------------------------------------------
    # Camera
    # --------------------------------------------------------

    ax.set_xlim(0, 2.3)
    ax.set_ylim(0, 2.3)
    ax.set_zlim(0, 2.3)

    ax.set_box_aspect(
        (1, 1, 1)
    )

    ax.view_init(
        elev=25,
        azim=35,
    )

    ax.set_axis_off()

    # --------------------------------------------------------
    # Title
    # --------------------------------------------------------

    ax.text2D(
        0.50,
        0.96,
        "Cube vertex",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=22,
        weight="bold",
    )


# ============================================================
# Figure 3 — Flat Vertex
# ============================================================

def draw_flat_vertex(ax):
    """
    Draw a flat 2*pi neighborhood split into three equal sectors.

    Each sector is:

        2*pi / 3 = 120 degrees
    """

    O = np.array([0.0, 0.0])

    R = 2.0

    angles = [
        0,
        120,
        240,
    ]

    # --------------------------------------------------------
    # Outer circle
    # --------------------------------------------------------

    circle = Circle(
        O,
        R,
        fill=False,
        linewidth=2,
    )

    ax.add_patch(circle)

    # --------------------------------------------------------
    # Three radial boundaries
    # --------------------------------------------------------

    for theta in angles:

        t = np.deg2rad(theta)

        endpoint = O + R * np.array([
            np.cos(t),
            np.sin(t),
        ])

        ax.plot(
            [O[0], endpoint[0]],
            [O[1], endpoint[1]],
            linewidth=2,
        )

    # --------------------------------------------------------
    # Three sectors
    # --------------------------------------------------------

    sector_ranges = [
        (0, 120),
        (120, 240),
        (240, 360),
    ]

    for theta1, theta2 in sector_ranges:

        ts = np.linspace(
            np.deg2rad(theta1),
            np.deg2rad(theta2),
            50,
        )

        pts = np.column_stack((
            R * np.cos(ts),
            R * np.sin(ts),
        ))

        pts = np.vstack((
            O,
            pts,
            O,
        ))

        sector = Polygon(
            pts,
            closed=True,
            alpha=0.08,
        )

        ax.add_patch(sector)

        # Angle arc and label
        draw_arc(
            ax=ax,
            center=O,
            radius=0.65,
            theta1=theta1,
            theta2=theta2,
            label=r"$\frac{2\pi}{3}$",
            label_radius=1.15,
            fontsize=20,
        )

    # --------------------------------------------------------
    # Central vertex
    # --------------------------------------------------------

    ax.scatter(
        [0],
        [0],
        s=90,
        zorder=10,
    )

    # --------------------------------------------------------
    # Appearance
    # --------------------------------------------------------

    ax.set_title(
        "Flat vertex",
        fontsize=22,
        weight="bold",
        pad=12,
    )

    ax.set_xlim(-2.4, 2.4)
    ax.set_ylim(-2.4, 2.4)

    ax.set_aspect("equal")

    ax.axis("off")


# ============================================================
# Generate Figure 1 — Cube
# ============================================================

fig1, ax1 = plt.subplots(
    figsize=(6, 6)
)

draw_cube(ax1)

save_figure(
    fig1,
    "cube",
)


# ============================================================
# Generate Figure 2 — Cube Vertex
# ============================================================

fig2 = plt.figure(
    figsize=(6, 6)
)

ax2 = fig2.add_subplot(
    111,
    projection="3d",
)

draw_cube_vertex(ax2)

save_figure(
    fig2,
    "cube_vertex",
)


# ============================================================
# Generate Figure 3 — Flat Vertex
# ============================================================

fig3, ax3 = plt.subplots(
    figsize=(6, 6)
)

draw_flat_vertex(ax3)

save_figure(
    fig3,
    "flat_vertex",
)


# ============================================================
# Display figures
# ============================================================

plt.show()