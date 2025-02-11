import os

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from System.system import System


def generate_cubic_lattice(n_points, spacing):
    """
    Generates coordinates for a simple cubic lattice.

    Parameters:
        n_points (int): Total number of points in the lattice.
        spacing (float): Distance between adjacent atoms.

    Returns:
        np.ndarray: Array of coordinates for each atom.
    """
    # Calculate the number of points along each axis
    points_per_axis = round(n_points ** (1 / 3))

    # Ensure that the total number of points matches the expected number
    # This step adjusts the number of points per axis if necessary
    while points_per_axis ** 3 < n_points:
        points_per_axis += 1

    # Generate grid of points
    x = np.linspace(0, spacing * (points_per_axis - 1), points_per_axis)
    y = np.linspace(0, spacing * (points_per_axis - 1), points_per_axis)
    z = np.linspace(0, spacing * (points_per_axis - 1), points_per_axis)

    # Create a meshgrid of points in 3D space
    xx, yy, zz = np.meshgrid(x, y, z)

    # Flatten the arrays and stack them into a (n_points, 3) array
    coordinates = np.vstack([xx.ravel(), yy.ravel(), zz.ravel()]).T

    # If the lattice is larger than requested due to cube root rounding, trim the excess
    return coordinates[:n_points]


# Example usage
n_points = 125  # Example: create 125 points
spacing = 1.5  # Spacing of 1.5 units between atoms

# Generate the lattice
lattice_coordinates = generate_cubic_lattice(n_points, spacing)

# # Print the coordinates
# print("Coordinates of atoms in the lattice:")
# print(lattice_coordinates)


def generate_cubic_lattice(n_points, spacing):
    """
    Generates coordinates for a simple cubic lattice.

    Parameters:
        n_points (int): Total number of points in the lattice.
        spacing (float): Distance between adjacent atoms.

    Returns:
        np.ndarray: Array of coordinates for each atom.
    """
    points_per_axis = round(n_points ** (1 / 3))
    while points_per_axis ** 3 < n_points:
        points_per_axis += 1

    x = np.linspace(0, spacing * (points_per_axis - 1), points_per_axis)
    y = np.linspace(0, spacing * (points_per_axis - 1), points_per_axis)
    z = np.linspace(0, spacing * (points_per_axis - 1), points_per_axis)

    xx, yy, zz = np.meshgrid(x, y, z)
    coordinates = np.vstack([xx.ravel(), yy.ravel(), zz.ravel()]).T
    return coordinates[:n_points]


def generate_2d_hexagonal_lattice(n, spacing):
    """
    Generates a 2D hexagonal lattice.

    Parameters:
        n (int): Approximate number of points (the actual number might be slightly different).
        spacing (float): Distance between adjacent points.

    Returns:
        np.ndarray: Coordinates of points in the hexagonal lattice.
    """
    # Approximate number of points along one axis
    num_per_row = int(np.sqrt(n))

    # Initialize list of coordinates
    coordinates = []

    # Adjust vertical spacing for hexagonal packing
    vertical_spacing = spacing * np.sqrt(3) / 2

    for row in range(num_per_row):
        for col in range(num_per_row):
            x = col * spacing + (row % 2) * (spacing / 2)
            y = row * vertical_spacing
            coordinates.append((x, y))

    return np.array(coordinates)


def generate_3d_hexagonal_lattice(layers, n_per_layer, spacing):
    """
    Generates a 3D hexagonal close-packed lattice by stacking 2D layers with offsets.

    Parameters:
        layers (int): Number of layers to generate.
        n_per_layer (int): Approximate number of points per layer.
        spacing (float): Distance between points in the same layer.

    Returns:
        np.ndarray: Coordinates of points in the 3D lattice.
    """
    # Generate one layer to start with
    layer = generate_2d_hexagonal_lattice(n_per_layer, spacing)

    # Stack layers
    coordinates = []
    for i in range(layers):
        # Determine z-coordinate based on layer index
        z = i * spacing * np.sqrt(6) / 3
        # Offset every other layer to achieve close-packing
        offset = (i % 2) * (spacing / 2)
        for (x, y) in layer:
            coordinates.append((x + offset, y, z))

    return np.array(coordinates)


def generate_concentric_spheres(center, num_layers, points_per_layer, radius_increment):
    """
    Generate points evenly distributed on concentric spheres around a center.

    Parameters:
        center (tuple): The center of the spheres (x, y, z).
        num_layers (int): Number of concentric spheres.
        points_per_layer (int): Number of points evenly distributed per sphere.
        radius_increment (float): Distance between each concentric sphere.

    Returns:
        list: A list of numpy arrays, each containing the points on one sphere.
    """
    layers = [[0, 0, 0]]
    golden_angle = np.pi * (3 - np.sqrt(5))  # This is approximately 2.39996323

    for i in range(1, num_layers + 1):
        radius = i * radius_increment
        points = []

        for j in range(points_per_layer):
            theta = golden_angle * j  # This spreads points out along the latitude
            z = radius - (2 * radius / points_per_layer) * j  # Linearly spaced along z-axis
            radius_at_z = np.sqrt(radius ** 2 - z ** 2)  # Radius of the horizontal slice at height z
            phi = np.arccos(z / radius)  # Angle from the vertical
            x = radius_at_z * np.cos(theta) + center[0]
            y = radius_at_z * np.sin(theta) + center[1]
            z = z + center[2]
            points.append([x, y, z])

        layers += points

    return np.array(layers)


# Example usage
center_point = (0, 0, 0)  # Center of the spheres
num_layers = 1  # Number of concentric spheres
points_per_layer = 100  # Points on each sphere
radius_increment = 1.0  # Spacing between spheres

# Generate the concentric spheres
spheres = generate_concentric_spheres(center_point, num_layers, points_per_layer, radius_increment)


def plot_3d_lattice(coordinates, colors=None):
    """
    Plots a 3D lattice of points.

    Parameters:
        coordinates (np.ndarray): Array of coordinates for each atom.
    """
    if colors is None:
        colors = 'blue'
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Scatter plot
    ax.scatter(coordinates[:, 0], coordinates[:, 1], coordinates[:, 2], color=colors, alpha=0.6, s=50)

    # Labeling axes
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_zlabel('Z Coordinate')

    # Set aspect of the plot to be equal
    ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio

    # Show the plot
    plt.title('3D Lattice Visualization')
    plt.show()


# Example usage
n_points = 125  # Create 125 points
spacing = 1.5  # Spacing of 1.5 units between atoms

# Generate the lattice
lattice_coordinates = generate_cubic_lattice(343, 5)

# Plot the lattice
plot_3d_lattice(spheres, ['red'] + ['blue' for _ in range(len(spheres) - 1)])
rads = [3] + [1 for _ in range(len(spheres) - 1)]

os.chdir('..')
# Create the
my_sys = System(atoms=[(spheres[i], rads[i]) for i in range(len(spheres))])

# Build the
my_sys.create_group(atoms=[0])

# Build the
my_sys.groups[0].build()


