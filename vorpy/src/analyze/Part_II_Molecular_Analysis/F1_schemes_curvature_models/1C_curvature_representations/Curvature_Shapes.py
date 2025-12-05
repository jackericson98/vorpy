import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plot_saddle_surface():
    # Create a grid of x and y values
    x = np.linspace(-5, 5, 100)
    y = np.linspace(-5, 5, 100)
    x, y = np.meshgrid(x, y)

    # Define the saddle surface equation
    z = x ** 2 - y ** 2

    # Set up the figure and axis
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Plot the surface
    ax.plot_surface(x, y, z)

    # Labels and title
    ax.set_xlabel('X coordinate')
    ax.set_ylabel('Y coordinate')
    ax.set_zlabel('Z coordinate')
    # ax.set_title('Saddle Surface Plot')
    ax.axis('off')
    # Show the plot
    plt.show()


def plot_sphere(radius=3):
    # Define the sphere's angles for the parameterization
    phi = np.linspace(0, np.pi, 100)
    theta = np.linspace(0, 2 * np.pi, 100)

    # Create the grid
    phi, theta = np.meshgrid(phi, theta)

    # Spherical coordinates
    x = radius * np.sin(phi) * np.cos(theta)
    y = radius * np.sin(phi) * np.sin(theta)
    z = radius * np.cos(phi)

    # Set up the figure and axis
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')


    # Plot the surface
    ax.plot_surface(x, y, z)

    # Labels and title
    ax.set_xlabel('X coordinate')
    ax.set_ylabel('Y coordinate')
    ax.set_zlabel('Z coordinate')
    # ax.set_title('Sphere Plot')
    ax.axis('off')
    # Show the plot
    plt.show()


def plot_cylinder(radius=2, height=5):
    # Angles for the circular base
    theta = np.linspace(0, 2 * np.pi, 100)
    # Z coordinates from -height/2 to height/2
    z = np.linspace(-height / 2, height / 2, 100)

    # Create the grid
    theta, z = np.meshgrid(theta, z)

    # Cylinder coordinates
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)

    # Set up the figure and axis
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Plot the surface
    ax.plot_surface(x, y, z)

    # Labels and title
    ax.set_xlabel('X coordinate')
    ax.set_ylabel('Y coordinate')
    ax.set_zlabel('Z coordinate')
    # ax.set_title('Cylinder Plot')
    ax.axis('off')
    # Show the plot
    plt.show()


def plot_flat_surface():
    # Define x and y range
    x = np.linspace(-5, 5, 100)
    y = np.linspace(-5, 5, 100)

    # Create a meshgrid
    x, y = np.meshgrid(x, y)

    # Define a constant z value
    z = np.zeros_like(x)  # This makes the surface flat at z = 0

    # Set up the figure and axis
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Plot the surface
    ax.plot_surface(x, y, z)  # 'm' stands for magenta

    # Labels and title
    ax.set_xlabel('X coordinate')
    ax.set_ylabel('Y coordinate')
    ax.set_zlabel('Z coordinate')
    # ax.set_title('Flat Surface Plot')
    ax.axis('off')
    # Show the plot
    plt.show()


# Call the function to plot the saddle surface
# plot_saddle_surface()
plot_sphere()
# plot_cylinder()
plot_flat_surface()
