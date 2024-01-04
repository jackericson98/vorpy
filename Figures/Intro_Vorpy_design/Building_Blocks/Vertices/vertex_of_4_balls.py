import numpy as np
from numpy import array as ar, sqrt as sqrt
import matplotlib.pyplot as plt
from Visualize.mpl_visualize import plot_atoms, plot_verts
from System.Network.verts.calc_vert import calc_vert

"""
Vertex Plotting: Set the 'vertex_type' Variable From the list below

    1. Equally Sized Balls - No Overlap
    2. Three Balls Equal One Different - No Overlap
    3. Two Balls Equal Two Different - No Overlap
    4. All Balls Different - No Overlap
    5. All Balls Equal - Two Overlap
    6. All Balls Different - Two Overlap
    7. All Balls Equal - Three Overlap
    8. All Balls Different - Three Overlap
    9. All Balls Equal - Positive Vertex, Four Overlap
    10. All Balls Different - Positive Vertex, Four Overlap
    11. All Balls Equal - Negative Vertex, Four Overlap
    12. All Balls Different - Negative Vertex, Four Overlap
"""

# Choose Vertex Type Below
vertex_type = 1

# Other settings
atom_alpha = 0.5
show_axes = False

# No Overlap Vertex 1
if vertex_type == 1:
    rads = [1.0, 1.0, 1.0, 1.0]
    locs = [3/np.sqrt(2), 0, -3/2], [-3/np.sqrt(2), 0, -3/2], [0, 3/np.sqrt(2), 3/2], [0, -3/np.sqrt(2), 3/2]
    title = 'Equally Sized Balls - No Overlap'

# No Overlap Vertex 2
elif vertex_type == 2:
    rads = [1.5, 1.5, 1.5, 2.5]
    locs = [3/np.sqrt(2), 0, -3/2], [-3/np.sqrt(2), 0, -3/2], [0, 3/np.sqrt(2), 3/2], [0, -3/np.sqrt(2), 3/2]
    title = 'Three Balls Equal One Ball Different - No Overlap'

# No Overlap Vertex 3
elif vertex_type == 3:
    rads = [1.5, 1.5, 1.75, 2.25]
    locs = [3/np.sqrt(2), 0, -3/2], [-3/np.sqrt(2), 0, -3/2], [0, 3/np.sqrt(2), 3/2], [0, -3/np.sqrt(2), 3/2]
    title = 'Two Balls Equal Two Different - No Overlap'

# No Overlap Vertex 4
elif vertex_type == 4:
    rads = [1.0, 1.25, 1.5, 2.0]
    locs = [3/np.sqrt(2), 0, -3/2], [-3/np.sqrt(2), 0, -3/2], [0, 3/np.sqrt(2), 3/2], [0, -3/np.sqrt(2), 3/2]
    title = 'All Balls Different - No Overlap'

# Two Overlap Vertex 1
elif vertex_type == 5:
    rads = [1.5, 1.5, 1.5, 1.5]
    locs = [1/np.sqrt(2), 0, -3/2], [-1/np.sqrt(2), 0, -3/2], [0, 3/np.sqrt(2), 3/2], [0, -3/np.sqrt(2), 3/2]
    title = 'All Balls Equal - Two Overlap'

# Two Overlap Vertex 2
elif vertex_type == 6:
    rads = [1.0, 1.5, 2.0, 2.5]
    locs = [3/np.sqrt(2), 0, -3/2], [-3/np.sqrt(2), 0, -3/2], [0, 3/np.sqrt(2), 3/2], [0, -3/np.sqrt(2), 3/2]
    title = 'All Balls Diffferent - Two Overlap'

# Three Overlap Vertex 1
elif vertex_type == 7:
    rads = [2.0, 2.0, 2.0, 2.0]
    locs = [1/np.sqrt(2), 0, -3/2], [-1/np.sqrt(2), 0, -3/2], [-1.5, 3/np.sqrt(2), 1/2], [0, -3/np.sqrt(2), 5/2]
    title = 'All Balls Equal - Three Overlap'


# Three Overlap Vertex 2
elif vertex_type == 8:
    rads = [2.5, 1.5, 2.0, 1.0]
    locs = [1/np.sqrt(2), 0, -3/2], [-3/np.sqrt(2), 0, -3/2], [-1.5, 3/np.sqrt(2), 1/2], [0, -3/np.sqrt(2), 5/2]
    title = 'All Balls Different - Three Overlap'

# All Overlap Vertex 1
elif vertex_type == 9:
    rads = [2.2, 2.2, 2.2, 2.2]
    locs = [3/np.sqrt(2), 0, -3/2], [-3/np.sqrt(2), 0, -3/2], [0, 3/np.sqrt(2), 3/2], [0, -3/np.sqrt(2), 3/2]
    title = 'All Balls Equal - All Overlapping, Positive'

# All Overlap Vertex 2
elif vertex_type == 10:
    rads = [1.0, 1.8, 2.0, 2.5]
    locs = [2/np.sqrt(2), 0, 1/2], [-3/np.sqrt(2), 0, -1/2], [0, 3/np.sqrt(2), 3/2], [0, -3/np.sqrt(2), 3/2]
    title = 'All Balls Different - Overlapping, Positive'

# All Overlap Vertex 3
elif vertex_type == 11:
    rads = [2.9, 2.9, 2.9, 2.9]
    locs = [3/np.sqrt(2), 0, -3/2], [-3/np.sqrt(2), 0, -3/2], [0, 3/np.sqrt(2), 3/2], [0, -3/np.sqrt(2), 3/2]
    title = 'All Balls Equal - All Overlapping, Negative'
    atom_alpha = 0.1

# All Overlap Vertex 4
elif vertex_type == 12:
    rads = [2.0, 2.5, 3.0, 3.5]
    locs = [1/np.sqrt(2), 0, -3/2], [-1/np.sqrt(2), 0, -3/2], [0, 3/np.sqrt(2), 3/2], [0, -3/np.sqrt(2), 3/2]
    title = 'All Balls Different - All Overlapping, Negative'
    atom_alpha = 0.1


# Calculate the vertex
my_vert = calc_vert(locs=ar([ar(_) for _ in locs]), rads=ar(rads))

# Make the plot
fig = plt.figure()
ax = fig.add_subplot(projection='3d')


# Plot the atoms
plot_atoms(locs, rads, fig=fig, ax=ax, res=10, alpha=atom_alpha)
# Plot the vertices
plot_verts([my_vert[0]], [abs(my_vert[1])], fig=fig, ax=ax, spheres=True, res=10)


# Set the axes lines
if show_axes:
    ax.plot([-3, -2], [-3, -3], [0, 0])
    ax.plot([-3, -3], [-2, -3], [0, 0])
    ax.plot([-3, -3], [-3, -3], [0, 1])

    # Set the axes labels
    ax.text(x=-2, y=-3, z=-0.25, s='x')
    ax.text(x=-3, y=-2, z=-0.25, s='y')
    ax.text(x=-3, y=-3, z=1, s='z')

# Set the scales for the figure
ax.set_xlim(-5, 5)
ax.set_ylim(-5, 5)
ax.set_zlim(-5, 5)

# Set the title
ax.set_title(title, font=dict(size=20, family='serif'))

# Show the plot
plt.show()
