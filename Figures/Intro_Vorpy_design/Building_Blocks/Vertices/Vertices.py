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
    6. All Balls Equal - Three Overlap
    7. All Balls Equal - All Overlapping, Positive Vertex
    8. All Balls Equal - All Overlapping, Negative Vertex
    9. 

"""

# Choose Vertex Type Below
vertex_type = 8

# No Overlap Vertex 1
if vertex_type == 1:
    rads = [1.0, 1.0, 1.0, 1.0]
    locs = [1.0, 0.0, -1.0], [-sqrt(2), -sqrt(2), -1.0], [-sqrt(2), sqrt(2), -1.0], [-sqrt(2)/2, 0.0, 1.0]
    title = 'Equally Sized Balls - No Overlap'

# No Overlap Vertex 2
elif vertex_type == 2:
    rads = [0.5, 1.0, 1.0, 1.0]
    locs = [1.0, 0.0, -1.0], [-sqrt(2), -sqrt(2), -1.0], [-sqrt(2), sqrt(2), -1.0], [-sqrt(2) / 2, 0.0, 1.0]
    title = 'Three Balls Equal One Ball Different - No Overlap'

# No Overlap Vertex 3
elif vertex_type == 3:
    rads = [0.5, 0.75, 1.0, 1.0]
    locs = [1.0, 0.0, -1.0], [-sqrt(2), -sqrt(2), -1.0], [-sqrt(2), sqrt(2), -1.0], [-sqrt(2) / 2, 0.0, 1.0]
    title = 'Two Balls Equal Two Different - No Overlap'

# No Overlap Vertex 4
elif vertex_type == 4:
    rads = [0.5, 0.75, 1.0, 1.25]
    locs = [1.0, 0.0, -1.0], [-sqrt(2), -sqrt(2), -1.0], [-sqrt(2), sqrt(2), -1.0], [-sqrt(2) / 2, 0.0, 1.0]
    title = 'All Balls Different - No Overlap'

# Two Overlap Vertex 1
elif vertex_type == 5:
    rads = [1.0, 1.0, 1.0, 1.0]
    locs = [0.2, 0.5, -1.0], [-sqrt(2), -sqrt(2), -1.0], [-sqrt(2), sqrt(2), -1.0], [-sqrt(2) / 2, 0.0, 1.0]
    title = 'All Balls Equal - Two Overlap'

# Three Overlap Vertex 1
elif vertex_type == 6:
    rads = [1.0, 1.0, 1.0, 1.0]
    locs = [0.5, 0.5, -1.0], [-0.5, -0.5, -1.0], [-0.5, 0.75, -1.0], [0.0, 0.2, 1.5]
    title = 'All Balls Equal - Three Overlap'

# All Overlap Vertex 1
elif vertex_type == 7:
    rads = [1.5, 1.5, 1.5, 1.5]
    locs = [1.0, 0.0, -1.0], [-sqrt(2), -sqrt(2), -1.0], [-sqrt(2), sqrt(2), -1.0], [-sqrt(2)/2, 0.0, 1.0]
    title = 'All Balls Equal - All Overlapping, Positive Vertex'

# All Overlap Vertex 2
elif vertex_type == 8:
    rads = [2.2, 2.2, 2.2, 2.2]
    locs = [1.0, 0.0, -1.0], [-sqrt(2), -sqrt(2), -1.0], [-sqrt(2), sqrt(2), -1.0], [-sqrt(2)/2, 0.0, 1.0]
    title = 'All Balls Equal - All Overlapping, Negative Vertex'


# Calculate the vertex
my_vert = calc_vert(locs=ar([ar(_) for _ in locs]), rads=ar(rads))

# Make the
fig = plt.figure()
ax = fig.add_subplot(projection='3d')


# Plot the atoms
plot_atoms(locs, rads, fig=fig, ax=ax, res=10, alpha=0.3)
# Plot the vertices
plot_verts([my_vert[0]], [abs(my_vert[1])], fig=fig, ax=ax, spheres=True, res=10)

# Set the axes lines
ax.plot([-5, -4], [-5, -5], [0, 0])
ax.plot([-5, -5], [-4, -5], [0, 0])
ax.plot([-5, -5], [-5, -5], [0, 1])

# Set the axes labels
ax.text(x=-4, y=-5, z=-0.25, s='x')
ax.text(x=-5, y=-4, z=-0.25, s='y')
ax.text(x=-5.25, y=-5.25, z=1, s='z')

# Set the scales for the figure
ax.set_xlim(-5, 5)
ax.set_ylim(-5, 5)
ax.set_zlim(-5, 5)

# Set the title
ax.set_title(title, font=dict(size=20, family='serif'))

# Show the plot
plt.show()
