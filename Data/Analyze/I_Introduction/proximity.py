# import numpy as np
# import matplotlib.pyplot as plt
#
# # Define circle parameters
# radius1 = 1
# radius2 = 3
# distance_between_circles = 5
#
# num_points = 2000
# markersize = 1
#
# # Sort the radii
# if radius1 > radius2:
#     radius1, radius2 = radius2, radius1
#
# # Get the areas
# s_ar = np.pi * radius1 ** 2  # Small Area
# l_ar = np.pi * radius2 ** 2  # Large Area
# f_ar = np.pi * (2 * radius2) ** 2  # Double the radius of the large area or the total area for both circles reach
#
# # Get the relative ratios
#
# half_points = num_points // 2
# small_in_num, small_out_num = int(half_points * s_ar / f_ar), int(half_points * (f_ar - s_ar) / f_ar)
# large_in_num, large_out_num = int(half_points * l_ar / f_ar), int(half_points * (f_ar - l_ar) / f_ar)
#
#
# # Generate points around the circles
# theta = np.linspace(0, 2*np.pi, num_points)
# x_circle1 = radius1 * np.cos(theta)
# y_circle1 = radius1 * np.sin(theta)
# x_circle2 = radius2 * np.cos(theta) + distance_between_circles
# y_circle2 = radius2 * np.sin(theta)
#
# # Generate random points both inside and outside the circles
# np.random.seed()  # for reproducibility
# angle_rand_inside = np.random.uniform(0, 2 * np.pi, small_in_num)
# radius_rand_inside = np.random.uniform(0, radius1, small_in_num)
# x1_rand_inside = np.concatenate([radius_rand_inside * np.cos(angle_rand_inside),
#                                 distance_between_circles + radius_rand_inside * np.cos(angle_rand_inside)])
# y1_rand_inside = np.concatenate([radius_rand_inside * np.sin(angle_rand_inside),
#                                 radius_rand_inside * np.sin(angle_rand_inside)])
#
# angle_rand_outside = np.random.uniform(0, 2 * np.pi, small_out_num)
# radius_rand_outside = np.random.uniform(radius1, 3 * radius2, small_out_num)
# x1_rand_outside = np.concatenate([radius_rand_outside * np.cos(angle_rand_outside),
#                                  distance_between_circles + radius_rand_outside * np.cos(angle_rand_outside)])
# y1_rand_outside = np.concatenate([radius_rand_outside * np.sin(angle_rand_outside),
#                                  radius_rand_outside * np.sin(angle_rand_outside)])
#
# # Generate random points both inside and outside the circles
# np.random.seed()  # for reproducibility
# angle_rand_inside = np.random.uniform(0, 2 * np.pi, large_in_num)
# radius_rand_inside = np.random.uniform(0, radius1, large_in_num)
# x2_rand_inside = np.concatenate([radius_rand_inside * np.cos(angle_rand_inside),
#                                 distance_between_circles + radius_rand_inside * np.cos(angle_rand_inside)])
# y2_rand_inside = np.concatenate([radius_rand_inside * np.sin(angle_rand_inside),
#                                 radius_rand_inside * np.sin(angle_rand_inside)])
#
# angle_rand_outside = np.random.uniform(0, 2 * np.pi, large_out_num)
# radius_rand_outside = np.random.uniform(radius1, 3 * radius2, large_out_num)
# x2_rand_outside = np.concatenate([radius_rand_outside * np.cos(angle_rand_outside),
#                                  distance_between_circles + radius_rand_outside * np.cos(angle_rand_outside)])
# y2_rand_outside = np.concatenate([radius_rand_outside * np.sin(angle_rand_outside),
#                                  radius_rand_outside * np.sin(angle_rand_outside)])
#
# # Create figure and axis
# fig, ax = plt.subplots()
#
# # Plot circles
# ax.plot(x_circle1, y_circle1, 'r')
# ax.plot(x_circle2, y_circle2, 'b')
#
# # Plot random points inside and outside the circles
# # ax.scatter(x_rand_inside, y_rand_inside, color='blue', label='Points Closer to Blue Circle')
# # ax.scatter(x_rand_outside, y_rand_outside, color='red', label='Points Closer to Red Circle')
#
# # Color the points based on the closest circle
# for x, y in zip(x1_rand_inside, y1_rand_inside):
#     distance_to_circle1 = np.sqrt((x - x_circle1)**2 + (y - y_circle1)**2)
#     distance_to_circle2 = np.sqrt((x - x_circle2)**2 + (y - y_circle2)**2)
#     if np.min(distance_to_circle1) < np.min(distance_to_circle2):
#         ax.scatter(x, y, color='r', marker='.', s=markersize)
#     else:
#         ax.scatter(x, y, color='b', marker='.', s=markersize)
#
# for x, y in zip(x2_rand_outside, y2_rand_outside):
#     distance_to_circle1 = np.sqrt((x - x_circle1)**2 + (y - y_circle1)**2)
#     distance_to_circle2 = np.sqrt((x - x_circle2)**2 + (y - y_circle2)**2)
#     if np.min(distance_to_circle1) < np.min(distance_to_circle2):
#         ax.scatter(x, y, color='r', marker='.', s=markersize)
#     else:
#         ax.scatter(x, y, color='b', marker='.', s=markersize)
#
# # Set aspect ratio to equal and add legend
# ax.set_aspect('equal', 'box')
# # ax.legend()
# plt.xticks([])
# plt.yticks([])
# plt.axis('off')
#
# # Show plot
# plt.show()
#
import numpy as np
import matplotlib.pyplot as plt
import random


def plot_proximity(radius1=1, radius2=3, distance_between_circles=5, num_points=2000, markersize=1,
                   distance_type="surface", bounding_box=[[-10, -5], [10, 5]], alpha=0.5):
    """
    Plots proximity of random points to two circles based on a chosen distance type.

    Parameters:
        radius1 (float): Radius of the smaller circle.
        radius2 (float): Radius of the larger circle.
        distance_between_circles (float): Distance between the centers of the two circles.
        num_points (int): Number of random points to generate.
        markersize (float): Size of the markers in the scatter plot.
        distance_type (str): Type of distance calculation ('surface', 'center', or 'power').

    """
    # Sort the radii
    if radius1 > radius2:
        radius1, radius2 = radius2, radius1

    # Generate points around the circles
    theta = np.linspace(0, 2 * np.pi, num_points)
    x_circle1 = radius1 * np.cos(theta) - 0.5 * distance_between_circles
    y_circle1 = radius1 * np.sin(theta)
    x_circle2 = radius2 * np.cos(theta) + 0.5 * distance_between_circles
    y_circle2 = radius2 * np.sin(theta)

    # Generate points at random
    points = []
    for _ in range(num_points):
        points.append((random.uniform(bounding_box[0][0], bounding_box[1][0]),
                       random.uniform(bounding_box[0][1], bounding_box[1][1])))

    # Create figure and axis
    fig, ax = plt.subplots()
    ax.plot(x_circle1, y_circle1, 'r')
    ax.plot(x_circle2, y_circle2, 'b')

    # Define the distance calculation based on selected type
    def calc_distance(x, y, x_circle, y_circle, radius):
        if distance_type == "aw":
            return np.sqrt((x - x_circle) ** 2 + (y - y_circle) ** 2) - radius
        elif distance_type == "prm":
            return np.sqrt((x - x_circle) ** 2 + (y - y_circle) ** 2)
        elif distance_type == "pow":
            return (x - x_circle) ** 2 + (y - y_circle) ** 2 - radius ** 2
        else:
            raise ValueError("Invalid distance type. Choose 'surface', 'center', or 'power'.")

    # Plot points with colors based on proximity
    for x, y in points:
        dist_to_circle1 = calc_distance(x, y, -0.5 * distance_between_circles, 0, radius1)
        dist_to_circle2 = calc_distance(x, y, 0.5 * distance_between_circles, 0, radius2)
        color = 'r' if dist_to_circle1 < dist_to_circle2 else 'b'
        ax.scatter(x, y, color=color, marker='.', s=markersize, alpha=alpha)

    # Finalize plot settings
    ax.set_aspect('equal', 'box')
    plt.xticks([])
    plt.yticks([])
    plt.axis('off')
    my_input = input('Save figure? (y/n)  >>>>    ')
    if my_input.lower() in {'y', 'yes'}:
        plt.savefig('C:/Users/Optiplex_7060/OneDrive - Georgia State University/GSU NSC/Manuscripts/'
                    'Ericson Voronoi DNA/Figures/P1/Figure1_Concepts_and_Scheme_Comparisons/PointProximity/PointProximityPow1.png', format='png', dpi=1200)


# Example usage:
plot_proximity(distance_type="pow", num_points=10000, alpha=0.6, markersize=0.6, radius1=1, radius2=np.sqrt(5),
               distance_between_circles=4, bounding_box=[[-5, -4], [5, 4]])
