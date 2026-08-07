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
from matplotlib.animation import FuncAnimation, PillowWriter
import random
import tkinter as tk
from tkinter import filedialog


def plot_proximity(radius1=1, radius2=3, distance_between_circles=5, num_points=2000, markersize=1,
                   distance_type="surface", bounding_box=[[-10, -5], [10, 5]], alpha=0.5, show_center=False,
                   animation=False, show_animation=True, save_animation=False,
                   animation_file="proximity.gif", points_per_frame=10, fps=30):
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
    # Plot the center point of the circles
    if show_center:
        ax.scatter([-0.5 * distance_between_circles], [0], c='r', s=20*markersize, marker='x', alpha=alpha)
        ax.scatter([0.5 * distance_between_circles], [0], c='b', s=20*markersize, marker='x', alpha=alpha)

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
    if animation:
        # Determine point colors based on proximity
        red_points = []
        blue_points = []

        for x, y in points:
            dist_to_circle1 = calc_distance(
                x, y, -0.5 * distance_between_circles, 0, radius1
            )
            dist_to_circle2 = calc_distance(
                x, y, 0.5 * distance_between_circles, 0, radius2
            )

            if dist_to_circle1 < dist_to_circle2:
                red_points.append((x, y))
            else:
                blue_points.append((x, y))

        # Empty scatter plots that will be filled during animation
        red_scatter = ax.scatter([], [], color='r', marker='.',
                                 s=markersize, alpha=alpha)

        blue_scatter = ax.scatter([], [], color='b', marker='.',
                                  s=markersize, alpha=alpha)

        # Keep the original random order
        colored_points = []

        for x, y in points:
            dist_to_circle1 = calc_distance(
                x, y, -0.5 * distance_between_circles, 0, radius1
            )
            dist_to_circle2 = calc_distance(
                x, y, 0.5 * distance_between_circles, 0, radius2
            )

            color = 'r' if dist_to_circle1 < dist_to_circle2 else 'b'
            colored_points.append((x, y, color))
    else:
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
    if animation:
        # Lock the plot limits so they do not change during animation
        ax.set_xlim(bounding_box[0][0], bounding_box[1][0])
        ax.set_ylim(bounding_box[0][1], bounding_box[1][1])

        def update(frame):
            # Build this frame from scratch so previewing and saving do not
            # interfere with one another.
            end = min((frame + 1) * points_per_frame, len(colored_points))
            visible_points = colored_points[:end]

            red_xy = [(x, y) for x, y, color in visible_points if color == 'r']
            blue_xy = [(x, y) for x, y, color in visible_points if color == 'b']

            red_scatter.set_offsets(
                np.asarray(red_xy).reshape(-1, 2) if red_xy else np.empty((0, 2))
            )
            blue_scatter.set_offsets(
                np.asarray(blue_xy).reshape(-1, 2) if blue_xy else np.empty((0, 2))
            )

            return red_scatter, blue_scatter

        num_frames = int(np.ceil(len(colored_points) / points_per_frame))

        anim = FuncAnimation(
            fig,
            update,
            frames=num_frames,
            interval=1000 / fps,
            blit=True,
            repeat=False
        )

        # Either option can be used independently, or both can be True.
        if save_animation:
            anim.save(
                animation_file,
                writer=PillowWriter(fps=fps)
            )

        if show_animation:
            plt.show()
        else:
            plt.close(fig)
    else:
        my_input = input('Save figure? (y/n)  >>>>    ')
        if my_input.lower() in {'y', 'yes'}:
            root = tk.Tk()
            root.withdraw()
            root.wm_attributes('-topmost', 1)
            location = filedialog.asksaveasfilename(defaultextension='.png',)
            plt.savefig(location, format='png', dpi=1200)
        elif my_input.lower() in {'n', 'no'}:
            plt.show()


# Example usage:
plot_proximity(distance_type="prm", num_points=10000, alpha=0.6, markersize=0.9, radius1=0, radius2=0,
               distance_between_circles=4, bounding_box=[[-5, -4], [5, 4]], show_center=True,
               animation=True, show_animation=True, save_animation=True,
               animation_file="proximity.gif", points_per_frame=10, fps=30)