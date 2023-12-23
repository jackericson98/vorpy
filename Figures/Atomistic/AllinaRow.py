import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D

# Points and labels
points = [1, 2, 4, 5, 6, 6.5, 7, 7.1, 8.5, 9.2, 10]
labels = ['one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine', 'ten', 'eleven']

# Plotting the points
plt.figure(figsize=(8, 4))  # Adjust the figure size as needed
plt.plot(points, [0] * len(points), 'ro')  # 'ro' means red color, circle markers

pad = 0.1
# Adding labels for each point with diagonal rotation
for point, label in zip(points, labels):
    text = plt.text(point + pad, 0.01, f'{label}', ha='right', va='bottom', rotation=90)

    # # Adjust the starting position of the labels
    # trans = Affine2D().translate(-12, 0) + plt.gca().transData
    # text.set_transform(trans)

# Drawing lines from points to labels
for point, label in zip(points, labels):
    plt.plot([point, point], [0, 0.01], 'k-', linewidth=0.7)

# Setting up plot limits and removing y-axis
plt.xlim(min(points) - 1, max(points) + 1)
plt.ylim(-0.1, 0.1)
plt.yticks([])  # Remove y-axis ticks and labels

# Display the plot
plt.show()
