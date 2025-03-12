from Visualize.mpl_visualize import plot_surfs, plot_balls
import matplotlib.pyplot as plt


ball1 = [-2, 0, 0], 1
ball2 = [3, 0, 0], 2

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.set_xlim([-5, 5])
ax.set_ylim([-5, 5])
ax.set_zlim([-5, 5])

plot_balls([ball1[0], ball2[0]], [ball1[1], ball2[1]], colors=['k', 'k'], alpha=0.5, fig=fig, ax=ax)
plot_surfs([[[0, 5, 5], [0, 5, -5], [0, -5, 5]]], [[[0, 1, 2]]], simps=True, fig=fig, ax=ax, colors=['k'], alpha=0.2)
plot_surfs([[[0, -5, -5], [0, 5, -5], [0, -5, 5]]], [[[0, 1, 2]]], simps=True, fig=fig, ax=ax, colors=['k'], alpha=0.2)
ax.plot([ball1[0][0], ball2[0][0]], [ball1[0][1], ball2[0][1]], [ball1[0][2], ball2[0][2]])
ax.scatter([0], [0], [0])
plt.show()
