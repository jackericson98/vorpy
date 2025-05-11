import matplotlib.pyplot as plt
from Visualize.mpl_visualize import plot_balls

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
plot_balls([[0.2, 0, 0.0], [0.0, 0.0, 0.0]], [1.0, 3.0], alpha=0.1, fig=fig, ax=ax, res=10)
plot_balls([[0.1, 0, 0]], [2.0], alpha=0.2, colors=['k'], fig=fig, ax=ax, Show=True, res=10)

