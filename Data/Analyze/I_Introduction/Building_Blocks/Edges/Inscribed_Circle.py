from System.sys_funcs.calcs.circle import calc_circ
from Visualize.mpl_visualize import plot_balls, plot_circles
import matplotlib.pyplot as plt


locs = [[0, 5, 0], [5, 0, 0], [0, 0, 0]]
rads = [1, 2, 1.5]
my_circ = calc_circ(*locs, *rads)

olap_loc = [4, 4, 0]
olap_rad = 1

cor_circ = calc_circ(locs[0], locs[2], olap_loc, rads[0], rads[2], olap_rad)


fig = plt.figure()
ax = fig.add_subplot(projection='3d')

plot_circles([cor_circ[0]], [cor_circ[1]], fig=fig, ax=ax)
plot_balls(locs + [olap_loc], rads + [olap_rad], colors=['r', 'b', 'g', 'y'], fig=fig, ax=ax)

ax.set_xlim([-10, 10])
ax.set_ylim([-10, 10])
ax.set_zlim([-10, 10])

plt.show()
