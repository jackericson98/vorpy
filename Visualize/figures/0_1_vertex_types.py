from Visualize.mpl_visualize import *
from System.sys_objs.atom import Atom
from System.Network.net_objs.vertex import Vertex

fig = plt.figure()
doublet = fig.add_subplot(131, projection='3d')

atoms1 = [Atom(location=[0, 0, 0], radius=1), Atom(location=[4, 0, .1], radius=2), Atom(location=[-2, 2, 0], radius=2), Atom(location=[-2, -2, 0], radius=2)]
vert1 = Vertex(atoms=atoms1)
vert1.calc_vert()
vert1.doublet = Vertex(atoms=vert1.atoms, location=vert1.loc2, radius=vert1.rad2)

doublet.set_title("Two Vertices")
plot_balls(atoms1, fig=fig, ax=doublet, colors=['blue'] * 4)
plot_verts([vert1, vert1.doublet], fig=fig, ax=doublet, spheres=True, colors=['g', 'g'])

singlet = fig.add_subplot(132, projection='3d')

atoms2 = [Atom(location=[1, 0, -1], radius=1), Atom(location=[4, 0, .1], radius=0.5), Atom(location=[-2, 2, 0], radius=2), Atom(location=[-2, -2, 0], radius=2)]
vert2 = Vertex(atoms=atoms2)
vert2.calc_vert()

singlet.set_title("Single Vertex")
plot_balls(atoms2, fig=fig, ax=singlet, colors=['blue'] * 4)
plot_verts([vert2], fig=fig, ax=singlet, spheres=True, colors=['g', 'g'])

nonelet = fig.add_subplot(133, projection='3d')

atoms3 = [Atom(location=[1, 0, -1], radius=1), Atom(location=[4, 0, -2], radius=0.5), Atom(location=[-2, 2, 0], radius=2), Atom(location=[-2, -2, 0], radius=2)]
vert3 = Vertex(atoms=atoms3)
vert3.calc_vert()

# Add the plot information
nonelet.set_title("No vertices")

plot_balls(atoms3, fig=fig, ax=nonelet, colors=['blue'] * 4, Show=True)
