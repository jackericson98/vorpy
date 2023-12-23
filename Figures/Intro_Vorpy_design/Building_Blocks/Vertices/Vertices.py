from Visualize.mpl_visualize import *
from System.sys_objs.atom import make_atom
from System.Network.verts.find_verts import calc_vert
"""Example for vertices"""

# With vertices, we have __ different cases
cases = []
# Case 0: No overlap
cases.append([make_atom(location=[1., 0., 0.5], radius=0.5), make_atom(location=[-1.0, 0., 0.5], radius=0.5), make_atom(location=[0., 1., -0.5], radius=0.5), make_atom(location=[0., -1., -0.5], radius=0.5)])
# Case 1: One overlapping ball
cases.append([make_atom(location=[.1, 0., 0.5], radius=0.5), make_atom(location=[-.1, 0., 0.5], radius=0.5), make_atom(location=[0., 1., -0.5], radius=0.5), make_atom(location=[0., -1., -0.5], radius=0.5)])
# Case 2: Two overlapping balls
cases.append([make_atom(location=[0.5, 0., 0.], radius=0.5), make_atom(location=[-0.5, 0., 0.], radius=0.5), make_atom(location=[-0.25, 0.25, 1.], radius=0.5), make_atom(location=[-0.25, -0.25, 1.], radius=0.75)])
# Case 3: Three overlapping balls
cases.append([make_atom(location=[0., 0., 1.], radius=0.25), make_atom(location=[0., 0., -1.], radius=0.25), make_atom(location=[1., 0.5, 0.], radius=0.25), make_atom(location=[-0.5, 0.5, 0.], radius=0.25)])
# Case 4: Four overlapping balls.
cases.append([make_atom(location=[0., 1., 1.], radius=1.3), make_atom(location=[0., 1., 0.], radius=1.3), make_atom(location=[-1., 0., 0.], radius=1.1), make_atom(location=[0., -1., 0.], radius=1.1)])


verts = []
for i in range(len(cases)):
    myVert = calc_vert(locs=np.array([np.array(_['loc']) for _ in cases[i]]), rads=np.array([_['rad'] for _ in cases[i]]))
    print(myVert)
    verts.append(myVert)


fig = plt.figure(figsize=(20, 40))
titles = ["Case 0: No Overlap", "Case 1: One overlapping set", "Case 2: Two overlapping balls",
          "Case 3: Three overlapping balls", "Case 4: Four overlapping balls"]
for i in range(len(verts)):
    axn = fig.add_subplot(int("23" + str(i + 1)), projection="3d", xlim=10)
    axn.set_title(titles[i])
    plot_atoms(cases[i], fig=fig, ax=axn, colors=['y', 'y', 'y', 'y'], alpha=.1, dfo=5)
    plot_verts([verts[i]], fig=fig, ax=axn, dfo=15, spheres=True)

plt.show()
