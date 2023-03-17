import sys

from vispy import scene
from vispy.visuals.transforms import STTransform
from System.system import System

canvas = scene.SceneCanvas(keys='interactive', bgcolor='white',
                           size=(800, 600), show=True)

view = canvas.central_widget.add_view()
view.camera = 'fly'

my_sys = System(file='../Data/test_data/Na5.pdb')


my_spheres = []
for atom in my_sys.atoms:

    sphere = scene.visuals.Sphere(radius=atom.rad, method='ico', parent=view.scene, edge_color='black')
    sphere.transform = STTransform(translate=atom.loc)


view.camera.set_range(x=[20, 60])

if __name__ == '__main__' and sys.flags.interactive == 0:
    canvas.app.run()
