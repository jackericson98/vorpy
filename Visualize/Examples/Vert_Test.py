from System.system import *

atoms = [[[27.314, 30.779, 9.417], 1.52], [[27.174, 30.469, 8.517], 1.2],
         [[26.454, 30.669, 9.837], 1.2], [[24.724, 30.199, 10.017], 1.52]]

sys = System(atoms)
sys.net.build()

print(sys.net.verts[0].rad)
