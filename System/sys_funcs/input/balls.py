import numpy as np
from System.sys_objs.atom import make_atom
import pandas as pd


def load_balls(sys):
    #
    balls = []
    for i, ball in enumerate(sys.atoms):
        loc = np.array([float(_) for _ in ball[0]])
        rad = float(ball[1])
        balls.append(make_atom(sys, loc, rad, i))

    sys.balls = pd.DataFrame(balls)
    sys.residues = []
    sys.chains = []
    sys.name = 'balls'
