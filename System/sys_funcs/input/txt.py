import pandas as pd
import numpy as np
from System.sys_objs.atom import make_atom


def read_txt(sys, file=None):
    if file is None:
        file = sys.files['base_file']
    with open(file, 'r') as read_file:
        balls = []
        for i, line in enumerate(read_file.readlines()):
            line = line.split(" ")
            line = [_ for _ in line if _ != ""]
            loc = np.array([float(_) for _ in line[:3]])
            rad = float(line[3])
            balls.append(make_atom(sys, loc, rad, i))
    sys.balls = pd.DataFrame(balls)
    sys.residues = []
    sys.chains = []

