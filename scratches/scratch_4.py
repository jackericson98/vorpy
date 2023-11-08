import plotly.graph_objects as go
import numpy as np

x = np.outer(np.linspace(-2, 2, 30), np.ones(30))

# transpose
y = x.copy().T
z = np.cos(x ** 2 + y ** 2)
print(x, y, z)
fig = go.Figure(data=[go.Surface(x=x, y=y, z=z)])

fig.update_traces(contours_z=dict(
    show=True, usecolormap=True,
    highlightcolor="limegreen",
    project_z=True))

fig.show()