import plotly.graph_objects as go
import plotly.express as px

import numpy as np

x = y = z = np.arange(1, 4)

fig = go.Figure()
fig.add_trace(go.Scatter3d(
        mode='markers',
        x=x,
        y=y,
        z=z,
        marker=dict(color=px.colors.qualitative.D3,size=[50, 5, 25],sizemode='diameter'
        )
    )
)