import plotly.express as px
import numpy as np
import scipy.stats as stats
from pandas import DataFrame as df
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize


def lognormal(r, mu, sd):
    sd = np.log(sd)
    return -(1/(r*sd*np.sqrt(2*np.pi))) * np.exp(-((np.log(r) - mu)**2/(2*sd**2)))


def gamma(r, a, B):
    return 0.25*stats.gamma.pdf(r, a, scale=1/B)


def physical_DeVries(r, p1=None, p2=None):
    return 2.082*r/(1+0.387*r**2)**4


def physical_Ranadive_Lemlich(r, p1=None, p2=None):
    return (32/np.pi**2)*r**2*np.exp(-(4/np.pi)*r**2)


def physical_GalOr_Hoelscher(r, p1=None, p2=None):
    return (16/np.pi)*r**2*np.exp(-(16/np.pi)**0.5*r**2)


def plot_function1(function, function_name="", p1=None, p2=None):
    my_x = np.linspace(0, 5, 1000)[1:]
    my_y = []
    if p1 is not None and type(p1) == list and p2 is not None and type(p2) == list:
        for i in range(len(p1)):
            my_y.append([function(_, p1[i], p2[i]) for _ in my_x])
        data = df(index=my_x, data=np.array(my_y).T, columns=p2)
        fig = px.line(data, title=function_name)
        fig.update_layout(dict1=dict(xaxis=dict(title='Bubble Radius', tickfont=dict(size=25), titlefont=dict(size=30)),
                                     yaxis=dict(title='Probability Distribution', tickfont=dict(size=25), titlefont=dict(size=30))),
                          title=dict(font=dict(size=50)),
                          legend=dict(title='CV', font=dict(size=25)))
        fig.show()
    elif p1 is not None and type(p1) == list:
        for val in p1:
            my_y.append([function(_, val, p2) for _ in my_x])
        data = df(data=np.array(my_y).T)
        fig = px.line(data, x='x', y='y', title=function_name)
        fig.show()
    elif p2 is not None and type(p2) == list:
        for val in p2:
            my_y.append([function(_, p1, val) for _ in my_x])

        data = df(index=my_x, data=np.array(my_y).T, columns=p2)

        fig = px.line(data, title=function_name)
        fig.update_layout(dict1=dict(xaxis=dict(title='Bubble Radius', tickfont=dict(size=25), titlefont=dict(size=30)),
                                     yaxis=dict(title='Probability Distribution', tickfont=dict(size=25), titlefont=dict(size=30))),
                          title=dict(font=dict(size=50)),
                          legend=dict(title='\u03B2', font=dict(size=25)))
        fig.show()
    else:
        my_x = np.linspace(0, 10, 1000)
        my_y = [function(_, p1, p2) for _ in my_x]
        data = df(dict(x=my_x, y=my_y))
        fig = px.line(data, x='x', y='y', title=function_name)
        fig.update_layout(dict1=dict(xaxis=dict(title='Bubble Radius', tickfont=dict(size=18), titlefont=dict(size=25)),
                                     yaxis=dict(title='Probability Distribution', tickfont=dict(size=18),
                                                titlefont=dict(size=25))),
                          title=dict(font=dict(size=40)))
        fig.show()


def plot_function(function, p1=None, p2=None, title='', x_label='', y_label='', legend_title='', max_x=5):
    my_x = np.linspace(0, max_x, 10000)[1:]  # Avoiding zero if function undefined at zero
    fig, ax = plt.subplots()

    # Determine the parameter array and colormap
    cmap = plt.cm.rainbow  # You can choose any colormap (e.g., viridis, plasma, inferno, magma)

    if p1 is not None and isinstance(p1, list) and p2 is not None and isinstance(p2, list):
        # Assuming p1 and p2 are related and both need gradient representation
        # For simplicity, we'll only create a gradient based on p2 values here
        norm = Normalize(vmin=min(p2), vmax=max(p2))
        sm = ScalarMappable(norm=norm, cmap=cmap)

        for i in range(len(p1)):
            my_y = [function(x, p1[i], p2[i]) for x in my_x]
            ax.plot(my_x, my_y, color=sm.to_rgba(p2[i]))

    elif p1 is not None and isinstance(p1, list):
        # Only p1 has multiple values, color gradient by p1
        norm = Normalize(vmin=min(p1), vmax=max(p1))
        sm = ScalarMappable(norm=norm, cmap=cmap)

        for val in p1:
            my_y = [function(x, val, p2) for x in my_x]
            ax.plot(my_x, my_y, color=sm.to_rgba(val))

    elif p2 is not None and isinstance(p2, list):
        # Only p2 has multiple values, color gradient by p2
        norm = Normalize(vmin=min(p2), vmax=max(p2))
        sm = ScalarMappable(norm=norm, cmap=cmap)

        for val in p2:
            my_y = [function(x, p1, val) for x in my_x]
            ax.plot(my_x, my_y, color=sm.to_rgba(val))

    else:
        # Both p1 and p2 are single parameters, plot a single line
        my_y = [function(x, p1, p2) for x in my_x]
        ax.plot(my_x, my_y, color='blue')  # Default color

    ax.set_title(title)
    axis_font = {'fontname': 'Arial', 'size': '20'}
    ax.set_xlabel(x_label, **axis_font)
    ax.set_ylabel(y_label, **axis_font)
    ax.tick_params(axis='both', which='major', width=2, length=12, labelsize=15)
    ax.set_yticks([])

    size = 20

    # # Create colorbar
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label(legend_title, **axis_font)
    cbar.ax.tick_params(labelsize=15, size=10, width=2, length=12)
    plt.tight_layout()
    # ax.axis('off')
    plt.show()


alphas = [100, 64, 44.4444, 32.65306, 25, 19.75309, 16, 13.22314, 11.11111, 9.46746, 8.16327, 7.11111, 6.25, 5.53633, 4.93827, 4.43213, 4]
betas = [0.00010, 0.00024, 0.00051, 0.00094, 0.00160, 0.00256, 0.00391, 0.00572, 0.00810, 0.01116, 0.01501, 0.01978, 0.02560, 0.03263, 0.04101, 0.05091, 0.06250]



# plot_function(lognormal, 1, [round((i+4)*0.025, 3) for i in range(17)],
#               x_label='Bubble Radius', y_label='Probability', legend_title='Coefficient of Variation')
# plot_function(lognormal, 1, [0.45], 'Lognormal Distributions',
#               x_label='Bubble Radius', y_label='Probability', legend_title='Coefficient of Variation', max_x=20)
plot_function(gamma, 10, [round((i+4)*0.5, 3) for i in range(17)],
              x_label='Bubble Radius', y_label='Probability', legend_title='\u03b2 Value', max_x=5)
# plot_function(physical_DeVries, "De Vries")
# plot_function(physical_Ranadive_Lemlich, "Ranadive Lemlich")
# plot_function(physical_GalOr_Hoelscher, "Gal-Or Hoelscher")
