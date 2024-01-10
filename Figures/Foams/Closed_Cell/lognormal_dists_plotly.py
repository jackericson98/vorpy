import plotly.express as px
import numpy as np
import scipy.stats as stats
from pandas import DataFrame as df


def lognormal(r, mu, sd):
    sd = np.log(sd)
    return -(1/(r*sd*np.sqrt(2*np.pi))) * np.exp(-((np.log(r) - mu)**2/(2*sd**2)))


def gamma(r, a, B):
    return 0.5*stats.gamma.pdf(r, a, scale=1/B)


def physical_DeVries(r, p1=None, p2=None):
    return 2.082*r/(1+0.387*r**2)**4


def physical_Ranadive_Lemlich(r, p1=None, p2=None):
    return (32/np.pi**2)*r**2*np.exp(-(4/np.pi)*r**2)


def physical_GalOr_Hoelscher(r, p1=None, p2=None):
    return (16/np.pi)*r**2*np.exp(-(16/np.pi)**0.5*r**2)


def plot_function(function, function_name="", p1=None, p2=None):
    if p1 is not None and type(p1) == list:
        my_x = np.linspace(0, 10, 10000)
        my_y = []
        for val in p1:
            my_y.append([function(_, val, p2) for _ in my_x])
        data = df(data=np.array(my_y).T)
        fig = px.line(data, x='x', y='y', title=function_name)
        fig.show()
    elif p2 is not None and type(p2) == list:
        my_y = []
        my_x = np.linspace(0, 5, 1000)[1:]
        for val in p2:
            my_y.append([function(_, p1, val) for _ in my_x])

        data = df(index=my_x, data=np.array(my_y).T, columns=p2)

        fig = px.line(data, title=function_name)
        fig.update_layout(dict1=dict(xaxis=dict(title='Bubble Radius', tickfont=dict(size=25), titlefont=dict(size=30)),
                                     yaxis=dict(title='Probability Distribution', tickfont=dict(size=25), titlefont=dict(size=30))),
                          title=dict(font=dict(size=50)),
                          legend=dict(title='CV', font=dict(size=25)))
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


plot_function(lognormal, 'Lognormal Distributions by Sigma Value', 1, [round((i+4)*0.025, 3) for i in range(17)])
# plot_function(gamma, 'Gamma Distributions by Beta Value - Alpha = 4', 4, [round((i+6)*0.25, 3) for i in range(21)])
# plot_function(physical_DeVries, "De Vries")
# plot_function(physical_Ranadive_Lemlich, "Ranadive Lemlich")
# plot_function(physical_GalOr_Hoelscher, "Gal-Or Hoelscher")
