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
    my_x = np.linspace(0, 10, 10000)[1:]
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

alphas = [100, 64, 44.4444, 32.65306, 25, 19.75309, 16, 13.22314, 11.11111, 9.46746, 8.16327, 7.11111, 6.25, 5.53633, 4.93827, 4.43213, 4]
betas = [0.00010, 0.00024, 0.00051, 0.00094, 0.00160, 0.00256, 0.00391, 0.00572, 0.00810, 0.01116, 0.01501, 0.01978, 0.02560, 0.03263, 0.04101, 0.05091, 0.06250]



# plot_function(lognormal, 'Lognormal Distributions by Sigma Value', 1, [round((i+4)*0.025, 3) for i in range(17)])
plot_function(gamma, 'Gamma Distributions by Beta Value - \u03B1 = 4', 4, [round((i+4)*0.5, 3) for i in range(17)])
# plot_function(physical_DeVries, "De Vries")
# plot_function(physical_Ranadive_Lemlich, "Ranadive Lemlich")
# plot_function(physical_GalOr_Hoelscher, "Gal-Or Hoelscher")
