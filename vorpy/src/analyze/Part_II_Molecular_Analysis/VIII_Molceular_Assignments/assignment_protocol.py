import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def new_function(random_number=0):
    if random_number % 2 == 1:
        plt.plot(np.random.randn(100))
    else:
        plt.plot(np.random.randn(1000))




def new_function_2(random_number=0):
    if random_number % 2 == 1:
        plt.plot(np.random.randn(100))
    else:
        plt.plot(np.random.randn(1000))
    print("Cumington Jones and company")
