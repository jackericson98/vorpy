import numpy as np

import matplotlib.pyplot as plt

from scipy import stats
# x=np.linspace(0.1,0.1,200)
# bubs = []
# while len(bubs) < 1000:
#     data = stats.lognorm(1, loc=0, scale=1).rvs(size=1)[0]
#     if data > 0:
#         bubs.append(data)
#
# plt.hist(bubs, bins=100)
# plt.show()
sds = [0.01, 0.05, 0.09]
numbers = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000]
mu = 0.1
num_its = 100
for sd in sds:
    sd1s, sd2s, sd3s = [], [], []
    sd1gsl, sd2gsl, sd3gsl = [], [], []
    print(sd)
    for i in numbers:
        sd1, sd2, sd3 = [], [], []
        sd1gs, sd2gs, sd3gs = 0, 0, 0
        for j in range(num_its):
            bubs = []
            while len(bubs) < i:
                data = stats.gamma(loc=mu, scale=sd, a=1).rvs(size=1)[0]
                if data > 0:
                    bubs.append(data)
            sd1_count, sd2_count, sd3_count = 0, 0, 0
            for bub in bubs:
                if abs(bub - mu) < sd:
                    sd1_count += 1
                if abs(bub - mu) < 2 * sd:
                    sd2_count += 1
                if abs(bub - mu) < 3 * sd:
                    sd3_count += 1
            sd1.append(sd1_count)
            sd2.append(sd2_count)
            sd3.append(sd3_count)
            sd1g, sd2g, sd3g = 0, 0, 0
            if sd1_count / i > 0.68:
                sd1g = 1
            if sd2_count / i > 0.95:
                sd2g = 1
            if sd3_count / i > 0.99:
                sd3g = 1
            sd1gs += sd1g
            sd2gs += sd2g
            sd3gs += sd3g
        # print("{} Done".format(i))
        print(i, (sum(sd1)/num_its) / i, (sum(sd2)/num_its) / i, (sum(sd3)/num_its) / i, sd1gs, sd2gs, sd3gs)
