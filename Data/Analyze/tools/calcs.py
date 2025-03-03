import numpy as np


def bin_data(points, xnum_bins, ynum_bins):
    # Get the maximum and minimum x and y
    maxx, minx = max([_[0] for _ in points]), min([_[0] for _ in points])
    maxy, miny = max([_[1] for _ in points]), min([_[1] for _ in points])
    # Get the x_step and y_step
    x_step, y_step = (maxx - minx) / xnum_bins, (maxy - miny) / ynum_bins
    # Create the bins
    bins = [[[] for _ in range(ynum_bins + 1)] for _ in range(xnum_bins + 1)]
    # Go through each point and assign them a bin
    for point in points:
        bins[int((point[0] - minx) / x_step)][int((point[1] - miny) / y_step)].append(point)
    # Create the count list and averages list
    counts, avgs = [], []
    # Average out the points
    for i in range(len(bins)):
        for j in range(len(bins[i])):
            counts.append(len(bins[i][j]))
            avgs.append(np.mean(bins[i][j]))
    # Return the data
    return bins, avgs, counts


def sort_lists(*lists):
    return [list(_) for _ in zip(*sorted(zip(*lists), key=lambda x: x[0]))]

