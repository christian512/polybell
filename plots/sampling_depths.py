""" Plots the r_max and calculation time for the sampling AD method """

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 12})

polys = [[4, 3, 2, 2], [5, 3, 2, 2], [3, 3, 2, 3], [6, 3, 2, 2], [4, 4, 2, 2]]
polys_str = ["(4,3,2,2)", "(5,3,2,2)", "(3,3,2,3)", "(6,3,2,2)", "(4,4,2,2)"]
dimensions = [19, 23, 27, 27, 24]
classes = [6, 7, 25, 7, 175]
reduced_sizes = np.array(dimensions) * np.array(classes)

r_max = [4, 6, 13, 7, 5]
times = [691, 1680, 2005, 34924, 125447]

r1 = np.arange(len(polys))

fig, ax = plt.subplots()
ax.scatter(r1, r_max, color='green', marker='x')
ax.set_ylabel('Maximum recursion depth', color='green')
ax.set_xlabel('Scenarios')
ax.set_ylim(0, 14)
ax.set_xticks(range(len(polys)), polys_str)

ax2 = ax.twinx()
ax2.scatter(r1, np.array(times) / 1000, color='blue')
ax2.set_ylabel('Run-time [s]', color='blue')
ax2.set_yscale('log')
plt.savefig('sampling_depths.png',dpi=300)
