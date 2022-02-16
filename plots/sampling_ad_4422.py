""" Plots the calculation time for different recursion levels of the sampling method """

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 12})
# Define quantities
recursion_levels = [20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5]
times = [277, 289, 318, 351, 10546, 11173, 11470, 11054, 11381, 11355, 13857, 14930, 19594, 27614, 58994, 125447]
num_classes_found = [2, 2, 2, 2, 162, 163, 161, 164, 162, 163, 164, 167, 171, 172, 174, 175]

assert len(recursion_levels) == len(times)
assert len(recursion_levels) == len(num_classes_found)

fig, ax = plt.subplots()
ax.scatter(recursion_levels, num_classes_found, color='green', marker='x')
ax.set_ylabel('Num. of found facet-classes', color='green')
ax.set_xlabel('Recursion depth')
ax2 = ax.twinx()

ax2.scatter(recursion_levels, np.array(times) / 1000, color='blue')
ax2.set_ylabel('Run-time [s]', color='blue')
ax2.set_ylim(-2, 135)
plt.savefig('4422_sampling_ad.png',dpi=300)
