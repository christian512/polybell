""" Plots calculation time for different recursion levels """
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 12})

# +++ PLOTS FOR 5322 +++
# Define times for 5322
recursion_levels_5322 = [0, 1, 2]
local_5322 = [1181236, 495837, 2297341]
stabilizer_5322 = [1181236, 427689, 1310019]
global_5322 = [1181236, 461629, 4162278]
dd_time_5322 = 55709719

# postions of bars
barWidth = 0.25
r1 = np.arange(len(recursion_levels_5322))
r2 = [x + barWidth for x in r1]
r3 = [x + barWidth for x in r2]

# plot the bars
plt.bar(r1, stabilizer_5322, width=barWidth, edgecolor='white', label='Stabilizer')
plt.bar(r2, global_5322, width=barWidth, edgecolor='white', label='Global')
plt.bar(r3, local_5322, width=barWidth, edgecolor='white', label='Local')
plt.axhline(dd_time_5322, color='red',label='DD runtime')

plt.xlabel('Recursion depth')
plt.xticks([r + barWidth for r in range(len(recursion_levels_5322))], [str(i) for i in recursion_levels_5322])
plt.ylabel('Run-time [ms]')
plt.yscale('log')
plt.legend()
plt.savefig('5322_recursion_depths.png', dpi=300)
plt.close()
# Define times for 4322
recursion_levels_4322 = [0, 1, 2]
local_4322 = [5713, 20687, 249762]
stabilizer_4322 = [5713, 16438, 120660]
global_4322 = [5713, 14159, 243867]
dd_time_4322 = 13574

# postions of bars
barWidth = 0.25
r1 = np.arange(len(recursion_levels_4322))
r2 = [x + barWidth for x in r1]
r3 = [x + barWidth for x in r2]

# plot the bars
plt.bar(r1, stabilizer_4322, width=barWidth, edgecolor='white', label='Stabilizer')
plt.bar(r2, global_4322, width=barWidth, edgecolor='white', label='Global')
plt.bar(r3, local_4322, width=barWidth, edgecolor='white', label='Local')
plt.axhline(dd_time_4322, color='red',label='DD runtime')


plt.xlabel('Recursion depth')
plt.xticks([r + barWidth for r in range(len(recursion_levels_4322))], [str(i) for i in recursion_levels_4322])
plt.ylabel('Run-time [ms]')
plt.yscale('log')
plt.legend()
plt.savefig('4322_recursion_depths.png', dpi=300)
