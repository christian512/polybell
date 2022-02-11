# Define input data
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 14})

# This is data with standard deviations, but I just plotted the data without, since the STD is relatively small.
scenarios = [[2, 2, 2, 2], [3, 2, 2, 2], [2, 2, 3, 2], [2, 2, 3, 3], [3, 3, 2, 2], [5, 2, 2, 2], [4, 3, 2, 2],
             [3, 2, 2, 3], [3, 2, 3, 2]]
pol_sizes = [320, 924, 1176, 28728, 11220, 4420, 239552, 25308, 7200]
symm_adm_times = [25.631, 38.272, 56.844, 1027.898, 285.158, 279.059, 5713.162, 648.102, 319.304]
symm_adm_std = [1.337, 4.155, 2.171, 11.399, 9.396, 1.939, 80.604, 5.265, 12.737]
dd_times = [5.27, 5.63, 6.121, 29.974, 21.263, 19.293, 13574.285, 36.359, 16.707]
dd_std = [0.688, 1.06, 0.875, 1.506, 1.25, 1.707, 122.824, 1.516, 2.028]

# Timings for DD method
scenarios_dd = [[2, 2, 2, 2], [3, 2, 2, 2], [2, 2, 3, 2], [2, 2, 3, 3], [3, 3, 2, 2], [5, 2, 2, 2], [4, 3, 2, 2],
                [3, 2, 2, 3], [3, 2, 3, 2], [5, 3, 2, 2]]
pol_sizes_dd = [320, 924, 1176, 28728, 11220, 4420, 239552, 25308, 7200, 1646708]
dd_times = [5.27, 5.63, 6.121, 29.974, 21.263, 19.293, 13574.285, 36.359, 16.707, 55709719]

# The symmetric ADM can run for a few higher cases in reasonable time too so I display them here as well
scenarios_adm = [[2, 2, 2, 2], [3, 2, 2, 2], [2, 2, 3, 2], [2, 2, 3, 3], [3, 3, 2, 2], [5, 2, 2, 2], [4, 3, 2, 2],
                 [3, 2, 2, 3], [3, 2, 3, 2], [5, 3, 2, 2], [3, 3, 2, 3]]
pol_sizes_adm = [320, 924, 1176, 28728, 11220, 4420, 239552, 25308, 7200, 1646708, 6824898]
symm_adm_times = [25.631, 38.272, 56.844, 1027.898, 285.158, 279.059, 5713.162, 648.102, 319.304, 1181236.984,
                  6053190.3014183]

# Check correct sizes of arrays
assert len(scenarios_dd) == len(pol_sizes_dd)
assert len(scenarios_dd) == len(dd_times)
assert len(scenarios_adm) == len(pol_sizes_adm)
assert len(scenarios_adm) == len(symm_adm_times)

# Errorbar plots
# plt.errorbar(pol_sizes, symm_adm_times, yerr=symm_adm_std, fmt='.', capsize=5,label='Symm. AD method')
# plt.errorbar(pol_sizes, dd_times, yerr=dd_std,fmt='.', capsize=5, label='DD method')

# Scatter plots
plt.scatter(pol_sizes_dd, dd_times, label='DD method')
plt.scatter(pol_sizes_adm, symm_adm_times, label='Symm. AD method')

plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Run-time [ms]')
plt.xlabel('Polytope size')
plt.savefig('output/symadm_vs_dd.png', dpi=200)
plt.close()
