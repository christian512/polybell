# Define input data
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})

scenarios = [[2, 2, 2, 2], [3, 2, 2, 2], [2, 2, 3, 2], [3, 3, 2, 2], [5, 2, 2, 2], [3, 2, 2, 3], [3, 2, 3, 2],
             [2, 2, 3, 3]]
pol_sizes = [320, 924, 1176, 11220, 4420, 25308, 7200, 28728]

ad_lex_asc = [4.1, 11.3, 17.7, 374.3, 227.0, 795.7, 258.0, 880.3]
dd_lex_asc = [3.1, 4.2, 4.5, 21.7, 17.3, 36.0, 14.2, 30.3]
ad_random = [3.9,12.4,19.0,902.3,47993.4,2811.1,26823.6,4137.6]
dd_random = [2.7,4.0,4.6,510.4,98932.7,727.7,18309.7,520.8]
ad_lex_desc = [3.8,11.9,19.5,420.2,437.2,1081.5,656.4,963.9]
dd_lex_desc = [3.3,4.1,4.3,43.4,79.7,86.2,65.5,36.1]


assert len(scenarios) == len(pol_sizes)
assert len(scenarios) == len(ad_lex_asc)
assert len(scenarios) == len(dd_lex_asc)
assert len(scenarios) == len(ad_random)
assert len(scenarios) == len(dd_random)
assert len(scenarios) == len(ad_lex_desc)



plt.scatter(pol_sizes, dd_random, label='DD - random')
plt.scatter(pol_sizes, dd_lex_desc, label='DD - lex. desc.')
plt.scatter(pol_sizes, dd_lex_asc, label='DD - lex. asc.')

plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Run-time [ms]')
plt.xlabel('Polytope size')
plt.savefig('output/dd_input_comparison.png',dpi=200)
plt.close()

plt.scatter(pol_sizes, ad_random, label='AD - random')
plt.scatter(pol_sizes, ad_lex_desc, label='AD - lex. desc.')
plt.scatter(pol_sizes, ad_lex_asc, label='AD - lex. asc.')

plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Run-time [ms]')
plt.xlabel('Polytope size')
plt.savefig('output/ad_input_comparison.png',dpi=200)
plt.close()