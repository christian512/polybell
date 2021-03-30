""" Calculates the vertices on the facet """

from linearbell.utils import get_deterministic_behaviors
from linearbell.panda_helper import read_inequalities
import numpy as np

relabels = np.loadtxt('../data/relabels_dets/4422.gz').astype(int)

print('TODO: Define Facet1 and Facet2 for 4422 case')
facet1 = [1, 2, 9, 11, 18, 20, 27, 31, 33, 34, 41, 45, 50, 54, 61, 63, 84, 88, 95, 96, 118, 120, 127, 128]
facet2 = [1, 2, 9, 11, 17, 25, 29, 34, 61, 66, 68, 75, 76, 79, 80, 93, 95, 98, 102, 110, 112, 118, 125, 126]

# subtract one because of pythonic index start (vs. GAP)
facet1 = [f-1 for f in facet1]
facet2 = [f-1 for f in facet2]
print('Starting equivalence test')
for r in relabels:
    tmp = np.copy(facet1)
    for i in range(len(facet1)):
        tmp[i] = r[tmp[i]]
    if np.all(tmp == facet2):
        print('equiv')
print('Not equiv')




