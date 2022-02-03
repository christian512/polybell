"""
This file checks the degeneracy of a polytope by checking how many vertices are on its facet.
"""

from linearbell.panda_helper import read_inequalities
from linearbell.utils import get_deterministic_behaviors_two_party
import numpy as np
import matplotlib.pyplot as plt

# set scenario
ma,mb,na,nb = 4,4,2,2
dim = ma*(na-1)*mb*(nb-1) + ma * (na-1) + mb*(nb-1)
print('Calculations for the {}{}{}{} scenario'.format(ma,mb,na,nb))

# set file
filename = '../facet_classes/{}{}{}{}.ine'.format(ma,mb,na,nb)

# read inequalities from file
ineqs = read_inequalities(filename)
print('Loaded {} facets'.format(ineqs.shape[0]))

# TODO: Generate vertices and check if the are on each facet
vertices = get_deterministic_behaviors_two_party(range(ma), range(mb), range(na), range(nb),)
print('Generated {} vertices'.format(vertices.shape[0]))

# number of vertices per facet
num_vert_per_facet = []
for facet in ineqs:
    counter = 0
    # check the correct dimensions of the facet and vertices
    assert len(facet) == vertices.shape[1] + 1
    for v in vertices:
        if v @ facet[:-1] == -1.0 * facet[-1]:
            counter += 1
    num_vert_per_facet.append(counter)
num_vert_per_facet = np.array(num_vert_per_facet)

print('Number of vertices on each facet: ', num_vert_per_facet)
print('Dimension of polytope: ', dim )
print('Maximal number of vertices on each facet: ', np.max(num_vert_per_facet))
print('Minimal number of vertices on each facet: ', np.min(num_vert_per_facet))
print('Average number of vertices on each facet: ', np.average(num_vert_per_facet))

if dim <= np.min(num_vert_per_facet):
    print('Polytope is non-simplical')


