"""
Generating a plot to show that Bell polytopes are non-simplical, i.e. facet enumeration is degenerate
"""

from linearbell.panda_helper import read_inequalities
from linearbell.utils import get_deterministic_behaviors_two_party
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


# set scenario
scenarios = [[3,2,2,2], [2,2,3,2], [2,2,3,3],[3,3,2,2],[5,2,2,2],[4,3,2,2],[5,3,2,2],[3,3,2,3]]


max_vert_per_facet_for_scenarios = []
min_vert_per_facet_for_scenarios = []
avg_vert_per_facet_for_scenarios = []
dimensions_for_scenarios = []

for scenario in scenarios:
    ma,mb,na,nb = scenario
    dim = ma*(na-1)*mb*(nb-1) + ma * (na-1) + mb*(nb-1)
    # set file
    filename = '../facet_classes/{}{}{}{}.ine'.format(ma,mb,na,nb)

    # read inequalities from file
    ineqs = read_inequalities(filename)

    # Generate vertices and check if the are on each facet
    vertices = get_deterministic_behaviors_two_party(range(ma), range(mb), range(na), range(nb),)

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

    # append to global documentation
    dimensions_for_scenarios.append(dim)
    max_vert_per_facet_for_scenarios.append(np.max(num_vert_per_facet))
    min_vert_per_facet_for_scenarios.append(np.min(num_vert_per_facet))
    avg_vert_per_facet_for_scenarios.append(np.average(num_vert_per_facet))

fig, ax = plt.subplots()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.scatter(dimensions_for_scenarios, max_vert_per_facet_for_scenarios, label='Maximal Incidence',zorder=1)
plt.scatter(dimensions_for_scenarios, min_vert_per_facet_for_scenarios, label='Minimal Incidence',zorder=2)
#plt.scatter(dimensions_for_scenarios, avg_vert_per_facet_for_scenarios, label='Avg incidence')
plt.scatter(dimensions_for_scenarios, dimensions_for_scenarios, label='Polytope dimension')
plt.xlabel('Polytope dimensions')
plt.legend()
plt.savefig('degeneracy.png',dpi=300)
