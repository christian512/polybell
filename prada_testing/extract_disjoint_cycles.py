from linearbell.utils import get_deterministic_behaviors
from linearbell.panda_helper import write_known_vertices, read_inequalities
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import time

relabels_dets = np.loadtxt('../data/relabels_dets/{}{}{}{}.gz'.format(4, 3, 2, 2)).astype(int)
print('loaded relabels, extracting disjoint cycles')

# storage for all cycles
all_cycles = []

for relabel in relabels_dets:
    # convert to dict
    d = {}
    for i in range(relabel.shape[0]):
        d[i] = relabel[i]
    # generate cycles
    cycles = []
    while len(d) > 0:
        cycle = []
        key = list(d)[0]
        val = d.pop(key)
        if val != key:
            cycle.append(key)
            cycle.append(val)
        while val in list(d):
            key = val
            val = d.pop(key)
            if val not in cycle:
                cycle.append(val)
        if cycle:
            cycles.append(cycle)
    if cycles:
        all_cycles.append(cycles)

# Store the cycles
out = ""
for cycles in all_cycles:
    for cycle in cycles:
        out = out + "("
        for val in cycle:
            out += str(val + 1) + ","
        out = out[:-1] + ")"
    out = out + ",\n"
out = out[:-2]
f = open('disjoint_cycles.ext', 'w+')
f.write(out)
f.close()
print('Number of exported cycles: ', relabels_dets.shape[0])

# store vertices for polytope
dets = get_deterministic_behaviors(range(4), range(3), range(2))
out = ""
for d in dets:
    out += str(list(d.astype(int))) + ",\n"
out = out[:-2]
print(out)
f = open('vertices.ext', 'w+')
f.write(out)
f.close()
print('Number of exported vertices: ', dets.shape[0])
