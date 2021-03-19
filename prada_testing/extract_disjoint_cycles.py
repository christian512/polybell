from linearbell.utils import get_deterministic_behaviors
from linearbell.panda_helper import write_known_vertices, read_inequalities
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import time


relabels_dets = np.loadtxt('../data/relabels_dets/{}{}{}{}.gz'.format(2, 2, 2, 2)).astype(int)

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
        for val in cycle:
            out += str(val) + " "
        out = out[:-1] + ","
    out = out[:-1] + "\n"
print(out)