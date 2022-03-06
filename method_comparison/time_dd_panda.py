from polybell.utils import get_relabelling_generators, get_deterministic_behaviors_two_party, get_relabels_dets
from polybell.gap_helper import relabels_dets_to_disjoint_cycles
from polybell.panda_helper import write_known_vertices
import argparse
import numpy as np
import time
import os

np.random.seed(42)

test_configs = [[2, 2, 2, 2], [3, 2, 2, 2], [2, 2, 3, 2], [3, 3, 2, 2], [5, 2, 2, 2], [3, 2, 2, 3], [3, 2, 3, 2],
                [2, 2, 3, 3]]

num_runs = 10
f = open('speed_dd_asc.txt', 'w+')
for tconfig in test_configs:
    ma, mb, na, nb = tconfig
    inputs_a, inputs_b, outputs_a, outputs_b = range(ma), range(mb), range(na), range(nb)

    # Write Vertices to a File
    vertices = get_deterministic_behaviors_two_party(inputs_a, inputs_b, outputs_a, outputs_b)
    avg_time = 0
    for _ in range(num_runs):
        np.random.shuffle(vertices)
        # shuffle the vertices randomly
        vertex_file = 'vertices/{}{}{}{}.ext'.format(ma, mb, na, nb)
        print('computing ' + vertex_file)
        write_known_vertices(vertices, vertex_file)
        # Run panda with given options
        print('Starting Panda')
        start = time.time()
        os.system("panda -t 1 -s lex_asc -m dd {}".format(vertex_file))
        exc_time = time.time() - start
        print('took: {} s'.format(exc_time))
        avg_time += exc_time

    avg_time = avg_time / num_runs
    f.write('{}{}{}{} : {} ms \n'.format(ma, mb, na, nb, avg_time * 1000))
f.close()
