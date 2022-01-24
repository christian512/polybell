""" Writing the vertices to a file and printing the generators of the automorphism group for a bell scenario """

from linearbell.utils import get_relabelling_generators, get_deterministic_behaviors_two_party, get_relabels_dets
from linearbell.gap_helper import relabels_dets_to_disjoint_cycles
from linearbell.panda_helper import write_known_vertices
import argparse
import numpy as np
import time
import os

# parse the inputs
# parser = argparse.ArgumentParser()
# parser.add_argument(dest='ma', help="number of inputs for ALICE")
# parser.add_argument(dest='mb', help='number of inputs for BOB')
# parser.add_argument(dest='na', help='number of outputs for ALICE')
# parser.add_argument(dest='nb', help='number of outputs for BOB')
# args = parser.parse_args()

# set the scenario

# test_configs = [[2,2,2,2], [3,2,2,2], [2,2,3,2], [3,3,2,2], [5,2,2,2], [4,3,2,2], [3,2,2,3], [5,3,2,2]]
test_configs = [[2,2,2,2], [3,2,2,2], [2,2,3,2], [3,3,2,2], [5,2,2,2], [4,3,2,2], [3,2,2,3]]
num_runs = 1
f = open('speed_dd_random.txt', 'w+')
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
        os.system("panda -t 1 -m dd {}".format(vertex_file))
        exc_time = time.time() - start
        print('took: {} s'.format(exc_time))
        avg_time += exc_time

    avg_time = avg_time / num_runs
    f.write('{}{}{}{} : {} \n'.format(ma,mb,na,nb, avg_time))
f.close()