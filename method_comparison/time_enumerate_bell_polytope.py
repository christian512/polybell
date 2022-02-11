""" Script for enumerating the facet-classes of a two-party Bell polytope """

import argparse
from gap_program_templates import *
from linearbell.utils import get_relabelling_generators, get_deterministic_behaviors_two_party, get_relabels_dets
from linearbell.gap_helper import relabels_dets_to_disjoint_cycles
from linearbell.panda_helper import write_known_vertices
import numpy as np
import subprocess
import time
import os
import signal

# parse the inputs
parser = argparse.ArgumentParser()
parser.add_argument(dest='ma', help="number of inputs for ALICE")
parser.add_argument(dest='mb', help='number of inputs for BOB')
parser.add_argument(dest='na', help='number of outputs for ALICE')
parser.add_argument(dest='nb', help='number of outputs for BOB')
parser.add_argument('-e', type=str, default='stabilizer', dest='equiv_method',
                    choices=('stabilizer', 'local', 'global'),
                    help='Choose wich type of equivalence test to perform. Default: stabilizer')
parser.add_argument('-r', type=int, default=0, dest='recursion_depth', help='Sets the recursion depth. Default: 0')
parser.add_argument('-s', action='store_true', dest='sampling_method', help='Flag to activate the sampling method')
parser.add_argument('-dd', action='store_true', dest='dd_method',
                    help='Flag to activate the DoubleDescription Method instead of Adjacency Decomposition.')
args = parser.parse_args()

# set the scenario
ma, mb, na, nb = int(args.ma), int(args.mb), int(args.na), int(args.nb)
inputs_a, inputs_b, outputs_a, outputs_b = range(ma), range(mb), range(na), range(nb)

# Write Vertices to a File
vertices = get_deterministic_behaviors_two_party(inputs_a, inputs_b, outputs_a, outputs_b)
vertices = vertices[np.lexsort(np.rot90(vertices))]
randa_filename = 'vertices/{}{}{}{}.ext'.format(ma, mb, na, nb)
write_known_vertices(vertices, randa_filename)

# Create the generators of the symmetry group
generators = get_relabelling_generators(inputs_a, inputs_b, outputs_a, outputs_b)
automorphisms = get_relabels_dets(vertices, generators, show_progress=0)
disjoint_cycles = relabels_dets_to_disjoint_cycles(automorphisms)

# Identify the equivalence test to perform and select the template script accordingly
script = ''
if args.equiv_method == 'stabilizer':
    print('Generating GAP program with stabilizer equivalence test')
    script = template_stabilizer_equiv_test % disjoint_cycles
if args.equiv_method == 'global':
    print('Generating GAP program with global equivalence test')
    script = template_global_equiv_test % disjoint_cycles
if args.equiv_method == 'local':
    print('Generating GAP program with local equivalence test')
    script = template_local_equiv_test % disjoint_cycles
gap_filename = 'gap_scripts/{}{}{}{}_{}.g'.format(ma, mb, na, nb, args.equiv_method)
f = open(gap_filename, 'w+')
f.write(script)
f.close()

# setup the FIFO PIPES
try:
    os.mkfifo('togap.pipe')
    os.mkfifo('fromgap.pipe')
except FileExistsError:
    pass

num_runs = 10
times = []
for i in range(num_runs):
    print('Starting run {} / {}'.format(i+1, num_runs))
    # Run the GAP Program using Popen.
    cmd = 'gap --quitonbreak ' + gap_filename
    gap_proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    # wait for GAP to finish initial launch
    time.sleep(5)

    # define behavior for receiving a SIGINT (e.g. pressing ctrl+c)
    def handler(signum, frame):
        print("\n Ctrl-c was pressed. Stopping GAP and RANDA.")
        os.killpg(os.getpgid(gap_proc.pid), signal.SIGTERM)
        exit(1)
    signal.signal(signal.SIGINT, handler)

    # Start the RANDA program
    sampling_flag = 0
    if args.sampling_method:
        sampling_flag = 1
    cmd = 'randa -t 1 -r ' + str(args.recursion_depth) + ' -p ' + str(
        sampling_flag) + ' ' + randa_filename


    if args.dd_method:
        cmd = 'randa -t 1 -m dd ' + randa_filename

    # execture command
    start = time.time()
    randa_proc = subprocess.run(cmd, shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    run_time = time.time() - start
    times.append(run_time)
    print('finished after {} ms'.format(run_time * 1000))

    # Wait for GAP to finish
    gap_proc.wait()

times = np.array(times)

print('Finished runs for {}{}{}{}'.format(ma,mb,na,nb))
print('Average run-time over {} runs: {} ms'.format(num_runs, np.average(times)*1000))
print('Standard deviation of runtime: {} ms'.format(np.std(times)*1000))