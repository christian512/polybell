from linearbell.utils import extremal_ns_binary_vertices, get_deterministic_behaviors, get_allowed_relabellings
import numpy as np
from linearbell.utils import find_local_weight, facet_inequality_check, check_equiv_bell
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(dest='ma', help="number of inputs for ALICE")
parser.add_argument(dest='mb', help='number of inputs for BOB')
parser.add_argument(dest='n', help='number of outputs')
parser.add_argument(dest='num_cpu', help='number of cpu cores to be used')
args = parser.parse_args()
ma = int(args.ma)
mb = int(args.mb)
n = int(args.n)
num_cpu = int(args.num_cpu)

# set inputs / outputs
inputs_a = range(ma)
inputs_b = range(mb)
outputs = range(n)

# files where to store the facets and equalizing deterministics for later
facets_file = '../data/facets/{}{}{}{}.txt'.format(len(inputs_a), len(inputs_b), len(outputs), len(outputs))
eq_dets_file = '../data/equalizing_deterministics/{}{}{}{}.txt'.format(len(inputs_a), len(inputs_b), len(outputs),
                                                                       len(outputs))
# get extremal points
extremals = extremal_ns_binary_vertices(inputs_a, inputs_b, outputs)

# get deterministic points
dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs)

# get allowed relabellings
allowed_relabellings = get_allowed_relabellings(inputs_a, inputs_b, outputs, outputs)

# open files to create them
f = open(facets_file, 'w+')
f.close()
f = open(eq_dets_file, 'w+')
f.close()


def get_string_from_numpy_array(arr):
    s = ""
    for x in arr:
        s += "%.18e" % x
        s += " "
    s = s[:-1]
    s += "\n"
    return s


def get_facet(extremal, facet_q, eqdets_q):
    # read current facets from file -> does not work with np.loadtxt
    lines = open(facets_file).readlines()
    curr_facets = []
    for line in lines:
        f = [float(x) for x in line.split(' ')]
        curr_facets.append(f)
    curr_facets = np.array(curr_facets)
    # calculate local weight
    bell_exp = find_local_weight(extremal, dets)
    # check if there is a facet
    isf, bell, eq_dets = facet_inequality_check(dets, bell_exp, len(inputs_a), len(inputs_b), len(outputs))
    if isf:
        for i in range(len(curr_facets) - 1, -1, -1):
            if check_equiv_bell(bell, curr_facets[i], allowed_relabellings, dets):
                isf = False
                break
    # write equalizing dets to equalizing dets file
    if eqdets_q:
        for eq_det in eq_dets:
            s = get_string_from_numpy_array(eq_det)
            eqdets_q.put(s)
        eqdets_q.put('+++++ \n')
    # if it's a facet, write facet to facet writer
    if isf:
        s = get_string_from_numpy_array(bell)
        # append string of facet to queue
        facet_q.put(s)
    return True


def facet_writer(queue):
    print('Facet writer started')
    new_facet = None
    # load current facets
    with open(facets_file, 'a+') as out:
        # infinity loop only stopped when kill is read
        while True:
            m = queue.get()
            # stop when kill is received
            if m == 'kill':
                print('facet writer killed')
                break
            # check if the facet that should be written was the last facet found
            # maybe the facet worker did not read it from the file yet
            # Maybe you can add here multiple new_facets if you have more workers available
            if not new_facet is None:
                # get the facet from the string
                f1 = np.array([[float(r) for r in m[:-2].split(' ')]])
                # check if the received facet is equal to the last newly found facet
                if not check_equiv_bell(new_facet, f1, allowed_relabellings, dets):
                    # if not
                    out.write(m)
                    out.flush()
                    new_facet = f1
                    print('wrote facet: {}'.format(new_facet))
            # if there is no new facet yet found
            if new_facet is None:
                out.write(m)
                out.flush()
                # new_facet = np.array([float(r) for r in m[:-2].split(' ')])
                print('wrote facet')
    return True


def eqdets_writer(queue):
    print('Extremal Writer Started')
    with open(eq_dets_file, 'a+') as out:
        while True:
            m = queue.get()
            if m == 'kill':
                print('extremal writer killed')
                break
            out.write(m)
            out.flush()
    return True


import multiprocessing as mp

# setup multiprocessing pool, queues and manager
manager = mp.Manager()
facet_q = manager.Queue()
extremal_q = manager.Queue()
pool = mp.Pool(num_cpu)
# setup watchers
watcher_facets = pool.apply_async(facet_writer, (facet_q,))
watcher_eqdets = pool.apply_async(eqdets_writer, (extremal_q,))

# start workers
jobs = []
for e in extremals:
    job = pool.apply_async(get_facet, (e, facet_q, extremal_q))
    jobs.append(job)
# run the jobs

for i, job in enumerate(jobs):
    if i % 20 == 0:
        print("{} / {}".format(i, len(jobs)))
    job.get()

facet_q.put('kill')
extremal_q.put('kill')
watcher_facets.get()
pool.close()
pool.join()
