import multiprocessing as mp
import numpy as np
from linearbell.utils import find_local_weight, facet_inequality_check, check_equiv_bell
from linearbell.utils import extremal_ns_binary_vertices, get_deterministic_behaviors, get_allowed_relabellings, \
    get_relabels_dets
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(dest='ma', help="number of inputs for ALICE")
parser.add_argument(dest='mb', help='number of inputs for BOB')
parser.add_argument(dest='n', help='number of outputs')
parser.add_argument(dest='num_cpu', help='number of cpu corse to be used')
parser.add_argument(dest='num_eps', help='number of epsilons to test')
args = parser.parse_args()
ma = int(args.ma)
mb = int(args.mb)
n = int(args.n)
num_cpu = int(args.num_cpu)
num_eps = int(args.num_eps)

# set inputs / outputs
inputs_a = range(ma)
inputs_b = range(mb)
outputs = range(n)
# get extremal points
extremals = extremal_ns_binary_vertices(inputs_a, inputs_b, outputs)

# get deterministic points
dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs)

# get allowed relabellings
allowed_relabellings = get_allowed_relabellings(inputs_a, inputs_b, outputs, outputs)

# get relabellings for deterministic points
file_relabels = '../data/relabels_dets/{}{}{}{}.gz'.format(ma, mb, n, n)
try:
    relabels_dets = np.loadtxt(file_relabels)
except IOError:
    print('Have to calculate the possible relabels before actual start')
    relabels_dets = get_relabels_dets(dets, allowed_relabellings)
    np.savetxt(file_relabels, relabels_dets)

# set the epsilons that we want to use
epsilons = np.linspace(0.5, 2 / 3 - 1e-2, num=num_eps)

# tolerance
tol = 1e-4

# array of facets
facets_folder = '../data/facets/{}{}{}{}/'.format(ma, mb, n, n)


# enumerate extremals
def find_facets_for_extremal(idx):
    # array for facets
    facets = []
    # file for the facets corresponding to this extremal
    facets_file = facets_folder + '{}.txt'.format(idx)
    bell_expression = find_local_weight(extremals[idx], dets)
    assert np.abs(bell_expression @ extremals[idx]) < tol, 'local weight of extremal is not zero. Problem with solver'
    # get the equalizing behaviors
    is_facet, bell_expression, eq_dets = facet_inequality_check(dets, bell_expression, ma, mb, n, tol)
    if is_facet:
        # check if it's equivalent to any other bell expression already found
        for j in range(len(facets) - 1, -1, -1):
            if check_equiv_bell(bell_expression, facets[j], relabels_dets, dets):
                is_facet = False
                break
    # if it's a new facet, append it to the facet array
    if is_facet:
        facets.append(bell_expression)
        print('num facets: {}'.format(len(facets)))
    count = 0
    # iterate through epsilons
    for m in range(num_eps):
        epsilon = epsilons[m]
        # TODO: This should be standard
        eq_dets = np.copy(dets)
        # define the new behavior
        for j in range(eq_dets.shape[0]):
            for k in range(eq_dets.shape[0]):
                print(' {} / {}  eq_dets || {} / {} extremals'.format(count, num_eps * eq_dets.shape[0] ** 2, idx,
                                                                      len(extremals)))
                count += 1
                # form new behavior
                e_new = (1 - 3 * epsilon / 2) * extremals[idx] + epsilon * eq_dets[j] + epsilon / 2 * eq_dets[k]
                # check again if we can find a facet
                bell_expression = find_local_weight(e_new, dets)
                is_facet, bell_expression, _ = facet_inequality_check(dets, bell_expression, ma, mb, n, tol)
                if not is_facet: continue
                for l in range(len(facets) - 1, -1, -1):
                    if check_equiv_bell(bell_expression, facets[l], relabels_dets, dets, tol=tol):
                        is_facet = False
                        break
                if not is_facet: continue
                facets.append(bell_expression)
                print('num facets: {}'.format(len(facets)))
    facets = np.array(facets)
    np.savetxt(facets_file, facets)


# setup parallel environment
manager = mp.Manager()
pool = mp.Pool(num_cpu)
jobs = []
for i in range(len(extremals)):
    job = pool.apply_async(find_facets_for_extremal, (i,))
    jobs.append(job)

print('Created all Jobs, start processing!')

for job in jobs:
    job.get()
