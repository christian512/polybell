import numpy as np
from linearbell.utils import find_local_weight, facet_inequality_check, check_equiv_bell
from linearbell.utils import extremal_ns_binary_vertices, get_deterministic_behaviors, get_allowed_relabellings
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
# set the epsilons that we want to use
epsilons = np.linspace(1 / 3, 2 / 3, num=num_eps)

# tolerance
tol = 1e-6

# array of facets
facets = []
facets_file = '../data/facets/{}{}{}{}.txt'.format(ma, mb, n, n)

# enumerate extremals
for i in range(len(extremals)):
    bell_expression = find_local_weight(extremals[i], dets)
    assert np.abs(bell_expression @ extremals[i]) < tol, 'local weight of extremal is not zero. Problem with solver'
    # get the equalizing behaviors
    is_facet, bell_expression, eq_dets = facet_inequality_check(dets, bell_expression, ma, mb, n, tol)
    if is_facet:
        # check if it's equivalent to any other bell expression already found
        for j in range(len(facets) - 1, -1, -1):
            if check_equiv_bell(bell_expression, facets[j], allowed_relabellings, dets):
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

        # define the new behavior
        for j in range(eq_dets.shape[0]):
            for k in range(eq_dets.shape[0]):
                print(' {} / {}  eq_dets || {} / {} extremals'.format(count, num_eps * eq_dets.shape[0] ** 2, i,
                                                                      len(extremals)))
                count += 1
                # form new behavior
                e_new = (1 - 3 * epsilon / 2) * extremals[i] + epsilon * eq_dets[j] + epsilon / 2 * eq_dets[k]
                # check again if we can find a facet
                bell_expression = find_local_weight(e_new, dets)
                is_facet, bell_expression, _ = facet_inequality_check(dets, bell_expression, ma, mb, n, tol)
                if not is_facet: continue
                for l in range(len(facets) - 1, -1, -1):
                    if check_equiv_bell(bell_expression, facets[l], allowed_relabellings, dets):
                        is_facet = False
                        break
                if not is_facet: continue
                facets.append(bell_expression)
                print('num facets: {}'.format(len(facets)))
facets = np.array(facets)
np.savetxt(facets_file, facets)
