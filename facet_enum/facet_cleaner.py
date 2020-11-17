from linearbell.utils import get_deterministic_behaviors, get_allowed_relabellings, check_equiv_bell, get_relabels_dets
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(dest='ma', help="number of inputs for ALICE")
parser.add_argument(dest='mb', help='number of inputs for BOB')
parser.add_argument(dest='n', help='number of outputs')
args = parser.parse_args()
ma = int(args.ma)
mb = int(args.mb)
n = int(args.n)

# set inputs / outputs
inputs_a = range(ma)
inputs_b = range(mb)
outputs = range(n)

# get deterministic points
dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs)

# get allowed relabellings
allowed_relabellings = get_allowed_relabellings(inputs_a, inputs_b, outputs, outputs)
# get relabellings for deterministic points
file_relabels = '../data/relabels_dets/{}{}{}{}.gz'.format(ma, mb, n, n)
try:
    relabels_dets = np.loadtxt(file_relabels, dtype=float).astype(int)
except IOError:
    print('Have to calculate the possible relabels before actual start')
    relabels_dets = get_relabels_dets(dets, allowed_relabellings)
    np.savetxt(file_relabels, relabels_dets)
print('Created all possible relabels')
# load facets from file
facets_file = '../data/facets/{}{}{}{}.txt'.format(ma, mb, n, n)
facets = np.loadtxt(facets_file)

# Round the facets
facets = np.round(facets, decimals=3)

# indices of facets to delete
del_facets = []
print('start iteration through facets')
# iterate through the facets
for i in range(facets.shape[0]):
    # if this facet can already be deleted -> continue
    print('facet: {} / {}'.format(i, facets.shape[0]))
    if i in del_facets: continue
    for j in range(i + 1, facets.shape[0]):
        # if this facet can already be deleted -> continue
        if j in del_facets: continue
        # check if the two facets are equivalent
        if check_equiv_bell(facets[i], facets[j], relabels_dets, dets, tol=1e-4):
            del_facets.append(j)
# store new facets
new_facets = np.delete(facets, del_facets, axis=0)
new_facets_file = '../data/facets/{}{}{}{}_filtered.txt'.format(ma, mb, n, n)
# np.savetxt(facets_file, new_facets)
print('Number of different facets: {}'.format(new_facets.shape[0]))
