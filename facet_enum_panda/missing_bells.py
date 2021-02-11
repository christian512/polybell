""" Here we analyse which Bell expressions we're missing out with the partial recursion strategy. """

from linearbell.utils import get_configs, get_deterministic_behaviors, equiv_check_adjacency_testing
from linearbell.representations import get_configs_mat, transform_vec_to_mat, transform_mat_to_vec
import numpy as np
from scipy.io import loadmat

# setup the configuration
inputs = range(4)
outputs = range(2)
configs_mat = get_configs_mat(inputs, inputs, outputs, outputs)
configs = get_configs(inputs, inputs, outputs, outputs)
dets = get_deterministic_behaviors(inputs, inputs, outputs)
relabels = np.loadtxt('../data/relabels/4422.gz').astype(int)

# load inequalities found by partial recursion
facets_tmp = np.loadtxt('../data/facets_panda/4422.txt')
facets_partial = []
for facet in facets_tmp:
    # Last entry gives the negative of the rhs, we want the RHS to be one
    # Thus divide by (-1 / last entry)
    # if the last entry is zero, we have to add 1 / (num_outputs ** 2) to each entry
    divisor = facet[-1]
    if divisor == 0:
        f = facet[:-1] + 1 / (len(outputs) ** 2)
    else:
        assert divisor < 0
        f = facet[:-1] / np.abs(divisor)
    facets_partial.append(f)
facets_partial = np.array(facets_partial)

# load inequalities found by TOM
facets_tmp = loadmat('../data/facets_known/Inequalities_Only_4422_Raw.mat')
facets_known = []
for k in facets_tmp.keys():
    mat = facets_tmp[k]
    f = transform_mat_to_vec(configs_mat, configs, mat)
    mat_copy = transform_vec_to_mat(configs, configs_mat, f)
    assert np.all(mat == mat_copy)
    # Transformation to get local bound 1
    facets_known.append(-1.0 * f + 2 / (len(outputs)**2))
facets_known = np.array(facets_known)

# Check the equivalence of all the known inequalities to the ones found by PANDA
facets_not_found = []
for i in range(facets_known.shape[0]):
    print('progress: {} / {}'.format(i, facets_known.shape[0]))
    equiv = False
    for facet_found in facets_partial:
        if equiv_check_adjacency_testing(facets_known[i], facet_found, relabels, dets) and not equiv:
            equiv = True
    if not equiv:
        facets_not_found.append(facets_known[i])
        print('Number of facets not found: ', len(facets_not_found))

facets_not_found = np.array(facets_not_found)
print('Number of facets not found by partial recursive algorithm: ', facets_not_found.shape[0])
