""" Here we analyse which Bell expressions we're missing out with the partial recursion strategy. """

from polybell.utils import get_configs, get_deterministic_behaviors, equiv_check_adjacency_testing
from polybell.representations import get_configs_mat, transform_vec_to_mat, transform_mat_to_vec
import numpy as np
from polybell.panda_helper import read_inequalities
from scipy.io import loadmat

# setup the configuration
inputs = range(4)
outputs = range(2)
configs_mat = get_configs_mat(inputs, inputs, outputs, outputs)
configs = get_configs(inputs, inputs, outputs, outputs)
dets = get_deterministic_behaviors(inputs, inputs, outputs)


# load inequalities found by partial recursion given as b * p + a <= 0
# facets_tmp = np.loadtxt('../data/facets_panda/4422.txt')
facets_tmp = read_inequalities('../randa_testing/4422.txt')
facets_partial = []
for facet in facets_tmp:
    # transform to b * p >= 1
    for d in dets:
        assert facet[:-1] @ d <= -1.0 * facet[-1], 'Assumption not correct'
    facet = -1 * facet
    divisor = -1 * facet[-1]
    f = facet[:-1] + (1 - divisor) / (len(outputs) ** 2)
    for d in dets:
        assert d @ f >= 1, 'Error in rescaling of inequality: d @ f = {} \n d = {} \n f = {} \n divisor = {}'.format(
            d @ f, d, f, divisor)

    facets_partial.append(f)
facets_partial = np.array(facets_partial)
# load inequalities found by TOM, Given in form Tr(b*p) >= 1
facets_tmp = loadmat('../data/facets/Inequalities_Only_4422_Raw.mat')
facets_known = []
for k in facets_tmp.keys():
    mat = facets_tmp[k]
    f = transform_mat_to_vec(configs_mat, configs, mat)
    mat_copy = transform_vec_to_mat(configs, configs_mat, f)
    assert np.all(mat == mat_copy)
    # Transformation to get local bound 1
    facets_known.append(f)
    # facets_known.append(-1.0 * f + 2 / (len(outputs)**2))
facets_known = np.array(facets_known)

relabels = np.loadtxt('../data/relabels/4422.gz').astype(int)
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

np.savetxt('../data/facets/heuristic_adm_4422_notfound.txt', facets_not_found)