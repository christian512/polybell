from linearbell.utils import get_deterministic_behaviors, find_local_weight_primal, get_configs, general_pr_box_extended
from linearbell.lrs_helper import polyhedra_h_representation
import numpy as np

# set inputs / outputs
inputs_a = range(2)
inputs_b = range(2)
outputs = range(3)
outputs_wo_failure = range(2)

# get deterministic points
dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs)

# setup pr box
configs = get_configs(inputs_a, inputs_b, outputs, outputs)
eta = 4 / (4 + len(inputs_a))
pr_box = [general_pr_box_extended(a, b, x, y, eta, outputs_wo_failure) for (a, b, x, y) in configs]
pr_box = np.array(pr_box)

rank = np.linalg.matrix_rank(dets)
print('rank: {}'.format(rank))

# choose independent deterministics
ind_dets = [dets[0]]
tmp_rank = np.linalg.matrix_rank(ind_dets)

for d in dets[1:]:
    tmp_dets = np.copy(np.array(ind_dets))
    tmp_dets = np.r_[tmp_dets, [d]]
    if np.linalg.matrix_rank(tmp_dets) > tmp_rank:
        ind_dets.append(d)
        tmp_rank = np.linalg.matrix_rank(tmp_dets)
ind_dets = np.array(ind_dets)
assert np.linalg.matrix_rank(ind_dets) == rank

# get the local weight using these independent deterministic points
weights = find_local_weight_primal(pr_box, dets)
print('number of nonzero weights: {}'.format(len(weights[weights > 1e-6])))
# set equalizing dets
eq_dets = dets[weights > 1e-6]

# write file for running the vertex enumeration
lhs = np.copy(dets)
lhs = np.r_[lhs, [pr_box]]
rhs = np.ones(lhs.shape[0])
lins = list(np.arange(weights.shape[0])[weights > 1e-6])
lins.append(lhs.shape[0] - 1)

repr = polyhedra_h_representation(lhs, rhs, lins, file='input.ine')
print(repr)


