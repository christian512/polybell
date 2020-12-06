from linearbell.utils import *
from linearbell.lrs_helper import polyhedra_h_representation, run_redund, run_lrs_h_repr

# number of inputs and outputs
ma = 3
mb = 3
n = 3

# set inputs / outputs
inputs_a = range(ma)
inputs_b = range(mb)
outputs = range(n)
outputs_wo_failure = range(n - 1)

# get deterministic points
dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs)

# setup pr box
configs = get_configs(inputs_a, inputs_b, outputs, outputs)
eta = 4 / (4 + len(inputs_a))
pr_box = [general_pr_box_extended(a, b, x, y, eta, outputs_wo_failure) for (a, b, x, y) in configs]
pr_box = np.array(pr_box)

# check which dets can actually contribute
epsilon = 0.01
dets_contrib = []
dets_contrib_idx = []
for i in range(dets.shape[0]):
    d = dets[i]
    pr_tmp = pr_box - epsilon * d
    try:
        lws = find_local_weight_primal(pr_tmp, dets)
        lw = np.sum(lws)
        if lw <= 1.0:
            dets_contrib_idx.append(i)
            dets_contrib.append(d)
    except:
        pass

print('dets that can contribute: {}'.format(len(dets_contrib_idx)))

# find out how many are acutally linear independent
dets_contrib = np.array(dets_contrib)
rank = np.linalg.matrix_rank(dets_contrib)
print('rank of contributing dets: {}'.format(rank))

# try to write finding all combinations as a vertex enumeration problem
lhs = np.transpose(dets_contrib)
rhs = np.copy(pr_box)

lhs = np.r_[lhs, [np.ones(lhs.shape[1])]]
rhs = np.r_[rhs, [1.0]]

# set these to be linear
lins = list(np.arange(lhs.shape[0]))

# set positivity of weights
lhs = np.r_[lhs, -1.0 * np.eye(lhs.shape[1])]
rhs = np.r_[rhs, np.zeros(lhs.shape[1])]

# setup the polyhedra file
polyhedra_h_representation(lhs, rhs, linearities=lins, file='input.ine')
run_redund('input.ine', 'input_redund.ine')
vertices, rays = run_lrs_h_repr('input_redund.ine', 'output.ext', nproc=30)

# store the factors of contribution for every deterministic point
factors = np.zeros((dets.shape[0], vertices.shape[1]))
for j in range(vertices.shape[0]):
    v = vertices[j]
    for i in range(dets.shape[0]):
        if i in dets_contrib_idx:
            idx = dets_contrib_idx.index(i)
            factors[j, i] = v[idx]

# store the factors of the deterministic points
file = '../data/pr_box_contrib_det/{}{}{}{}.gz'.format(ma, mb, n, n)
np.savetxt(file, factors)
