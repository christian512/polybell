from linearbell.utils import *
from linearbell.lrs_helper import polyhedra_h_representation, run_redund, run_lrs_h_repr

# set inputs / outputs
inputs_a = range(3)
inputs_b = range(3)
outputs = range(3)
outputs_wo_failure = range(2)

# get deterministic points
dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs)

# setup pr box
configs = get_configs(inputs_a, inputs_b, outputs, outputs)
eta = 4 / (4 + len(inputs_a))
pr_box = [general_pr_box_extended(a, b, x, y, eta, outputs_wo_failure) for (a, b, x, y) in configs]
pr_box = np.array(pr_box)

# check how many deterministic points we have in the decomposition to local points
lws = find_local_weight_primal(pr_box, dets)
decomp_counter = len(lws[lws > 0.0])
print('dets in decomp: {}'.format(decomp_counter))

# check which dets can actually contribute
epsilon = 0.01
no_solution_counter = 0
locals_counter = 0
non_local_counter = 0
dets_contrib = []
for d in dets:
    pr_tmp = pr_box - epsilon * d
    try:
        lws = find_local_weight_primal(pr_tmp, dets)
        lw = np.sum(lws)
        if lw <= 1.0:
            locals_counter += 1
            dets_contrib.append(d)
        else:
            non_local_counter += 1
    except:
        no_solution_counter += 1

print('dets that can contribute: {}'.format(locals_counter))

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
vertices, rays = run_lrs_h_repr('input_redund.ine', 'output.ext')





