import argparse
from linearbell.utils import get_deterministic_behaviors, get_configs, general_pr_box_extended, \
    get_parametrisation_configs, parametrise_behavior
from pypoman import compute_polytope_vertices, compute_polytope_halfspaces
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument(dest='ma', help="number of inputs for ALICE")
parser.add_argument(dest='mb', help='number of inputs for BOB')
parser.add_argument(dest='n', help='number of outputs')
args = parser.parse_args()
ma = int(args.ma)
mb = int(args.mb)
n = int(args.n)
assert n == 3, 'This program currently only works for n = 3 '

# set inputs / outputs
inputs_a = range(ma)
inputs_b = range(mb)
outputs = range(n)
outputs_wo_failure = range(n - 1)

# get deterministic points
dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs)

# get the configurations
configs = get_configs(inputs_a, inputs_b, outputs, outputs)
# setup eta and epsilon
assert ma == mb, 'You must choose the same input dimension for Alice and Bob'
epsilon = 0.2
tol = 1e-3
eta = 4 / (ma + 4)
print('eta: {} || epsilon: {}'.format(eta, epsilon))

# define the pr-boxes
pr_eff_high = np.array(
    [general_pr_box_extended(a, b, x, y, eta + epsilon, outputs_wo_failure) for a, b, x, y in configs])

pr_eff_low = np.array(
    [general_pr_box_extended(a, b, x, y, eta - epsilon, outputs_wo_failure) for a, b, x, y in configs])

# Parametrise to get lower dimensional representations
configs_param = get_parametrisation_configs(inputs_a, inputs_b, outputs, outputs)
# get parametrisations of deterministic points
dets_param = []
for d in dets:
    d_param = parametrise_behavior(d, configs, configs_param, inputs_a, inputs_b, outputs, outputs)
    dets_param.append(d_param)
dets_param = np.array(dets_param)
dets = np.copy(dets_param)
# parametristations of pr boxes
pr_eff_high = parametrise_behavior(pr_eff_high, configs, configs_param, inputs_a, inputs_b, outputs, outputs)
pr_eff_low = parametrise_behavior(pr_eff_low, configs, configs_param, inputs_a, inputs_b, outputs, outputs)

# setup the constraints: A * x <= b

# determinstics should be local
lhs_ineq = np.copy(dets)
rhs_ineq = np.ones(dets.shape[0])

# positive entries
lhs_ineq = np.r_[lhs_ineq, -1.0 * np.eye(dets.shape[1])]
rhs_ineq = np.r_[rhs_ineq, np.zeros(dets.shape[1])]

# low efficiency pr box should be local
lhs_ineq = np.r_[lhs_ineq, [pr_eff_low]]
rhs_ineq = np.r_[rhs_ineq, [1.0 - tol]]

# high efficiency pr box shoud be non local
lhs_ineq = np.r_[lhs_ineq, [-1.0 * pr_eff_high]]
rhs_ineq = np.r_[rhs_ineq, [-1.0 + tol]]

# compute vertices of the polytope
print('Start computation of vertices')
vertices = compute_polytope_vertices(lhs_ineq, rhs_ineq)

# store the vertices
vertices = np.array(vertices)
file = '../data/vertex_enum_pr_box/{}{}{}{}.txt'.format(ma, mb, n, n)
np.savetxt(file, vertices)
print('done')
print(vertices)
