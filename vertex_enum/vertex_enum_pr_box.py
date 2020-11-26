import argparse
from linearbell.utils import get_deterministic_behaviors, get_configs, general_pr_box_extended, \
    get_parametrisation_configs, parametrise_behavior, check_equiv_bell
from lrs_helper import polytope_v_representation, run_lrs_v_repr, polyhedra_h_representation, run_lrs_h_repr
import numpy as np
from itertools import product

parser = argparse.ArgumentParser()
parser.add_argument(dest='ma', help="number of inputs for ALICE")
parser.add_argument(dest='mb', help='number of inputs for BOB')
parser.add_argument(dest='n', help='number of outputs')
parser.add_argument(dest='nproc', help='Number of processes')
parser.add_argument(dest='outfile', help='Output file for lrs computation')
args = parser.parse_args()
ma = int(args.ma)
mb = int(args.mb)
n = int(args.n)
nproc = int(args.nproc)
outfile = str(args.outfile)

# set outfile with path
outfile = '../data/vertex_enum_pr_box/' + outfile

# set inputs / outputs
inputs_a = range(ma)
inputs_b = range(mb)
outputs = range(n)
outputs_wo_failure = range(n - 1)

# values for pr boxes
epsilon = 1e-2
tol = 1e-3
assert ma == mb
eta = 4 / (ma + 4)

# get deterministic points
dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs)

# Define new origin, such that the 0 origin is inside of the polytope (for dual representation)
p_origin = np.sum(dets, axis=0) / dets.shape[0]

# shift the origin of the deterministic points
dets = dets - p_origin

# get configurations
configs = get_configs(inputs_a, inputs_b, outputs, outputs)

# define the pr box for high and low efficiency
pr_eff_high = np.array(
    [general_pr_box_extended(a, b, x, y, eta + epsilon, outputs_wo_failure) for a, b, x, y in configs])
pr_eff_low = np.array(
    [general_pr_box_extended(a, b, x, y, eta - epsilon, outputs_wo_failure) for a, b, x, y in configs])
pr_eff_exact = np.array(
    [general_pr_box_extended(a, b, x, y, eta, outputs_wo_failure) for a, b, x, y in configs])

# Parametrise to get lower dimensional representations
configs_param = get_parametrisation_configs(inputs_a, inputs_b, outputs, outputs)

# Parametrise deterministic points
dets_param = []
for d in dets:
    d_param = parametrise_behavior(d, configs, configs_param, inputs_a, inputs_b, outputs, outputs)
    dets_param.append(d_param)
dets_param = np.array(dets_param)
dets = np.copy(dets_param)

# parametrise pr boxes
pr_eff_high = parametrise_behavior(pr_eff_high, configs, configs_param, inputs_a, inputs_b, outputs, outputs)
pr_eff_low = parametrise_behavior(pr_eff_low, configs, configs_param, inputs_a, inputs_b, outputs, outputs)
pr_eff_exact = parametrise_behavior(pr_eff_exact, configs, configs_param, inputs_a, inputs_b, outputs, outputs)

# setup the constraints: dets @ bell <= 1
lhs_ineq = np.copy(dets)
rhs_ineq = np.ones(dets.shape[0])

# constraint for low efficiency pr box
# lhs_ineq = np.r_[lhs_ineq, [pr_eff_low]]
# rhs_ineq = np.r_[rhs_ineq, [1.0]]

# constraint for high efficiency pr box
lhs_ineq = np.r_[lhs_ineq, [-1.0 * pr_eff_high]]
rhs_ineq = np.r_[rhs_ineq, [-1.0]]

# constraint for exact pr box
# lhs_ineq = np.r_[lhs_ineq, [pr_eff_exact]]
# rhs_ineq = np.r_[rhs_ineq, [1.0]]

# set the inequalities that should be equalities
linearities = []

# get the h representation of the polyhedra
hrepr = polyhedra_h_representation(lhs_ineq, rhs_ineq,
                                   linearities=linearities,
                                   file='input_h.ine')

# run the h representation of the polyhedra
vertices, rays = run_lrs_h_repr('input_h.ine', output_file=outfile, nproc=nproc)

# print the numbers
print('number of vertices: {}'.format(vertices.shape[0]))
print('number of rays: {}'.format(rays.shape[0]))
