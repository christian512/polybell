import argparse
from polybell.utils import get_deterministic_behaviors, get_configs, general_pr_box_extended, \
    get_parametrisation_configs, parametrise_behavior, check_equiv_bell
from lrs_helper import polytope_v_representation, run_lrs_v_repr, polyhedra_h_representation, run_lrs_h_repr
import numpy as np
from itertools import product

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

# Define new origin, such that the 0 origin is inside of the polytope (for dual representation)
p_origin = np.sum(dets, axis=0) / dets.shape[0]

# shift the origin of the deterministic points
dets = dets - p_origin

# get configurations
configs = get_configs(inputs_a, inputs_b, outputs, outputs)

# write v-representation file
file = 'input_v.ext'
vrepr = polytope_v_representation(dets, file=file)

# run the file
facets, lins, lins_values = run_lrs_v_repr(file)

# print the numbers
print('number facets: {}'.format(facets.shape[0]))
print('number linearities: {}'.format(lins.shape[0]))

# now we want to check that we find these facets as vertices of the dual

# setup the constraints

# det @ bell <= 1
lhs_ineq = np.copy(dets)
rhs_ineq = np.ones(dets.shape[0])

# equalizing one determinstic point
# -det[0] @ bell <= -1
# lhs_ineq = np.r_[lhs_ineq, [-1.0 * dets[0]]]
# rhs_ineq = np.r_[rhs_ineq, [-1.0]]

# 1 @ b >= 0 equiv to -1.0 @ b <= 0
#lhs_ineq = np.r_[lhs_ineq, -1.0 * np.eye(dets.shape[1])]
#rhs_ineq = np.r_[rhs_ineq, np.zeros(dets.shape[1])]

# (det + ray) @ bell <= 1

for comb in product([0, 1], repeat=lins.shape[0]):
    comb = np.array(comb)
    r = lins * comb[:, np.newaxis]
    r = np.sum(r, axis=0)
    for d in dets:
        continue
        lhs_ineq = np.r_[lhs_ineq, [r]]
        lhs_ineq = np.r_[lhs_ineq, [-r]]
        rhs_ineq = np.r_[rhs_ineq, [0.0]]
        rhs_ineq = np.r_[rhs_ineq, [0.0]]


# get the h representation of the polyhedra
hrepr = polyhedra_h_representation(lhs_ineq, rhs_ineq, file='input_h.ine')

# run the h representation of the polyhedra
vertices, rays = run_lrs_h_repr('input_h.ine')

# print the numbers
print('H polytope run')
print('Number of vertices: {}'.format(vertices.shape[0]))
print('number of rays: {}'.format(rays.shape[0]))
