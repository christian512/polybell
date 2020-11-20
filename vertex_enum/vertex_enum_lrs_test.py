import argparse
from linearbell.utils import get_deterministic_behaviors, get_configs, general_pr_box_extended, \
    get_parametrisation_configs, parametrise_behavior
import numpy as np
from lrs_helper import polyhedra_h_representation, run_lrs_polytope

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

# get the configurations
configs = get_configs(inputs_a, inputs_b, outputs, outputs)

# Parametrise to get lower dimensional representations
configs_param = get_parametrisation_configs(inputs_a, inputs_b, outputs, outputs)
# get parametrisations of deterministic points
dets_param = []
for d in dets:
    d_param = parametrise_behavior(d, configs, configs_param, inputs_a, inputs_b, outputs, outputs)
    dets_param.append(d_param)
dets_param = np.array(dets_param)
#dets = np.copy(dets_param)

# setup the constraints

# dets @ b <= 1
lhs_ineq = np.copy(dets)
rhs_ineq = np.ones(dets.shape[0])

# -1.0 * b <= 0 equiv to b >= 0
lhs_ineq = np.r_[lhs_ineq, -1.0 * np.eye(dets.shape[1])]
rhs_ineq = np.r_[rhs_ineq, np.zeros(dets.shape[1])]

# compute vertices of the polytope
print('Start computation of vertices')
polyhedra_h_representation(lhs_ineq, rhs_ineq, file='lrs.ine')
vertices = run_lrs_polytope('lrs.ine', 'out.ext')

#vertices = compute_polytope_vertices(lhs_ineq, rhs_ineq)

# store the vertices
vertices = np.array(vertices)
file = '../data/vertex_enum/{}{}{}{}.txt'.format(ma, mb, n, n)
np.savetxt(file, vertices)
print('done')
print('number of vertices: {}'.format(vertices.shape[0]))
print(vertices)