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
dets = np.copy(dets_param)

# setup the constraints
lhs_ineq = np.copy(dets)
lhs_ineq = np.r_[
    lhs_ineq, -1.0 * np.eye(dets.shape[1])]  # constraint that each entry of the bell inequality should be positive

rhs_ineq = np.ones(dets.shape[0])
rhs_ineq = np.r_[rhs_ineq, np.zeros(dets.shape[1])]  # positivity constraint

# compute vertices of the polytope
print('Start computation of vertices')
vertices = compute_polytope_vertices(lhs_ineq, rhs_ineq)

# store the vertices
vertices = np.array(vertices)
file = '../data/vertex_enum/{}{}{}{}.txt'.format(ma, mb, n, n)
np.savetxt(file, vertices)
print('done')
print('number of vertices: {}'.format(vertices.shape[0]))