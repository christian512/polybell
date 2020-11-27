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
args = parser.parse_args()
ma = int(args.ma)
mb = int(args.mb)
n = int(args.n)
nproc = int(args.nproc)

# set inputs / outputs
inputs_a = range(ma)
inputs_b = range(mb)
outputs = range(n)

# setup output file
outfile = '../data/vertex_enum/{}{}{}{}.ext'.format(ma, mb, n, n)

# get deterministic points
dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs)

# Define new origin, such that the 0 origin is inside of the polytope (for dual representation)
p_origin = np.sum(dets, axis=0) / dets.shape[0]

# shift the origin of the deterministic points
dets = dets - p_origin

# get configurations
configs = get_configs(inputs_a, inputs_b, outputs, outputs)

# Parametrise to get lower dimensional representations
configs_param = get_parametrisation_configs(inputs_a, inputs_b, outputs, outputs)

# Parametrise behaviors
dets_param = []
for d in dets:
    d_param = parametrise_behavior(d, configs, configs_param, inputs_a, inputs_b, outputs, outputs)
    dets_param.append(d_param)
dets_param = np.array(dets_param)
dets = np.copy(dets_param)

# setup the constraints: dets @ bell <= 1
lhs_ineq = np.copy(dets)
rhs_ineq = np.ones(dets.shape[0])

# get the h representation of the polyhedra
hrepr = polyhedra_h_representation(lhs_ineq, rhs_ineq, linearities=[0], file='input_h.ine')

# run the h representation of the polyhedra
vertices, rays = run_lrs_h_repr('input_h.ine', output_file=outfile, nproc=nproc)

# print the numbers
print('number of vertices: {}'.format(vertices.shape[0]))
print('number of rays: {}'.format(rays.shape[0]))
