"""
This script tests panda. It will find Bell Inequalities for finite efficiency PR-Boxes.
"""

import argparse
from linearbell.utils import get_deterministic_behaviors, get_configs, general_pr_box_extended, check_equiv_bell_vertex_enum_non_rescale
from linearbell.panda_helper import write_panda_input_inequalities, run_panda, read_vertices_rays
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument(dest='ma', help="number of inputs for ALICE")
parser.add_argument(dest='mb', help='number of inputs for BOB')
parser.add_argument(dest='n', help='number of outputs')
parser.add_argument(dest='threads', help='number of threads')
args = parser.parse_args()
ma = int(args.ma)
mb = int(args.mb)
n = int(args.n)
threads = int(args.threads)

# set inputs / outputs
inputs_a = range(ma)
inputs_b = range(mb)
outputs = range(n)
outputs_wo_failure = range(n - 1)

# setup output file
outfile = '../data/vertex_enum/{}{}{}{}.ext'.format(ma, mb, n, n)

# get deterministic points
dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs)
dets_unshifted = np.copy(dets)

relabels = np.loadtxt('../data/relabels/{}{}{}{}.gz'.format(ma, mb, n, n)).astype(int)

# setup pr box
configs = get_configs(inputs_a, inputs_b, outputs, outputs)
eta = 4 / (4 + len(inputs_a))
pr_box = [general_pr_box_extended(a, b, x, y, eta, outputs_wo_failure) for (a, b, x, y) in configs]
pr_box = np.array(pr_box)

# Define new origin, such that the 0 origin is inside of the polytope (for dual representation)
p_origin = np.sum(dets, axis=0) / dets.shape[0]

# shift the origin of the deterministic points
dets = dets - p_origin
pr_box = pr_box - p_origin

# get configurations
configs = get_configs(inputs_a, inputs_b, outputs, outputs)

# load the possible combinations to build the PR box from file
file = '../data/pr_box_contrib_det/{}{}{}{}.gz'.format(ma, mb, n, n)
combinations = np.loadtxt(file)
if len(combinations.shape) == 1:
    # reshape to 2D array
    combinations = combinations.reshape(1, combinations.shape[0])

# select one combination
c = combinations[0]

# setup the constraints: dets @ bell = 1 for dets contributing
lhs = dets[c > 1e-4]
lhs = np.r_[lhs, [pr_box]]
rhs = np.ones(lhs.shape[0])
lins = list(np.arange(lhs.shape[0]))
# setup constraints: dets @ bell <= 1 for dets not contributing
lhs = np.r_[lhs, dets[c <= 1e-4]]
rhs = np.r_[rhs, np.ones(len(c[c <= 1e-4]))]

# write input file for panda
hrepr = write_panda_input_inequalities(lhs, rhs, idx_equalities=lins, symmetries=relabels,dets=dets_unshifted, file='input.ine')

# run the file
run_panda('input.ine', outfile='out.ext', threads=threads)

# read vertices and rays
vertices, rays = read_vertices_rays('out.ext')

# clean the vertices for classes
# check how many classes
classes = [vertices[0]]
for b in vertices:
    equiv = False
    for c in classes:
        if check_equiv_bell_vertex_enum_non_rescale(c, b, relabels, dets_unshifted):
            equiv = True
            continue
    if not equiv:
        classes.append(b)
classes = np.array(classes)
print('Found {} classes.'.format(classes.shape[0]))

