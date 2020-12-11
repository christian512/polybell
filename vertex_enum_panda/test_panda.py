import argparse
from linearbell.utils import get_deterministic_behaviors, get_configs, general_pr_box_extended, get_allowed_relabellings
from linearbell.panda_helper import write_panda_input_inequalities, run_panda
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
outputs_wo_failure = range(n-1)

# setup output file
outfile = '../data/vertex_enum/{}{}{}{}.ext'.format(ma, mb, n, n)

# get deterministic points
dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs)

relabels = np.loadtxt('../data/relabels/{}{}{}{}.gz'.format(ma,mb,n,n)).astype(int)

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

for i in range(combinations.shape[0]):
    c = combinations[i]

    # setup the constraints: dets @ bell <= 1
    lhs = dets[c > 1e-4]
    lhs = np.r_[lhs, [pr_box]]
    rhs = np.ones(lhs.shape[0])
    lins = list(np.arange(lhs.shape[0]))
    lhs = np.r_[lhs, dets[c < 1e-4]]
    rhs = np.r_[rhs, np.ones(len(c[c < 1e-4]))]
# write the file
hrepr = write_panda_input_inequalities(lhs, rhs, idx_equalities=lins, symmetries=relabels[:2], file='input.ine')

# run the file
run_panda('input.ine')

