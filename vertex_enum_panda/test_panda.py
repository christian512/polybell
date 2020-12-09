import argparse
from linearbell.utils import get_deterministic_behaviors, get_configs, get_allowed_relabellings
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



# setup the constraints: dets @ bell <= 1
lhs_ineq = np.copy(dets)
rhs_ineq = np.ones(dets.shape[0])

# write the file
hrepr = write_panda_input_inequalities(lhs_ineq, rhs_ineq, file='input.ine')

# run the file
run_panda('input.ine')

