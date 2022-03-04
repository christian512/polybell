import argparse
from polybell.utils import get_deterministic_behaviors, get_allowed_relabellings, get_relabels_dets
import numpy as np

# Get arguments
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

# get allowed relabellings
allowed_relabellings = get_allowed_relabellings(inputs_a, inputs_b, outputs, outputs)
# store the relabels
file = '../data/relabels/{}{}{}{}.gz'.format(ma, mb, n, n)
np.savetxt(file, allowed_relabellings)
