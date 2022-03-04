"""
Convert facets obtained from: https://www-users.york.ac.uk/~rc973/Bell
to the representation that we use within this repository
"""

import argparse
import scipy.io
import numpy as np
from polybell.utils import get_configs
from polybell.panda_helper import write_known_inequalities
from polybell.representations import get_configs_mat, transform_mat_to_vec

# parse the inputs
parser = argparse.ArgumentParser()
parser.add_argument(dest='ma', help="number of inputs for ALICE")
parser.add_argument(dest='mb', help='number of inputs for BOB')
parser.add_argument(dest='na', help='number of outputs for ALICE')
parser.add_argument(dest='nb', help='number of outputs for BOB')
parser.add_argument(dest='input_file', type=str,
                    help='Path to input file (relative to this script). Input needs to be given in the standard PANDA/PORTA format.')
parser.add_argument(dest='output_file', type=str, help='Output filename')

args = parser.parse_args()

# import facets
try:
    matrix_ineqs = scipy.io.loadmat(args.input_file)
except:
    print('Error loading input file -> did you give the correct path and is it a matlab file?')
    exit(1)

print('Loaded matrices from Matlab file.')

# get configurations
ma, mb, na, nb = int(args.ma), int(args.mb), int(args.na), int(args.nb)
configs = get_configs(range(ma), range(mb), range(na), range(nb))
configs_mat = get_configs_mat(range(ma), range(mb), range(na), range(nb))

# transform each matrix to a vectorial representation
vecs = []
for mat in matrix_ineqs.values():
    vec = transform_mat_to_vec(configs_mat, configs, mat)
    # apply transformations that it is in the sme format as PANDA output
    vec = -1.0 * vec
    # add vector to all vecs
    vecs.append(vec)

# the vectors form the left hand side of the inequalities
lhs = np.array(vecs)
# rhs is just the -1
rhs = -1.0 * np.ones(lhs.shape[0])

# write out known inequalities
write_known_inequalities(lhs,rhs,args.output_file)
print('Wrote Inequalities in vetor form to: {}'.format(args.output_file))