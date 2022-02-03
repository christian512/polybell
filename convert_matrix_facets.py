"""
Convert facets obtained from: https://www-users.york.ac.uk/~rc973/Bell
to the representation that we use within this repository
"""

import scipy.io
import numpy as np
from linearbell.utils import get_configs
from linearbell.panda_helper import write_known_inequalities
from linearbell.representations import get_configs_mat, transform_mat_to_vec

# import facets
matrix_ineqs = scipy.io.loadmat('Inequalities_Only_4422_Affine.mat')



# get configurations
ma,mb,na,nb = 4,4,2,2
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
write_known_inequalities(lhs,rhs,'4422.ine')