from lrs_helper import polyhedra_h_representation, run_lrs_polytope
import numpy as np

# setup the problem I want to solve
lhs = np.eye(3)
rhs = np.ones(3) / 2
lhs = np.r_[lhs, -1.0 * np.eye(3)]
rhs = np.r_[rhs, 1.0 * np.ones(3) / 3]

# write out the string
polyhedra_h_representation(lhs, rhs, file='lrs_test.ine')

# start the calculation with lrs
vertices = run_lrs_polytope('lrs_test.ine', 'out.ext')
print(vertices)
