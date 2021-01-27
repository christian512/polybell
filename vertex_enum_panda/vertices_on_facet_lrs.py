"""
In this script I am trying to identify all ridges of a facet, that was calculated by adjacency decomposition.
The data was extracted manually from PANDA and usually Fourier-Motzkin-Elimination is used to find the
ridges. However Fourier-Motzkin can not be parallized and generate a lot of redundant facets.
This is why we will use LRS to check if it is faster to run it.
"""
import numpy as np
from linearbell.lrs_helper import polytope_v_representation, run_lrs_v_repr, run_redund, polyhedra_h_representation, run_lrs_h_repr

raw_file = '../data/example_vertices_on_facet/raw3322.txt'
vertex_file = 'vertices.ext'
ineq_file = 'ineq.ine'

# Read vertices from raw data
raw_data = np.loadtxt(raw_file)
vertices = []
lhs = []
rhs = []
for row in raw_data:
    v = row[:-1] / row[-1]
    lhs.append(row[:-1])
    rhs.append(row[-1])
    vertices.append(v)
vertices = np.array(vertices)
lhs = np.array(lhs)
rhs = np.array(rhs)


# write the v-representation for lrs
polytope_v_representation(vertices, file=vertex_file)

# write the h-representation
polyhedra_h_representation(lhs, rhs, file=ineq_file)

# run redund

# run lrs (with timing)
run_lrs_h_repr(ineq_file, output_file='ridges.ext')

# read the resulting ridges