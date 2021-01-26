"""
In this script I am trying to identify all ridges of a facet, that was calculated by adjacency decomposition.
The data was extracted manually from PANDA and usually Fourier-Motzkin-Elimination is used to find the
ridges. However Fourier-Motzkin can not be parallized and generate a lot of redundant facets.
This is why we will use LRS to check if it is faster to run it.
"""
import numpy as np
from linearbell.lrs_helper import polytope_v_representation, run_lrs_v_repr, run_redund

raw_file = '../data/example_vertices_on_facet/raw4422.txt'
vertex_file = 'vertices.ext'

# Read vertices from raw data
raw_data = np.loadtxt(raw_file)
vertices = []
for row in raw_data:
    v = row[:-1] / row[-1]
    vertices.append(v)
vertices = np.array(vertices)
print(vertices)

# write the v-representation for lrs
polytope_v_representation(vertices, file=vertex_file)

# run redund

# run lrs (with timing)
run_lrs_v_repr(vertex_file, 'ridges.ine', nproc=6)

# read the resulting ridges