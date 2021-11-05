""" Calculates the vertices on the facet """

from linearbell.utils import get_deterministic_behaviors, get_configs
from linearbell.representations import get_configs_mat, transform_mat_to_vec

import numpy as np
from scipy.io import loadmat

configs = get_configs(range(4), range(4), range(2), range(2))
configs_mat = get_configs_mat(range(4), range(4), range(2), range(2))

vertices = get_deterministic_behaviors(range(4), range(4), range(2))
facets_mat = loadmat("../data/facets_known/Inequalities_Only_4422_Raw.mat")

facets = []
for k in facets_mat.keys():
    mat = facets_mat[k]
    f = transform_mat_to_vec(configs_mat, configs, mat)
    facets.append(f)

# as minimal value in this representation is one we check for equality with one
facet_vertices = []
for f in facets:
    curr_vertices = []
    for i in range(len(vertices)):
        if vertices[i] @ f == 1.0:
            curr_vertices.append(i+1)
    facet_vertices.append(curr_vertices)
    print('number vertices:  ', len(curr_vertices))

f = open("4422facet_vertices.txt", "w+")
f.write(str(facet_vertices))
f.close()











