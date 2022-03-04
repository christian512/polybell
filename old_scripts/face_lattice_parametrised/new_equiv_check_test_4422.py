""" Checking if the new equivalence check works for 4422 """

from polybell.utils import get_deterministic_behaviors, get_configs, get_parametrisation_configs, parametrise_behavior
from polybell.representations import transform_mat_to_vec, get_configs_mat
import numpy as np
from polytope import ParamPolytope, parampolytope_from_inequality
from scipy.io import loadmat

inputs = range(4)
outputs = range(2)
configs = get_configs(inputs, inputs, outputs, outputs)
configs_mat = get_configs_mat(inputs, inputs, outputs, outputs)

# load inequalities found by TOM
facets_tmp = loadmat('../data/facets_known/Inequalities_Only_4422_Raw.mat')
facets_known = []
for k in facets_tmp.keys():
    mat = facets_tmp[k]
    f = transform_mat_to_vec(configs_mat, configs, mat)
    f = np.r_[f, -1.0]
    facets_known.append(f)
print('Generated known facets')

# Generate Bell Polytope
vertices = get_deterministic_behaviors(inputs, inputs, outputs)
permutations_vertices = np.loadtxt('../data/relabels_dets/{}{}{}{}.gz'.format(4, 4, 2, 2)).astype(int)
print('Loaded permutations of vertices')
bell_polytope = ParamPolytope(vertices, permutations_vertices)

# Creation of ParamPolytopes
faces = []
for f in facets_known:
    vertices_on_face = []
    indices_vertices_on_face = []
    for i in range(vertices.shape[0]):
        if vertices[i] @ f[:-1] == -1.0 * f[-1]:
            vertices_on_face.append(vertices[i])
            indices_vertices_on_face.append(i)
    # generate polytope
    p = ParamPolytope(np.array(vertices_on_face), np.array([]), creating_face=f,
                      indices_vertices=indices_vertices_on_face, parent=bell_polytope, initial_polytope=bell_polytope)
    assert not p.equiv_under_parent(faces)
    faces.append(p)
    print(p)

