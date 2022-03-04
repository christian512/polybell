"""
This file includes functions to represent bell inequalities in different ways.
For example we need the matrix form to compare to the previous results.
"""

from itertools import product
import numpy as np

def get_configs_mat(inputs_a, inputs_b, outputs_a, outputs_b):
    """ Returns the configurations in a matrix form """
    configs = []
    for x, a in product(inputs_a, outputs_a):
        c = []
        for y, b in product(inputs_b, outputs_b):
            c.append((a, b, x, y))
        configs.append(c)
    return configs


def transform_vec_to_mat(configs_vec, configs_mat, vec):
    """ Transforms a vector to matrix representation """
    mat = np.zeros((len(configs_mat), len(configs_mat[0])))
    for l, c in enumerate(configs_vec):
        for i in range(len(configs_mat)):
            try:
                j = configs_mat[i].index(c)
                # if not fails -> set the entry in matrix
                mat[i, j] = vec[l]
            except ValueError as e:
                pass
    return mat

def transform_mat_to_vec(configs_mat, configs_vec, mat):
    """ Transforms a matrix to vector representation """
    vec = np.zeros(len(configs_vec))
    for i, line in enumerate(configs_mat):
        for j, c in enumerate(configs_mat[i]):
            k = configs_vec.index(c)
            vec[k] = mat[i,j]
    return vec