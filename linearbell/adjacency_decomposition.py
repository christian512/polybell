"""
This file defines functions for adjacency decomposition
"""
import numpy as np
from fractions import Fraction

def distance(vertex, facet):
    assert vertex.shape[0] == facet.shape[0]
    return -1.0 * (vertex @ facet)

def furthest_vertex(vertices, facet):
    """ Returns furthest vertex """
    assert vertices.shape[1] == facet.shape[0]
    best_vertex = vertices[0]
    d = distance(best_vertex, facet)
    for v in vertices[1:]:
        if distance(v, facet) > d:
            d = distance(v, facet)
            best_vertex = v
    return best_vertex

def nearest_vertex(vertices, facet):
    """ Returns nearest vertex """
    assert vertices.shape[1] == facet.shape[0]
    best_vertex = vertices[0]
    d = distance(best_vertex, facet)
    for v in vertices[1:]:
        if distance(v, facet) < d:
            d = distance(v, facet)
            best_vertex = v
    return best_vertex


def rotate(vertices, vertex, facet, ridge):
    # TODO: Try to Implement rotation not as in PANDA, but as you researched it.
    """ Rotates a facet around a ridge """
    counter = 0
    d_r = 1
    while not d_r == 0:
        d_f = distance(vertex, facet)
        d_r = distance(vertex, ridge)
        assert d_f % 1 == 0
        assert d_r % 1 == 0
        d_f = int(d_f)
        d_r = int(d_r)
        if counter > 1000:
            print('Ran over 1000 iterations in rotation')
            print('input facet:', facet)
            return False
        gcd_ds = np.gcd(d_f, d_r)
        if gcd_ds > 1:
            d_f /= gcd_ds
            d_r /= gcd_ds
        ridge = d_f * ridge - d_r * facet
        for val in ridge:
            assert val % 1.0 == 0, val
        gcd_value = np.abs(np.gcd.reduce(ridge.astype(int)))
        if gcd_value > 1:
            ridge /= gcd_value
        vertex = nearest_vertex(vertices, ridge)
        #print('d_f: ', d_f)
        #print('d_r: ', d_r)
        counter += 1

    return ridge


