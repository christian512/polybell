"""
This file defines functions for adjacency decomposition
"""
import numpy as np
from fractions import Fraction
from typing import Optional

def distance(vertex: np.ndarray, facet: np.ndarray) -> np.ndarray:
    assert vertex.shape[0] == facet.shape[0]
    return -1.0 * (vertex @ facet)

def furthest_vertex(vertices: np.ndarray, facet: np.ndarray) -> np.ndarray:
    """ Returns furthest vertex """
    assert vertices.shape[1] == facet.shape[0]
    best_vertex = vertices[0]
    d = distance(best_vertex, facet)
    for v in vertices[1:]:
        if distance(v, facet) > d:
            d = distance(v, facet)
            best_vertex = v
    return best_vertex

def nearest_vertex(vertices: np.ndarray, facet: np.ndarray) -> np.ndarray:
    """ Returns nearest vertex """
    assert vertices.shape[1] == facet.shape[0]
    best_vertex = vertices[0]
    d = distance(best_vertex, facet)
    for v in vertices[1:]:
        if distance(v, facet) < d:
            d = distance(v, facet)
            best_vertex = v
    return best_vertex


def rotate(vertices: np.ndarray, vertex: np.ndarray, facet: np.ndarray, ridge: np.ndarray) -> np.ndarray:
    """ Rotates a facet around a ridge """
    d_f = distance(vertex, facet)
    d_r = distance(vertex, ridge)
    first_round = True
    while d_r != 0 or first_round:
        first_round = False
        assert d_f % 1 == 0
        assert d_r % 1 == 0
        d_f = int(d_f)
        d_r = int(d_r)
        gcd_ds = np.gcd(d_f, d_r)
        if gcd_ds > 1:
            d_f /= gcd_ds
            d_r /= gcd_ds
        ridge = d_f * ridge - d_r * facet
        gcd_value = np.abs(np.gcd.reduce(ridge.astype(int)))
        if gcd_value > 1:
            ridge /= gcd_value
        vertex = nearest_vertex(vertices, ridge)
        d_f = distance(vertex, facet)
        d_r = distance(vertex, ridge)
    return ridge


