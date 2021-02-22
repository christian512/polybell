""" Generates the 2222-case face lattice """

from linearbell.utils import get_deterministic_behaviors
import numpy as np
from face import Polytope
import matplotlib.pyplot as plt
import networkx as nx

G = nx.Graph()

inputs = range(2)
outputs = range(2)

dets = get_deterministic_behaviors(inputs, inputs, outputs)
relabels = np.loadtxt('../data/relabels/{}{}{}{}.gz'.format(2, 2, 2, 2)).astype(int)

faces_per_level = {}

def recursive_polytope_finder(p, parent_poly=None, level=0, all_polys={}):
    """ Finds all subpolytopes for face lattice structure """
    # G.add_node(p.id)
    num_dets = len(p.deterministics)
    print(p.id)
    # add polytope to all polytopes dict
    if num_dets not in all_polys.keys():
        all_polys[num_dets] = [p]

    else:
        # Check if polytope was already calculated
        for other_p in all_polys[num_dets]:
            if other_p == p:
                G.add_edge(parent_poly.id, other_p.id)
                return all_polys
        all_polys[num_dets].append(p)
    if level in faces_per_level.keys():
        faces_per_level[level] += 1
    else:
        faces_per_level[level] = 1

    # add face to graph
    G.add_node(p.id, pos=(faces_per_level[level], -level))
    if parent_poly:
        G.add_edge(parent_poly.id, p.id)
    # check if there will be subpolytopes
    if len(p.deterministics) <= 1:
        return all_polys

    # get the subpolytopes
    sub_polys = p.get_classes()
    for f in sub_polys:
        all_polys = recursive_polytope_finder(f, p, level + 1, all_polys)
    return all_polys


# create initial polytope
bell_polytope = Polytope(dets, relabels)
all_polytopes = recursive_polytope_finder(bell_polytope)
pos = nx.get_node_attributes(G, 'pos')
nx.draw(G, pos=pos, with_labels=False)
plt.savefig('hierarchy.png')
