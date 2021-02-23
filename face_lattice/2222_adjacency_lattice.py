""" Generates the 2222-case face lattice """

from linearbell.utils import get_deterministic_behaviors, equiv_check_adjacency_testing
from linearbell.adjacency_decomposition import rotate, furthest_vertex
import numpy as np
from face import Polytope
import matplotlib.pyplot as plt
from bokeh.io import show, save
from bokeh.models import Circle, MultiLine, Range1d
from bokeh.plotting import figure, from_networkx
import networkx as nx

G = nx.Graph()

inputs = range(2)
outputs = range(2)

dets = get_deterministic_behaviors(inputs, inputs, outputs)
relabels = np.loadtxt('../data/relabels/{}{}{}{}.gz'.format(2, 2, 2, 2)).astype(int)

faces_per_level = {}
all_polys = {}


def recursive_polytope_finder(p, parent_poly=None, level=0):
    """ Finds all subpolytopes for face lattice structure """
    # G.add_node(p.id)
    print(p.id)
    # add polytope to all polytopes dict
    if level not in all_polys.keys():
        all_polys[level] = [p]

    else:
        # Check if polytope was already calculated
        for other_p in all_polys[level]:
            if other_p == p:
                G.add_edge(parent_poly.id, other_p.id)
                return True
        all_polys[level].append(p)
    if level in faces_per_level.keys():
        faces_per_level[level] += 1
    else:
        faces_per_level[level] = 1

    # add face to graph
    G.add_node(p.id, pos=(faces_per_level[level], -level), ndets=len(p.deterministics), nrel=len(p.poss_relabellings))
    if parent_poly:
        G.add_edge(parent_poly.id, p.id)
    # check if there will be subpolytopes
    if len(p.deterministics) <= 1:
        return True

    # get the subpolytopes
    sub_polys = p.get_classes()
    for f in sub_polys:
        recursive_polytope_finder(f, p, level + 1)
    return True


# create initial polytope
bell_polytope = Polytope(dets, relabels)
recursive_polytope_finder(bell_polytope)

# Check which faces can be pivoted into another face and display in the Graph.
G.remove_edges_from(list(G.edges()))
# iterate through the levels
for level in all_polys.keys():
    if level + 2 not in all_polys.keys():
        continue
    poly = all_polys[level][-1]
    # Choose the first polytope of this level
    face = all_polys[level + 1][-1]
    # take it's first child -> in this case is the first polytope of the level below
    ridge = all_polys[level + 2][-1]
    # rotate the face around the ridge
    vertex = furthest_vertex(poly.deterministics, face.face[:-1])
    new_face = rotate(poly.deterministics, vertex, face.face[:-1], ridge.face[:-1])
    # check if the new face is equivalent to any on the level of the original face.
    equiv = False
    for test_face in all_polys[level + 1]:
        equiv = equiv_check_adjacency_testing(test_face.face[:-1], new_face, test_face.poss_relabellings,
                                              test_face.deterministics)
        if equiv:
            G.add_edge(face.id, ridge.id)
            G.add_edge(ridge.id, test_face.id)
            break
    if not equiv:
        print('Not equivalent <- this should not happen')

pos = nx.get_node_attributes(G, 'pos')
network_graph = from_networkx(G, pos)
# Set node size and color
network_graph.node_renderer.glyph = Circle(size=15, fill_color='skyblue')

# Set edge opacity and width
network_graph.edge_renderer.glyph = MultiLine(line_alpha=0.5, line_width=1)
HOVER_TOOLTIPS = [("Number Deterministics", "@ndets"), ("Number relabels", "@nrel")]
plot = figure(tooltips=HOVER_TOOLTIPS, x_range=Range1d(0, 20), y_range=Range1d(-8, 2),
              title='Adjacency-Classes-Lattice for 2222 case')
plot.renderers.append(network_graph)
show(plot)
save(plot, 'adjacency_classes_lattice_2222.html')
