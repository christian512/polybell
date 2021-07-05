""" Generates the 2222-case face lattice """

from linearbell.utils import get_deterministic_behaviors, equiv_check_adjacency_panda
from linearbell.adjacency_decomposition import rotate
import numpy as np
from face import Polytope, polytope_from_face
import matplotlib.pyplot as plt
from bokeh.io import show, save
from bokeh.models import Circle, MultiLine, Range1d
from bokeh.plotting import figure, from_networkx
import networkx as nx

G = nx.Graph()

inputs = range(3)
outputs = range(2)

dets = get_deterministic_behaviors(inputs, inputs, outputs)
relabels = np.loadtxt('../data/relabels/{}{}{}{}.gz'.format(3, 3, 2, 2)).astype(int)
relabels_dets = np.loadtxt('../data/relabels_dets/{}{}{}{}.gz'.format(3, 3, 2, 2)).astype(int)

all_polys = {}


def polytope_finder(polys, level=0):
    """ Finds all subpolytopes for face lattice structure """
    print('level: ', level)
    # nothing to do if no polys given
    if len(polys) == 0:
        return True
    # set all polytopes of this level
    all_polys[level] = polys
    # now for every polytope, calculate the faces
    new_polys = []
    for p in all_polys[level]:
        if len(p.deterministics) == 1:
            continue
        # TODO: Here you can change if classes or faces should be taken
        faces = p.get_classes()
        for f in faces:
            new_polys.append(f)
    if len(new_polys) == 0:
        return True
    # check the equivalence on the level
    new_polys_classes = []
    for p in new_polys:
        equiv = False
        for c in new_polys_classes:
            tmp_dets = p.initial_polytope.deterministics
            tmp_relabels = p.initial_polytope.poss_relabellings
            if equiv_check_adjacency_panda(c.creating_face, p.creating_face, relabels=tmp_relabels,
                                           dets=tmp_dets):
                equiv = True
                # draw an edge from parent of p to c
                G.add_edge(p.parent.id, c.id)
                break
        if not equiv:
            # add to class
            new_polys_classes.append(p)
            # add node
            G.add_node(p.id, pos=(len(new_polys_classes) - 1, -level), ndets=len(p.deterministics),
                       nrel=len(p.poss_relabellings),
                       dims=p.dims, dets_indices=str(p.indices_deterministics), face=str(p.creating_face))
            # add edge
            G.add_edge(p.parent.id, p.id)

    return polytope_finder(new_polys_classes, level + 1)


# create initial polytope
bell_polytope = Polytope(dets, relabels)
p = bell_polytope
# add node for original polytpe
G.add_node(p.id, pos=(0, 0), ndets=len(p.deterministics),
           nrel=len(p.poss_relabellings),
           dims=p.dims, dets_indices=str(p.indices_deterministics))
polytope_finder([bell_polytope], level=1)
print('Creating Adjacency Graph')

# Create a copy with removed edges
G_adj = G.copy()
G_adj.remove_edges_from(list(G_adj.edges()))
# Go through the levels to determine adjacencies
for level in all_polys.keys():
    if not level + 2 in all_polys.keys():
        break
    # get polytope, faces and ridges -> by edges in the original tree
    for poly in all_polys[level]:
        for f1 in all_polys[level + 1]:
            for f2 in all_polys[level + 1]:
                if f1.id == f2.id:
                    continue

                for r in all_polys[level + 2]:
                    for g1 in relabels_dets:
                        # build the union of f1 and g(f2)
                        intersect = np.intersect1d([f1.indices_deterministics], g1[f2.indices_deterministics])
                        for g2 in relabels_dets:
                            if np.all(intersect == g2[r.indices_deterministics]):
                                G_adj.add_edge(f1.id, r.id)
                                G_adj.add_edge(r.id, f2.id)
                                break


pos = nx.get_node_attributes(G_adj, 'pos')
network_graph = from_networkx(G_adj, pos)
# Set node size and color
network_graph.node_renderer.glyph = Circle(size=15, fill_color='skyblue')

# Set edge opacity and width
network_graph.edge_renderer.glyph = MultiLine(line_alpha=0.5, line_width=1)
HOVER_TOOLTIPS = [("Number Deterministics", "@ndets"), ("Indices Deterministics", "@dets_indices"),
                  ("Number relabels", "@nrel"), ("Dimensions", "@dims"), ("Face", "@face")]
plot = figure(tooltips=HOVER_TOOLTIPS, x_range=Range1d(0, 20), y_range=Range1d(-8, 2),
              title='Face-Classes-Lattice for 3322 case')
plot.renderers.append(network_graph)
show(plot)
save(plot, '3322_class_lattice_adjacency.html')
