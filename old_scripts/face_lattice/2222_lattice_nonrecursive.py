""" Generates the 2222-case face lattice """

from linearbell.utils import get_deterministic_behaviors, equiv_check_adjacency_testing
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

all_polys = {}


def recursive_polytope_finder(polys, level=0):
    """ Finds all subpolytopes for face lattice structure """
    print('level: ',level)
    if len(polys) == 0:
        return True
    # set all polytopes of this level
    all_polys[level] = polys
    # now for every polytope, calculate the faces
    new_polys = []
    for p in all_polys[level]:
        if len(p.deterministics) == 1:
            continue
        faces = p.get_faces()
        for f in faces:
            new_polys.append(f)
    # set first classes representative
    if len(new_polys) == 0:
        return True
    p = new_polys[0]
    new_polys_classes = [p]
    G.add_node(p.id, pos=(len(new_polys_classes) - 1, -level), ndets=len(p.deterministics), nrel=len(p.poss_relabellings),
               dims=p.dims, dets_indices=str(p.indices_deterministics))
    # check the equivalence of the others to this, draw edges
    for p in new_polys:
        equiv = False
        for c in new_polys_classes:
            tmp_dets = p.initial_polytope.deterministics
            tmp_relabels = p.initial_polytope.poss_relabellings
            if equiv_check_adjacency_testing(c.creating_face[:-1], p.creating_face[:-1], relabels=tmp_relabels,
                                             dets=tmp_dets):
                equiv = True
                # draw an edge from parent of p to c
                G.add_edge(p.parent.id, c.id)
                break
        if not equiv:
            # add to class
            new_polys_classes.append(p)
            # add node
            G.add_node(p.id, pos=(len(new_polys_classes)-1, -level), ndets=len(p.deterministics),
                       nrel=len(p.poss_relabellings),
                       dims=p.dims, dets_indices=str(p.indices_deterministics))
            # add edge
            G.add_edge(p.parent.id, p.id)
    return recursive_polytope_finder(new_polys_classes, level + 1)



# create initial polytope
bell_polytope = Polytope(dets, relabels)
p = bell_polytope
# add node for original polytpe
G.add_node(p.id, pos=(0, 0), ndets=len(p.deterministics),
           nrel=len(p.poss_relabellings),
           dims=p.dims, dets_indices=str(p.indices_deterministics))
all_polytopes = recursive_polytope_finder([bell_polytope], level=1)
pos = nx.get_node_attributes(G, 'pos')
network_graph = from_networkx(G, pos)
# Set node size and color
network_graph.node_renderer.glyph = Circle(size=15, fill_color='skyblue')

# Set edge opacity and width
network_graph.edge_renderer.glyph = MultiLine(line_alpha=0.5, line_width=1)
HOVER_TOOLTIPS = [("Number Deterministics", "@ndets"), ("Indices Deterministics", "@dets_indices"),
                  ("Number relabels", "@nrel"), ("Dimensions", "@dims")]
plot = figure(tooltips=HOVER_TOOLTIPS, x_range=Range1d(0, 20), y_range=Range1d(-8, 2),
              title='Face-Classes-Lattice for 2222 case')
plot.renderers.append(network_graph)
show(plot)
# save(plot, 'face_classes_lattice_2222.html')
