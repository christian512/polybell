""" Iterates the class-lattice on the already existing lattice """

from polybell.utils import get_deterministic_behaviors, equiv_check_adjacency_panda
import numpy as np
from polytope import Polytope, polytope_from_inequality, create_nx_graph, create_bokeh_plot
from bokeh.io import show, save

inputs = range(2)
outputs = range(2)

# dict to store all polytopes
all_polys = {}

# Generate Bell Polytope
dets = get_deterministic_behaviors(inputs, inputs, outputs)
relabels = np.loadtxt('../data/relabels/{}{}{}{}.gz'.format(2, 2, 2, 2)).astype(int)
bell_polytope = Polytope(deterministics=dets, relabellings=relabels)
all_polys[0] = [bell_polytope]


def recursive_classes_lattice(level):
    """ Recursive construction of a classes face lattice """
    print('Level: ', level)
    if len(all_polys[level]) == 0:
        return True
    all_polys[level + 1] = []
    for p in all_polys[level]:
        if p.dims <= 0:
            continue
        # iterate through each face-class of this polytope
        for c in p.get_all_faces():
            # check if the class is equivalent under the bell polytope to already found polytopes
            if not c.equiv_under_bell(all_polys[level + 1]):
                all_polys[level + 1].append(c)
    recursive_classes_lattice(level + 1)


# build up the lattice
recursive_classes_lattice(0)

G = create_nx_graph(all_polys)

# dict of polytopes enumerated by the recursive strategy
all_polys_new = {}
all_polys_new[0] = [bell_polytope]

cut_off_level = 8

def add_polytope_to_level(poly, level):
    """ Add a polytope to a level """
    print(level)
    # check that level is in the storage
    if level not in all_polys_new.keys():
        all_polys_new[level] = []
    if level + 1 not in all_polys_new.keys():
        all_polys_new[level+1] = []
    # Temporary check due to invalid rotations
    o = poly.equiv_under_bell(all_polys[level])
    if not o:
        return False
    # STEP 1: Check if poly is equivalent to any polytope on that level -> then do nothing
    if poly.equiv_under_bell(all_polys_new[level]):
        return True
    else:
        all_polys_new[level].append(poly)
    # STEP 2: rotate above polytopes around that poly -> recursive call on the newly found polytopes
    if level > 1:
        for upper_poly in all_polys_new[level-1]:
            ridge = upper_poly.get_valid_face_relabelling(poly)
            if ridge:
                new_poly = upper_poly.rotate_polytope(ridge)
                add_polytope_to_level(new_poly, level - 1)

    # STEP 3: rotate poly around all polytopes below that level -> recursive call on the newly found polytopes
    if poly.dims == 0:
        return True
    for ridge in all_polys_new[level + 1]:
        new_ridge = poly.get_valid_face_relabelling(ridge)
        if new_ridge:
            new_poly = poly.rotate_polytope(new_ridge)
            add_polytope_to_level(new_poly, level)


    # STEP 4: Find a face of this and perform recursive call
    face = poly.get_single_face()
    add_polytope_to_level(face, level + 1)
    # add_polytope_to_level(all_polys[level+1][0], level + 1)

initial_face = bell_polytope.get_all_faces()[0]
add_polytope_to_level(initial_face, 1)

G = create_nx_graph(all_polys_new)
plot = create_bokeh_plot(G)
show(plot)