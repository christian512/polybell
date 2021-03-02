""" Checks if the relabelling around the lowest polytope (and all it's relabellings, results in the other class) """


from linearbell.utils import get_deterministic_behaviors, equiv_check_adjacency_panda
import numpy as np
from polytope import Polytope, polytope_from_inequality

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
    if len(all_polys[level]) == 0:
        return True
    all_polys[level+1] = []
    for p in all_polys[level]:
        if p.dims <= 0:
            continue
        # iterate through each face-class of this polytope
        for c in p.get_all_classes():
            # check if the class is equivalent under the bell polytope to already found polytopes
            if not c.equiv_under_bell(all_polys[level+1]):
                all_polys[level+1].append(c)
    recursive_classes_lattice(level+1)


recursive_classes_lattice(0)

# take the first class from level 6
f1 = all_polys[7][0]
ridge = all_polys[8][0]

# Since this is equivalent we want to see if we can rotate around a class representative that results in another class
for relabel in f1.relabellings:
    new_ineq = ridge.creating_face[relabel]
    new_ineq = np.r_[new_ineq, ridge.creating_face[-1]]
    # polytope from face
    new_ridge = polytope_from_inequality(new_ineq, f1)
    # rotate f1 around the new polytope
    f2 = f1.rotate_polytope(new_ridge)
    if np.all(f2.creating_face == 0):
        continue
    if not f1.equiv_under_parent(f2):
        if not f1.equiv_under_bell(f2):
            print('found new class under parent and bell')
            if not f2.equiv_under_bell(all_polys[6]):
                print('Found completly new class -> This case should not happen')
        else:
            print('equiv under bell')
    else:
        print('equiv under parent')



# plot = print_face_lattice(all_polys)
# from bokeh.io import show, save
# show(plot)

