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
    all_polys[level + 1] = []
    for p in all_polys[level]:
        if p.dims <= 0:
            continue
        # iterate through each face-class of this polytope
        for c in p.get_all_classes():
            # check if the class is equivalent under the bell polytope to already found polytopes
            if not c.equiv_under_bell(all_polys[level + 1]):
                all_polys[level + 1].append(c)
    recursive_classes_lattice(level + 1)


# build up the lattice
recursive_classes_lattice(0)

# TODO: Get network graph

# select faces and ridges to test
faces = all_polys[7]
ridges = all_polys[8]
print('Number of faces: ', len(faces))
print('number of ridges: ', len(ridges))

# iterate through each face
for i in range(len(faces)):
    f1 = faces[i]
    # indicator if a new face was found when rotating
    found_new = False
    # iterate through each ridge
    for j in range(len(ridges)):
        ridge = ridges[j]
        # identifier string
        prepend = 'face {}, ridge {}: '.format(i, j)
        # check if the ridge is valid and get the correct polytope representation
        new_ridge = f1.get_valid_face(ridge)
        if not new_ridge:
            print(prepend + 'Not a valid ridge')
            continue
        # set the ridge to be the valid ridge
        ridge = new_ridge
        # iterate through the relabellings of the face and rotate face around each relabelled version of the ridge
        for relabel in f1.relabellings:
            # relabel ridge
            new_ineq = ridge.creating_face[relabel]
            new_ineq = np.r_[new_ineq, ridge.creating_face[-1]]
            new_ridge = polytope_from_inequality(new_ineq, f1)
            new_ridge = f1.get_valid_face(new_ridge)
            assert new_ridge
            # rotate f1 around the new ridge
            f2 = f1.rotate_polytope(new_ridge)
            if np.all(f2.creating_face == 0):
                print('all zeros')
                continue
            if not f1.equiv_under_parent(f2):
                if not f1.equiv_under_bell(f2):
                    found_new = True
                    print(prepend + 'Found new face!')
                    # TODO: Add edges to the Graph
                    assert f2.equiv_under_bell(faces), 'Found completely new class'

# TODO: Draw the graph