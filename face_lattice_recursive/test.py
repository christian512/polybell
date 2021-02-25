from linearbell.utils import get_deterministic_behaviors, equiv_check_adjacency_panda
import numpy as np
from polytope import Polytope

inputs = range(3)
outputs = range(2)

# Generate Bell Polytope
dets = get_deterministic_behaviors(inputs, inputs, outputs)
relabels = np.loadtxt('../data/relabels/{}{}{}{}.gz'.format(3, 3, 2, 2)).astype(int)
bell_polytope = Polytope(dets, relabels)

# get first face of bell polytope
face = bell_polytope.find_single_face()
bell_polytope.add_child(face)

# Storage of all polytopes for each level
all_polys = {}
all_polys[0] = [bell_polytope]
all_polys[1] = [face]


def add_polytope_at_level(poly, k):
    """ Adds a polytope to the level k"""
    print('level: ',k )
    # +++ PHASE 1: Equivalence checking +++
    if k not in all_polys.keys():
        all_polys[k] = []
    # check if polytope is equivalent to any child of the parent
    if not poly.parent:
        print('Polytope has no parent')
        return False
    equiv = False
    for child in poly.parent.children:
        break
        if equiv_check_adjacency_panda(poly.creating_face, child.creating_face, poly.parent.relabellings,
                                       poly.parent.deterministics):
            equiv = True
            break
    if equiv:
        print('Polytope is equivalent to a child')
        return False
    # check if the polytope is equivalent to any polytope on that level
    for p in all_polys[k]:
        if equiv_check_adjacency_panda(poly.creating_face, p.creating_face, bell_polytope.relabellings, bell_polytope.deterministics):
            equiv = True
            break
    if equiv:
        print('Polytope is equivalent to another poly on that level')
        return False
    # add the polytope
    print('adding polytope at level:', k)
    all_polys[k].append(poly)
    # +++ PHASE 2 : Rotate above faces around this new polytope
    if k > 1:

        for p in all_polys[k-1]:
            # check if poly is a valid face of p
            if not p.add_child(poly):
                continue

            # rotate p around polytope
            new_p = p.rotate_polytope(poly)
            if np.all(new_p.creating_face == 0):
                continue
            # recursive call
            print('phase 2 -> try adding: ', new_p.creating_face)
            add_polytope_at_level(new_p, k-1)
    # +++ PHASE 3: Rotate poly around all lower dimensional faces

    if poly.deterministics.shape[0] > 1:
        pass
        if k + 1 not in all_polys.keys():
            all_polys[k + 1] = []
        for p in all_polys[k+1]:
            # try to add child
            if not poly.add_child(p):
                continue
            new_poly = poly.rotate_polytope(p)
            if np.all(new_poly.creating_face == 0):
                continue
            # recursive polytope
            print('phase 3 -> try adding')
            add_polytope_at_level(new_poly, k+1)
    # +++ PHASE 4: CALL THE HEURISTIC TO FIND A NEW LOWER DIMENSIONAL POLYTOPE +++
    if poly.deterministics.shape[0] > 1:
        new_poly = poly.find_single_face()
        add_polytope_at_level(new_poly, k+1)
    return True


ridge = face.find_single_face()
add_polytope_at_level(ridge, 2)