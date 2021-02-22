""" Generates the 2222-case face lattice """

from linearbell.utils import get_deterministic_behaviors
import numpy as np
from face import Polytope

inputs = range(2)
outputs = range(2)

dets = get_deterministic_behaviors(inputs, inputs, outputs)
relabels = np.loadtxt('../data/relabels/{}{}{}{}.gz'.format(2, 2, 2, 2)).astype(int)


def recursive_polytope_finder(p, level=0, all_polys={}):
    """ Finds all subpolytopes for face lattice structure """

    num_dets = len(p.deterministics)
    print(p.id)
    # add polytope to all polytopes dict
    if num_dets not in all_polys.keys():
        all_polys[num_dets] = [p]
    else:
        # TODO: Iterate through this
        for other_p in all_polys[num_dets]:
            if other_p == p:
                return all_polys
        all_polys[num_dets].append(p)
    # check if there will be subpolytopes
    if len(p.deterministics) <= 10:
        return all_polys

    # get the subpolytopes
    sub_polys = p.get_faces()
    for f in sub_polys:
        all_polys = recursive_polytope_finder(f, level + 1, all_polys)
    return all_polys


# create initial polytope
bell_polytope = Polytope(dets, relabels)
all_polytopes = recursive_polytope_finder(bell_polytope)
