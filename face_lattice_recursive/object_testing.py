""" This performs a small test if everything works """

from linearbell.utils import get_deterministic_behaviors, equiv_check_adjacency_panda
import numpy as np
from polytope import Polytope

inputs = range(2)
outputs = range(2)

# Generate Bell Polytope
dets = get_deterministic_behaviors(inputs, inputs, outputs)
relabels = np.loadtxt('../data/relabels/{}{}{}{}.gz'.format(2, 2, 2, 2)).astype(int)

bell_polytope = Polytope(dets, relabels)
faces = bell_polytope.get_all_faces()
inequalties = np.array([f.creating_face for f in faces])
ridges = []
for f in faces:
    for r in f.get_all_faces():
        ridges.append(r)
        # rotate f around this ridge
        new_f = f.rotate_polytope(r)
        if not new_f in faces:
            print('rotated not in faces')
        # check if rotation resulted in the same
        if new_f == f:
            print('results in same')
        # check if ridge is in faces of the new face
        if r not in new_f.get_all_faces():
            print('ridge is not a ridge of the new face')

