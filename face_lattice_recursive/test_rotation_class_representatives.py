""" This tests how the rotation of faces around representative of a class works """

from linearbell.utils import get_deterministic_behaviors, equiv_check_adjacency_panda
import numpy as np
from polytope import Polytope

inputs = range(2)
outputs = range(2)

# Generate Bell Polytope
dets = get_deterministic_behaviors(inputs, inputs, outputs)
relabels = np.loadtxt('../data/relabels/{}{}{}{}.gz'.format(2, 2, 2, 2)).astype(int)
bell_polytope = Polytope(deterministics=dets, relabellings=relabels)

# get faces
faces = bell_polytope.get_all_faces()

# get all unequivalent ridges (both under face and bell polytope)
ridges = []
for f in faces:
    for tmp_ridge in f.get_all_classes():
        if not tmp_ridge.equiv_under_bell(ridges):
            ridges.append(tmp_ridge)

# rotate each face around each ridge class
for i in range(len(faces)):
    for j in range(len(ridges)):
        f = faces[i]
        r = ridges[j]
        new_r = f.get_valid_face_relabelling(r)
        if not new_r:
            print('Face {} , ridge {} : Ridge is not valid'.format(i, j))
            new_r = r
        if new_r != r:
            print('Face {} , ridge {} : Using other representation of ridge'.format(i, j))
        new_f = f.rotate_polytope(new_r)
        if new_f != f:
            print('Face {} , ridge {} : Lead to new face'.format(i, j))
        if new_f not in faces:
            print('Face {} around ridge {} leads to non valid face'.format(i, j))
