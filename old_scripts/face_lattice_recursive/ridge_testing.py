""" This performs tests on the generation of the new face by rotation and the shared ridge """

from polybell.utils import get_deterministic_behaviors, equiv_check_adjacency_panda
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
        assert r.dims == f.dims - 1
        ridges.append(r)
        # rotate f around this ridge
        new_f = f.rotate_polytope(r)
        # check if new face is in faces
        assert new_f in faces
        # check if new face was found
        assert new_f != f
        # check that all deterministics of the ridge are also deterministics of the new face
        assert np.all(np.isin(r.indices_deterministics, new_f.indices_deterministics))
        # check that the dimensions of the ridge is als
        assert r.dims == new_f.dims - 1
        # check that intersection of F and H is R
        intersection_indices = np.intersect1d(r.indices_deterministics, new_f.indices_deterministics)
        assert np.all(intersection_indices == r.indices_deterministics)

        # check if ridge is in faces of the new face
        new_ridges = new_f.get_all_faces()
        if r not in new_ridges:
            print('ridge is not a ridge of the new face')

