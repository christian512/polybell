""" A class which handles a face """

from linearbell.panda_helper import write_known_vertices, read_inequalities
from linearbell.utils import equiv_check_adjacency_testing
from linearbell.adjacency_decomposition import distance
import numpy as np
import subprocess
import string
import random


class Polytope():
    """ Object Describing a polytope """

    def __init__(self, deterministics, poss_relabellings, creating_face=np.array([]), parent=None, initial_polytope=None):
        """ initialize the polytope """
        self.id = ''.join(random.choices(string.digits + string.ascii_letters, k=30))
        self.deterministics = deterministics
        # set parent polytope (self if toplevel)
        self.parent = parent
        if not self.parent:
            self.parent = self
        # set the initial polytope (self if toplevel)
        self.initial_polytope = initial_polytope
        if not self.initial_polytope:
            self.initial_polytope = self
        self.poss_relabellings = poss_relabellings
        self.faces = []
        self.classes = []
        # this is the face of the parent polytope that was used to create this polytope
        self.creating_face = creating_face

    def __str__(self):
        out = 'polytope with {} vertices and {} allowed relabellings'.format(len(self.deterministics),
                                                                             len(self.poss_relabellings))
        if len(self.faces) > 0:
            out += '/ {} faces and {} classes'.format(len(self.faces), len(self.classes))
        return out

    def get_classes(self):
        """ Calculates all faces of a polytope and returns them as polytope objects """
        # run Panda
        write_known_vertices(self.deterministics, file='knownvertices.ext')
        cmd = 'panda_org knownvertices.ext -t 1 --method=dd > out.ine'
        out = subprocess.run(cmd, shell=True)
        self.faces = read_inequalities('out.ine').astype(float)
        # rescale the faces
        for i in range(self.faces.shape[0]):
            if self.faces[i, -1] == -1:
                continue
            n = np.sum(self.deterministics[0])
            self.faces[i, :-1] = self.faces[i, :-1] + (self.faces[i, -1] + 1) / n
            self.faces[i, -1] = -1
        # calculate the classes
        self.__reduce_to_inequiv()
        # for every face generate a subpolytope object
        subpolytopes = []
        for f in self.classes:
            # find the deterministic points on the face
            sub_dets = np.array([v for v in self.deterministics if v @ f[:-1] == -1.0 * f[-1]])
            # find allowed relabellings
            sub_poss_relabellings = []
            for r in self.poss_relabellings:
                allowed = False
                for sd in sub_dets:
                    if np.any(np.sum((sub_dets - sd[r]) ** 2, axis=1) < 1e-6):
                        allowed = True
                        break
                if allowed:
                    sub_poss_relabellings.append(r)
            sub_poss_relabellings = np.array(sub_poss_relabellings)
            # create new polytope

            sub_p = Polytope(sub_dets, sub_poss_relabellings, creating_face=f,
                             parent=self, initial_polytope=self.initial_polytope)
            subpolytopes.append(sub_p)
        return subpolytopes

    def __reduce_to_inequiv(self):
        """ Reduces the own faces to inequivalent classes """
        if len(self.faces) == 0:
            print('Running reduce to inequiv_classes with no faces given')
            return False

        ineq_bells = [self.faces[0]]
        for i, b in enumerate(self.faces):
            equiv = False
            for c in ineq_bells:
                # TODO: HERE WE CAN EITHER CHOOSE ALL DETRMINISTICS OR ONLY THE ONE of the current face
                if equiv_check_adjacency_testing(b[:-1], c[:-1], self.poss_relabellings, self.deterministics,
                                                 tol=1e-6):
                    equiv = True
                    break
            if not equiv:
                ineq_bells.append(b)
        self.classes = ineq_bells

    def __eq__(self, other):
        """ Equality definition dependent on equivalence """
        if not isinstance(other, self.__class__):
            return False

        if len(other.creating_face) == 0 or len(self.creating_face) == 0:
            print('Comparison of polytopes that have no face')
            return False
        # Equivalence checking with both polytopes dets and relabellings
        equiv1 = False
        equiv2 = False
        # Perform the equivalence check with the parent deterministic points, as we check the polytopes themselve for equivalence.
        equiv1 = equiv_check_adjacency_testing(other.creating_face[:-1], self.creating_face[:-1],
                                               relabels=self.parent.poss_relabellings,
                                               dets=self.parent.deterministics)
        equiv2 = equiv_check_adjacency_testing(other.creating_face[:-1], self.creating_face[:-1],
                                               relabels=other.parent.poss_relabellings,
                                               dets=other.parent.deterministics)
        return equiv1 or equiv2
