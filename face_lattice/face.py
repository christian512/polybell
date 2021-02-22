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

    def __init__(self, deterministics, poss_relabellings):
        """ initialize the polytope """
        self.id = ''.join(random.choices(string.digits + string.ascii_letters, k=30))
        self.deterministics = deterministics
        self.poss_relabellings = poss_relabellings
        self.faces = []
        self.classes = []

    def __str__(self):
        out = 'polytope with {} vertices and {} allowed relabellings'.format(len(self.deterministics), len(self.poss_relabellings))
        if len(self.faces) > 0:
            out += '/ {} faces and {} classes'.format(len(self.faces), len(self.classes))
        return out

    def get_faces(self):
        """ Calculates all faces of a polytope and returns them as polytope objects """
        # run Panda
        write_known_vertices(self.deterministics, file='knownvertices.ext')
        cmd = 'panda_org knownvertices.ext -t 1 --method=dd > out.ine'
        out = subprocess.run(cmd, shell=True)
        self.faces = read_inequalities('out.ine')
        # calculate the classes
        self.__reduce_to_inequiv()
        # for every face generate a subpolytope object
        subpolytopes = []
        for f in self.faces:
            # find the deterministic points on the face
            sub_dets = np.array([v for v in self.deterministics if distance(v, f[:-1]) == -1.0 * f[-1]])
            # find allowed relabellings
            sub_poss_relabellings = []
            for r in self.poss_relabellings:
                allowed = False
                for sd in sub_dets:
                    if np.any(np.sum((sub_dets - sd[r])**2, axis=1) < 1e-6):
                        allowed = True
                        break
                if allowed:
                    sub_poss_relabellings.append(r)
            sub_poss_relabellings = np.array(sub_poss_relabellings)
            # create new polytope

            sub_p = Polytope(sub_dets, sub_poss_relabellings)
            subpolytopes.append(sub_p)
        return subpolytopes

    def get_classes(self):
        """ Calculates all faces of a polytope and returns them as polytope objects """
        # run Panda
        write_known_vertices(self.deterministics, file='knownvertices.ext')
        cmd = 'panda_org knownvertices.ext -t 1 --method=dd > out.ine'
        out = subprocess.run(cmd, shell=True)
        self.faces = read_inequalities('out.ine')
        # calculate the classes
        self.__reduce_to_inequiv()
        # for every face generate a subpolytope object
        subpolytopes = []
        for f in self.classes:
            # find the deterministic points on the face
            sub_dets = np.array([v for v in self.deterministics if distance(v, f[:-1]) == -0.0 * f[-1]])
            # find allowed relabellings
            sub_poss_relabellings = []
            for r in self.poss_relabellings:
                allowed = False
                for sd in sub_dets:
                    if np.any(np.sum((sub_dets - sd[r])**2, axis=1) < 1e-6):
                        allowed = True
                        break
                if allowed:
                    sub_poss_relabellings.append(r)
            sub_poss_relabellings = np.array(sub_poss_relabellings)
            # create new polytope

            sub_p = Polytope(sub_dets, sub_poss_relabellings)
            subpolytopes.append(sub_p)
        return subpolytopes



    def __reduce_to_inequiv(self):
        """ Reduces the own faces to inequivalent classes """
        if len(self.faces) == 0:
            print('Running reduce to inequiv_classes with no faces given')

        ineq_bells = [self.faces[0]]
        for i, b in enumerate(self.faces):
            # print('equiv check progress : {} / {}'.format(i, bells.shape[0]))
            equiv = False
            for c in ineq_bells:
                if equiv_check_adjacency_testing(b[:-1], c[:-1], self.poss_relabellings, self.deterministics, tol=1e-6):
                    equiv = True
                    break
            if not equiv:
                ineq_bells.append(b)
        self.classes = ineq_bells

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if len(other.deterministics) == len(self.deterministics):
                return np.all(self.deterministics == other.deterministics)
        return False
