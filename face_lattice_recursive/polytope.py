""" A class to describe a polytope """
import numpy as np
import subprocess
import string
import random
from linearbell.panda_helper import run_panda_vertices_on_facet, read_inequalities, write_known_vertices
from linearbell.adjacency_decomposition import furthest_vertex, rotate
from linearbell.utils import equiv_check_adjacency_panda


class Polytope():
    """ Describing a polytope """

    def __init__(self, deterministics, relabellings, creating_face=np.array([]), parent=None, initial_polytope=None):
        """ initialize polytope """
        self.id = ''.join(random.choices(string.digits + string.ascii_letters, k=30))
        self.deterministics = deterministics
        self.relabellings = relabellings
        # set the face that created this
        self.creating_face = creating_face
        # if no initial polytope is given
        self.initial_polytope = initial_polytope
        if not self.initial_polytope:
            self.initial_polytope = self
        # set parent
        self.parent = parent
        if not self.parent:
            self.parent = self
        # dimensions
        self.dims = np.linalg.matrix_rank(self.deterministics) - 1
        # empty faces
        self.__faces = []
        # empty classes
        self.__classes = []
        # indices of the deterministics
        self.indices_deterministics = []
        for det in self.deterministics:
            equal = det == self.initial_polytope.deterministics
            equal = np.all(equal, axis=1)
            idx = np.where(equal)[0][0]
            self.indices_deterministics.append(idx)
        self.indices_deterministics = np.array(self.indices_deterministics)

    def __str__(self):
        return "polytope {} with {} deterministics  ".format(self.id, self.deterministics.shape[0])

    def __eq__(self, other):
        if np.all(self.deterministics == other.deterministics):
            return True
        return False

    def get_all_faces(self):
        """ Finds all faces of the polytope """
        if len(self.__faces) != 0:
            return self.__faces
        # Run panda to get faces
        write_known_vertices(self.deterministics, file='knownvertices.ext')
        cmd = 'panda_org knownvertices.ext -t 1 --method=dd > out.ine'
        out = subprocess.run(cmd, shell=True)
        inequalities = read_inequalities('out.ine').astype(float)
        # generate polytopes from faces
        for ineq in inequalities:
            p = polytope_from_inequality(ineq, self)
            self.__faces.append(p)
        return self.__faces

    def get_all_classes(self):
        """ Finds all classes of faces of a polytope """
        if len(self.__classes) != 0:
            return self.__classes
        if len(self.__faces) == 0:
            self.get_all_faces()
        # Reduce the faces to inequivalent ones
        # TODO: Maybe a bit more advanced thing to choose the class representative?
        self.__classes = [self.__faces[0]]
        for i, f in enumerate(self.__faces):
            equiv = False
            for c in self.__classes:
                if equiv_check_adjacency_panda(f.creating_face, c.creating_face, self.relabellings,
                                               self.deterministics):
                    equiv = True
                    break
            if not equiv:
                self.__classes.append(f)
        return self.__classes

    def get_single_face(self):
        """ Finds a single face of the polytope, creates a new polytope and returns it """
        # run the modified version of panda
        write_known_vertices(self.deterministics, file='input.ext')
        run_panda_vertices_on_facet('input.ext', outfile='out.ine')
        ineq = read_inequalities('out.ine').astype(float)[0]
        return polytope_from_inequality(ineq, self)

    def rotate_polytope(self, ridge):
        """ Rotates the polytope around a ridge """
        if self.parent == self:
            print('Can not rotate as no parent polytope.')
            return False
        # set structure
        polytope = self.parent
        f = np.copy(self.creating_face[:-1])
        f0 = np.copy(-1.0 * self.creating_face[-1])
        g = np.copy(ridge.creating_face[:-1])
        g0 = np.copy(-1.0 * ridge.creating_face[-1])
        assert g.shape[0] == f.shape[0]

        # find the initial point
        assert polytope.deterministics.shape[1] == f.shape[0]
        products = polytope.deterministics @ f
        assert products.shape[0] == polytope.deterministics.shape[0]
        idx = np.where(products == np.amin(products))[0][0]
        vertex = polytope.deterministics[idx]
        assert vertex @ f < f0
        # start the iteration
        counter = 0
        while counter < 1000:
            counter += 1
            # rotate
            h = (f0 - vertex @ f) * g - (g0 - vertex @ g) * f
            h0 = (f0 - vertex @ f) * g0 - (g0 - vertex @ g) * f0
            # update g
            g = np.copy(h)
            g0 = np.copy(h0)
            # find the new vertex
            products = polytope.deterministics @ g
            idx = np.where(products == np.amax(products))[0][0]
            vertex = polytope.deterministics[idx]
            # stop iteration
            if vertex @ g == g0:
                new_ineq = np.r_[h, -1.0 * h0]
                return polytope_from_inequality(new_ineq, polytope)
        # maximum number of rotations reached
        print('Rotation got stuck')
        return False


def polytope_from_inequality(ineq, poly):
    """ Generates a polytope from an inequality """
    sub_dets = np.array([v for v in poly.deterministics if v @ ineq[:-1] == -1.0 * ineq[-1]])
    # find possible relabelling
    sub_relabellings = []
    for r in poly.relabellings:
        allowed = False
        for sd in sub_dets:
            if np.any(np.sum((sub_dets - sd[r]) ** 2, axis=1) < 1e-6):
                allowed = True
                break
        if allowed:
            sub_relabellings.append(r)
    sub_relabellings = np.array(sub_relabellings)
    # create new polytope
    return Polytope(sub_dets, sub_relabellings, creating_face=ineq, parent=poly, initial_polytope=poly.initial_polytope)
