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

    def __init__(self, deterministics, relabellings, creating_face=np.array([]), parent=None):
        """ initialize polytope """
        self.id = ''.join(random.choices(string.digits + string.ascii_letters, k=30))
        self.deterministics = deterministics
        self.relabellings = relabellings
        # set the face that created this
        self.creating_face = creating_face
        # set parent
        self.parent = parent
        if not self.parent:
            self.parent == self
        # dimensions
        self.dims = np.linalg.matrix_rank(self.deterministics) - 1
        # empty faces
        self.__faces = []
        # empty classes
        self.__classes = []

    def __str__(self):
        return "polytope {} with {} deterministics and {} children ".format(self.id, self.deterministics.shape[0],
                                                                            len(self.children))

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
        """ Rotates the polytope around a child.This polytope is the face, and it's parent is the 'real' polytope """
        if self.parent == self:
            print('Can not rotate as this polytope does not have a parent')
        assert ridge in self.children
        assert self in self.parent.children
        # set the structure
        polytope = self.parent
        face = self
        # do rotation
        vertices = np.c_[polytope.deterministics, np.ones(polytope.deterministics.shape[0])]
        vertex = furthest_vertex(vertices, face.creating_face)
        new_face = rotate(vertices, vertex, face.creating_face, ridge.creating_face)
        # create a new polytope
        return polytope_from_inequality(new_face, polytope)

def polytope_from_inequality(ineq, poly):
    """ Generates a polytope from an inequality """
    sub_dets = np.array([v for v in poly.deterministics if v @ ineq[:-1] == -1.0 * ineq[-1]])
    # find possible relabelling
    sub_poss_relabellings = []
    for r in poly.poss_relabellings:
        allowed = False
        for sd in sub_dets:
            if np.any(np.sum((sub_dets - sd[r]) ** 2, axis=1) < 1e-6):
                allowed = True
                break
        if allowed:
            sub_poss_relabellings.append(r)
    sub_poss_relabellings = np.array(sub_poss_relabellings)
    # create new polytope
    return Polytope(sub_dets, sub_poss_relabellings, creating_face=ineq, parent=poly)
