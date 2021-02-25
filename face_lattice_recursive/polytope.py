""" A class to describe a polytope """
import numpy as np
import string
import random
from linearbell.panda_helper import run_panda_vertices_on_facet, read_inequalities, write_known_vertices
from linearbell.adjacency_decomposition import furthest_vertex, rotate


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
        # Children of the polytope, which are faces with respect to some symmetry
        self.children = []
        # dimensions
        self.dims = np.linalg.matrix_rank(self.deterministics) - 1
        # indices deterministics

    def __str__(self):
        return "polytope {} with {} deterministics and {} children ".format(self.id, self.deterministics.shape[0],
                                                                            len(self.children))

    def add_child(self, poly):
        """ Adds a child as polytope """
        # check if the deterministic points of child are subgroup of the ones that we have
        if not self.check_valid_face(poly):
            return False
        if poly in self.children:
            print('Child already stored')
            return False
        self.children.append(poly)
        return True

    def find_single_face(self):
        """ Finds a single face of the polytope, creates a new polytope and returns it """
        # run the modified version of panda
        write_known_vertices(self.deterministics, file='input.ext')
        run_panda_vertices_on_facet('input.ext', outfile='out.ine')
        face = read_inequalities('out.ine')[0]
        # get vertices that are on this facet
        face_dets = np.array([v for v in self.deterministics if v @ face[:-1] == -1.0 * face[-1]])
        sub_relabellings = []
        for r in self.relabellings:
            allowed = False
            for sd in face_dets:
                if np.any(np.sum((face_dets - sd[r]) ** 2, axis=1) < 1e-6):
                    allowed = True
                    break
            if allowed:
                sub_relabellings.append(r)
        sub_relabellings = np.array(sub_relabellings)
        # create new polytope
        return Polytope(face_dets, sub_relabellings, creating_face=face, parent=self)

    def rotate_polytope(self, ridge):
        """ Rotates the polytope around a child.This polytope is the face, and it's parent is the 'real' polytope """
        if self.parent == self:
            print('Can not rotate as this polytope does not have a parent')
        assert ridge in self.children
        assert self in self.parent.children
        # set the structure
        polytope = self.parent
        face = self
        # get furthest vertex
        vertices = np.c_[polytope.deterministics, np.ones(polytope.deterministics.shape[0])]
        vertex = furthest_vertex(vertices, face.creating_face)
        new_face = rotate(vertices, vertex, face.creating_face, ridge.creating_face)
        # create a new polytope
        face_dets = np.array([v for v in polytope.deterministics if v @ new_face[:-1] == -1.0 * new_face[-1]])
        sub_relabellings = []
        for r in polytope.relabellings:
            allowed = False
            for sd in face_dets:
                if np.any(np.sum((face_dets - sd[r]) ** 2, axis=1) < 1e-6):
                    allowed = True
                    break
            if allowed:
                sub_relabellings.append(r)
        sub_relabellings = np.array(sub_relabellings)
        return Polytope(face_dets, sub_relabellings, creating_face=new_face, parent=polytope)

    def check_valid_face(self, face):
        """ Checks if a given face is valid for this polytope """
        for det in face.deterministics:
            equal = det == self.deterministics
            equal = np.all(equal, axis=1)
            if np.all(equal == False):
                print('Not a valid face, due to deterministics')
                return False
        return True
