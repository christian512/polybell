""" A class to describe a polytope """
import numpy as np
import subprocess
import string
import random
from linearbell.panda_helper import run_panda_vertices_on_facet, read_inequalities, write_known_vertices
from linearbell.utils import equiv_check_adjacency_panda
from bokeh.models import Circle, MultiLine, Range1d
from bokeh.plotting import figure, from_networkx
import networkx as nx


class ParamPolytope():
    """ Describing a polytope """

    def __init__(self, vertices, permutations_vertices, permutations_coordinates, creating_face=np.array([]), indices_vertices=[], parent=None,
                 initial_polytope=None):
        """

        Parameters
        ----------
        vertices : Vertices of the polytope
        permutations_vertices : Permutation which vertex can be changed to which under relabelling
        creating_face : Ray description of the face that created this polytope
        parent: Parent polytope
        initial_polytope: first polytope
        """
        self.id = ''.join(random.choices(string.digits + string.ascii_letters, k=30))
        self.vertices = vertices
        self.permutations_vertices = permutations_vertices
        self.permutations_coordinates = permutations_coordinates
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
        self.dims = np.linalg.matrix_rank(self.vertices)
        # empty faces
        self._faces = []
        # empty classes
        self._classes = []
        # indices of the vertices
        self.indices_vertices = indices_vertices
        if not self.indices_vertices:
            for det in self.vertices:
                equal = det == self.initial_polytope.vertices
                equal = np.all(equal, axis=1)
                idx = np.where(equal)[0][0]
                self.indices_vertices.append(idx)

    def __str__(self):
        return "polytope {} with {} vertices  ".format(self.id, self.vertices.shape[0])

    def __eq__(self, other):
        if np.all(self.vertices == other.vertices):
            return True
        return False

    def add_faces(self,faces):
        """ Adds faces to the polytope """
        for f in faces:
            if f not in self._faces:
                self._faces.append(f)
                if not f.equiv_under_parent(self._classes):
                    self._classes.append(f)
        return True


    def get_single_face(self):
        """ Returns a single face """
        if len(self._faces) != 0:
            return self._faces[0]
        # Run panda to get faces
        write_known_vertices(self.vertices, file='knownvertices.ext')
        cmd = 'panda knownvertices.ext -t 1 > out.ine'
        out = subprocess.run(cmd, shell=True)
        inequality = read_inequalities('out.ine').astype(float)[0]
        return parampolytope_from_inequality(inequality, self)

    def get_all_faces(self):
        """ Finds all faces of the polytope """
        if len(self._faces) != 0:
            return self._faces
        # Run panda to get faces
        write_known_vertices(self.vertices, file='knownvertices.ext')
        cmd = 'panda_org knownvertices.ext -t 1 --method=dd > out.ine'
        out = subprocess.run(cmd, shell=True)
        inequalities = read_inequalities('out.ine').astype(float)
        # generate polytopes from faces
        for ineq in inequalities:
            p = parampolytope_from_inequality(ineq, self)
            self._faces.append(p)
        # also set classes
        self.get_all_classes()
        return self._faces

    def get_all_classes(self):
        """ Finds all classes of faces of a polytope """
        if len(self._classes) != 0:
            return self._classes
        if len(self._faces) == 0:
            self.get_all_faces()
        # Reduce the faces to inequivalent ones
        self._classes = [self._faces[0]]
        for f in self._faces:
            if not f.equiv_under_parent(self._classes):
                self._classes.append(f)
        return self._classes

    def rotate(self, ridge):
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
        assert polytope.vertices.shape[1] == f.shape[0]
        products = polytope.vertices @ f
        assert products.shape[0] == polytope.vertices.shape[0]
        idx = np.where(products == np.amin(products))[0][0]
        vertex = polytope.vertices[idx]
        assert vertex @ f < f0
        # start the iteration
        counter = 0
        while counter < 1000:
            counter += 1
            # rotate
            h = (f0 - vertex @ f) * g - (g0 - vertex @ g) * f
            h0 = (f0 - vertex @ f) * g0 - (g0 - vertex @ g) * f0
            # update g and divide by GCD
            gcd = np.gcd.reduce(np.r_[h.astype(int), int(h0)])
            g = h / gcd
            g0 = h0 / gcd
            # find the new vertex
            products = polytope.vertices @ g
            products = np.round(products, decimals=5)
            idx = np.where(products == np.amax(products))[0][0]
            vertex = polytope.vertices[idx]
            # stop iteration
            if vertex @ g == g0:
                new_ineq = np.r_[h, -1.0 * h0]
                new_poly = parampolytope_from_inequality(new_ineq, polytope)
                return new_poly
        # maximum number of rotations reached
        print('Rotation got stuck')
        return False

    def equiv_under_parent(self, other):
        """ Checks if two polytopes are equivalent under the symmetry of the parent """
        if type(other) == list:
            for o in other:
                if self.equiv_under_parent(o):
                    return o
            return False
        assert self.parent == other.parent, 'Trying to check equivalence for polytopes with unequal parents'
        if self.vertices.shape[0] != other.vertices.shape[0]:
            return False
        for r in self.parent.permutations_vertices:
            # check if the polytopes have the same vertices under the relabelling
            d = dict(enumerate(r))
            new_vertices_indices = sorted([d[x] for x in self.indices_vertices])
            if sorted(other.indices_vertices) == new_vertices_indices:
                return True
        return False

    def equiv_under_initial(self, other):
        """ Checks if two polytopes are equivalent under the symmetry of the parent """
        if type(other) == list:
            for o in other:
                res = self.equiv_under_initial(o)
                if res:
                    return o
            return False
        assert self.initial_polytope == other.initial_polytope, 'Trying to check equivalence for polytopes with unequal parents'
        if self.vertices.shape[0] != other.vertices.shape[0]:
            return False
        for r in self.initial_polytope.permutations_vertices:
            # check if the polytopes have the same vertices under the relabelling
            d = dict(enumerate(r))
            new_vertices_indices = sorted([d[x] for x in other.indices_vertices])
            if sorted(self.indices_vertices) == new_vertices_indices:
                return True
        return False


def parampolytope_from_inequality(ineq, poly):
    """ Generates a polytope from an inequality """
    vertices_on_face = []
    indices_vertices_on_face = []
    # Find the vertices that equalise the face
    for i in range(poly.vertices.shape[0]):
        if poly.vertices[i] @ ineq[:-1] == -1.0 * ineq[-1]:
            vertices_on_face.append(poly.vertices[i])
            indices_vertices_on_face.append(poly.indices_vertices[i])
    # Find which permutations of vertices are still possible
    permutation_vertices_on_face = []
    for r in poly.permutations_vertices:
        d = dict(enumerate(r))
        valid_permutation = True
        for i in indices_vertices_on_face:
            if not d[i] in indices_vertices_on_face:
                valid_permutation = False
                break
        if not valid_permutation:
            continue
        permutation_vertices_on_face.append(r)
    vertices_on_face = np.array(vertices_on_face)
    if len(vertices_on_face) == 0:
        print('Error no vertex on face')
        return None
    permutation_vertices_on_face = np.array(permutation_vertices_on_face)
    permutation_coordinates_on_face = np.array([])
    # create new polytope
    return ParamPolytope(vertices_on_face, permutation_vertices_on_face, permutation_coordinates_on_face, creating_face=ineq,
                         indices_vertices=indices_vertices_on_face, parent=poly,
                         initial_polytope=poly.initial_polytope)
