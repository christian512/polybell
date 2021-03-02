""" A class to describe a polytope """
import numpy as np
import subprocess
import string
import random
from linearbell.panda_helper import run_panda_vertices_on_facet, read_inequalities, write_known_vertices
from linearbell.utils import equiv_check_adjacency_panda
from bokeh.io import show, save
from bokeh.models import Circle, MultiLine, Range1d
from bokeh.plotting import figure, from_networkx
import networkx as nx


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
                new_poly = polytope_from_inequality(new_ineq, polytope)
                return new_poly
        # maximum number of rotations reached
        print('Rotation got stuck')
        return False

    def equiv_under_parent(self, other):
        """ Checks if two polytopes are equivalent under the symmetry of the parent """
        assert self.parent == other.parent, 'Parents of the polytopes to compare are different'
        return equiv_check_adjacency_panda(self.creating_face, other.creating_face, self.parent.relabellings,
                                           self.parent.deterministics)

    def equiv_under_bell(self, other):
        """ Checks if two polytopes are equivalent under the symmetry of the Bell Polytope """
        # if multiple other polytopes are given
        if type(other) == list:
            for o in other:
                if self.equiv_under_bell(o):
                    return True
            return False
        # Here starts if other is an actual polytope
        assert self.initial_polytope == other.initial_polytope, 'Initial Polytope of both polytopes is not the same'
        return equiv_check_adjacency_panda(self.creating_face, other.creating_face, self.initial_polytope.relabellings,
                                           self.initial_polytope.deterministics)

    def get_valid_face(self, poss_face):
        """ Checks if a given face is a valid face """
        # Check by Parent
        if poss_face.parent == self:
            return poss_face
        # the dimensions is one less than the self dimension
        if self.dims - 1 != poss_face.dims:
            print('Dimensions of Possible face are not matching')
            return False
        # get the deterministics that equalize the creating face
        ineq = poss_face.creating_face
        sub_dets = np.array([v for v in self.deterministics if v @ ineq[:-1] == -1.0 * ineq[-1]])
        if self.dims - 1 == np.linalg.matrix_rank(sub_dets) - 1:
            return polytope_from_inequality(ineq, self)
        return False

    def get_valid_face_relabelling(self, poss_face):
        """
        Checks if a given face is a valid face under every relabelling of the polytope
        For the relabelling we can use the Relabellings of the Bell Polytope
         """
        if poss_face.parent == self:
            return poss_face
        if self.dims - 1 != poss_face.dims:
            print('Dimensions of possible face are not matching')
            return False
        for r in self.relabellings:
            # relabel creating face of the possible face
            relabelled_ineq = poss_face.creating_face[r]
            relabelled_ineq = np.r_[relabelled_ineq, poss_face.creating_face[-1]]
            # take the deterministic points that equalize this -> this way the deterministic points are only from polytope
            sub_dets = np.array(
                [v for v in self.deterministics if v @ relabelled_ineq[:-1] == -1.0 * relabelled_ineq[-1]])
            # check dimensions -> if they are correct the face contains only valid dets -> is a valid face
            if self.dims - 1 == np.linalg.matrix_rank(sub_dets) - 1:
                return polytope_from_inequality(relabelled_ineq, self)
        return False


def polytope_from_inequality(ineq, poly):
    """ Generates a polytope from an inequality """
    sub_dets = np.array([v for v in poly.deterministics if v @ ineq[:-1] == -1.0 * ineq[-1]])
    # find possible relabelling
    sub_relabellings = []
    for r in poly.relabellings:
        allowed = True
        for sd in sub_dets:
            if np.all(np.sum((sub_dets - sd[r]) ** 2, axis=1) > 1e-6):
                allowed = False
                break
        if allowed:
            sub_relabellings.append(r)
    sub_relabellings = np.array(sub_relabellings)
    # create new polytope
    return Polytope(sub_dets, sub_relabellings, creating_face=ineq, parent=poly, initial_polytope=poly.initial_polytope)

def print_face_lattice(all_polys):
    """ Prints all polytopes given in a dict, where the keys are the levels """
    # Create networkx graph
    G = nx.Graph()
    # add the nodes
    for level in all_polys.keys():
        for j in range(len(all_polys[level])):
            p = all_polys[level][j]
            G.add_node(p.id, pos=(j, -level), ndets=len(p.deterministics),
                       nrel=len(p.relabellings),
                       dims=p.dims, dets_indices=str(p.indices_deterministics), face=str(p.creating_face))
    pos = nx.get_node_attributes(G, 'pos')
    network_graph = from_networkx(G, pos)
    # Set node size and color
    network_graph.node_renderer.glyph = Circle(size=15, fill_color='skyblue')

    # Set edge opacity and width
    network_graph.edge_renderer.glyph = MultiLine(line_alpha=0.5, line_width=1)
    HOVER_TOOLTIPS = [("Number Deterministics", "@ndets"), ("Indices Deterministics", "@dets_indices"),
                      ("Number relabels", "@nrel"), ("Dimensions", "@dims"), ("Face", "@face")]
    plot = figure(tooltips=HOVER_TOOLTIPS, x_range=Range1d(0, 20), y_range=Range1d(-8, 2),
                  title='Face-Classes-Lattice for 2222 case')
    plot.renderers.append(network_graph)
    return plot