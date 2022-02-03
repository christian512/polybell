""" Here I collect helper functions for the panda vertex enumeration code """
import numpy as np
from fractions import Fraction
import subprocess


def write_panda_input_inequalities(lhs, rhs, idx_equalities=[], dets=None, symmetries=None, file='', denom_limit=10000):
    """
    the inequality is interpreted as lhs * x <= rhs
    lhs: left hand side of inequalities
    rhs: righ hand side of inequalities
    idx_equalities: which entries of the inequality system should be equalities
    dets: deterministic points for equivalence checking in panda
    symmetries: the symmetries that we want to consider given in a list of lists
                each list tells us how the variables are interchanged.
    """
    # check dimensions
    if type(symmetries) != type(None):
        symmetries = np.array(symmetries)
    assert lhs.shape[0] == rhs.shape[0]
    # string to write out
    string = ""
    # variable names (only needed for symmetries)
    string += 'Names:\n'
    for i in range(lhs.shape[1]):
        string += 'a' + str(i) + ' '
    string += '\n'
    # add symmetries
    if type(symmetries) != type(None):
        string += 'Maps:\n'
        for symm in symmetries:
            for i in range(symm.shape[0]):
                string += 'a' + str(symm[i]) + ' '
            string += '\n'
    # equalities
    if idx_equalities:
        string += 'Equations:\n'
        for i in idx_equalities:
            arr = [Fraction(x).limit_denominator(denom_limit).denominator for x in lhs[i]]
            arr.append(Fraction(rhs[i]).limit_denominator(denom_limit).denominator)
            factor = np.product(np.unique(arr))
            for j, x in enumerate(lhs[i]):
                string += str(Fraction(factor * x).limit_denominator(denom_limit)) + 'a' + str(j) + ' '
            # TODO: check if this = is correct or can I just drop it?
            string += '= ' + str(Fraction(factor * rhs[i]).limit_denominator(denom_limit)) + '\n'
    # inequalities -> as we can not give fractions we multiply the inequality by the product of unique denominators
    string += 'Inequalities:\n'
    for i in range(lhs.shape[0]):
        if i in idx_equalities:
            continue
        # find factor to multiply inequality with, such that we have a only integers inequality
        arr = [Fraction(x).limit_denominator(denom_limit).denominator for x in lhs[i]]
        arr.append(Fraction(rhs[i]).limit_denominator(denom_limit).denominator)
        factor = np.product(np.unique(arr))
        # write the inequalitiy
        for j in range(lhs.shape[1]):
            val = int(lhs[i, j] * factor)
            string += str(val) + 'a' + str(j) + ' '
        string += '<= ' + str(int(rhs[i] * factor)) + '\n'

    # deterministics given to file
    if type(dets) != type(None):
        string += 'Deterministics:'
        for d in dets:
            string += '\n'
            for val in d:
                string += str(int(val)) + ' '

    # write file
    if file:
        f = open(file, 'w+')
        f.write(string)
        f.close()
    return string


def write_known_vertices(vertices, file='knownvertices.ext', relabellings_vertices=[], denom_limit=9999999):
    """ Writes known vertices to a file """
    s = ""
    if len(relabellings_vertices) > 0:
        s += "VERTEXMAPS: \n"
        for r in relabellings_vertices:
            for val in r:
                s += str(int(val)) + ' '
            s += '\n'
    s += "Vertices: \n"
    for i in range(vertices.shape[0]):
        arr = [Fraction(x).limit_denominator(denom_limit) for x in vertices[i]]
        for val in arr:
            s += str(val) + ' '
        s += '\n'
    f = open(file, 'w+')
    f.write(s)
    f.close()
    return s


def write_known_inequalities(lhs, rhs, file='knowninequalities.ine', denom_limit=100000):
    """ Writes known inequalities to a file """
    s = "Inequalities: \n"
    assert lhs.shape[0] == rhs.shape[0]
    for i in range(lhs.shape[0]):
        # find factor to multiply inequality with, such that we have a only integers inequality
        arr = [Fraction(x).limit_denominator(denom_limit).denominator for x in lhs[i]]
        arr.append(Fraction(rhs[i]).limit_denominator(denom_limit).denominator)
        factor = np.product(np.unique(arr))
        s += '  '
        # write the inequalitiy
        for j in range(lhs.shape[1]):
            val = int(lhs[i, j] * factor)
            s += str(val) + '  '
        s += str(int(-1 * rhs[i] * factor)) + '\n'
    f = open(file, 'w+')
    f.write(s)
    f.close()
    return s


def read_vertices_rays(file):
    """ Reads vertices and rays from a PANDA output file """
    f = open(file, 'r')
    string = f.read()
    f.close()
    vertices_string = string.split('\n')[1:-1]
    vertices = []
    rays = []
    for v_string in vertices_string:
        arr = np.fromstring(v_string[1:], sep=' ', dtype=float)
        if arr[-1] == 0:
            rays.append(arr[:-1])
        else:
            v = arr[:-1] / arr[-1]
            vertices.append(v)
    return np.array(vertices), np.array(rays)


def run_panda(file, threads=4, outfile='', known_data=''):
    """ Runs panda on a file"""
    cmd = 'panda ' + file + ' -t ' + str(threads)
    if known_data:
        cmd += ' -k ' + known_data
    if outfile:
        cmd += ' > ' + outfile
    out = subprocess.run(cmd, shell=True)
    return out


def run_panda_vertices_on_facet(file, outfile='', known_data=''):
    """ Catches the output given about the vertices on facet """
    cmd = 'panda ' + file + ' -t 1'
    if known_data:
        cmd += ' -k ' + known_data
    if outfile:
        cmd += ' > ' + outfile
    out = subprocess.run(cmd, shell=True, capture_output=True)
    lines = out.stderr.decode('utf-8').split('\n')
    counter = 0
    while not "Vertices on facet" in lines[counter]:
        counter += 1
    vertices = []
    counter += 1
    while not "Stopping" in lines[counter]:
        line = lines[counter].strip()
        vertex = np.fromstring(line, sep=' ', dtype=float)
        vertices.append(vertex)
        counter += 1
    vertices = np.array(vertices)
    return vertices

def read_inequalities(file):
    """ READS INEQUALITIES FROM A FILE """
    f = open(file, 'r')
    lines = f.readlines()
    f.close()
    # run until inequalities section
    counter = 0
    while not 'Inequalities' in lines[counter]:
        counter += 1
    counter += 1
    facets = []
    for line in lines[counter:]:
        arr = np.fromstring(line[1:], sep=' ', dtype=float)
        facets.append(arr)
    return np.array(facets, dtype=int)

