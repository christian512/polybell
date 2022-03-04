import argparse
from polybell.utils import get_deterministic_behaviors, check_equiv_bell_vertex_enum_non_rescale, \
    equiv_check_adjacency_testing
from polybell.panda_helper import write_known_vertices, read_inequalities, run_panda_vertices_on_facet
from polybell.adjacency_decomposition import rotate, furthest_vertex, distance
import numpy as np
import subprocess
import time

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument(dest='ma', help="number of inputs for ALICE")
parser.add_argument(dest='mb', help='number of inputs for BOB')
parser.add_argument(dest='n', help='number of outputs')
args = parser.parse_args()
ma = int(args.ma)
mb = int(args.mb)
n = int(args.n)
inputs_a = range(ma)
inputs_b = range(mb)
outputs = range(n)

# INITIALIZE BELL POLYTOPE
dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs)
# dets = dets - p_origin
dets = dets[np.lexsort(np.rot90(dets))]

# add a one column to vertices to match panda format
dets_extended = np.c_[dets, np.ones(dets.shape[0])]

relabels = np.loadtxt('../data/relabels/{}{}{}{}.gz'.format(ma, mb, n, n)).astype(int)

test_vertices = np.array(
    [[1, 1, 1], [-1, 1, 1], [1, -1, 1], [-1, -1, 1], [1, 1, -1], [-1, 1, -1], [1, -1, -1], [-1, -1, -1]])
# append array of ones to vertices to match PANDA description
test_vertices_extended = np.c_[test_vertices, np.ones(test_vertices.shape[0])]

# recursion depth counter
recursion_depth = -1


def reduce_to_inequiv(bells, relabels, dets, tol=1e-6):
    """
    Returns an array of inequivalent Bell expressions
    Parameters
    ----------
    bells
    relabels
    dets
    tol

    Returns
    -------
    """
    ineq_bells = [bells[0]]
    for i, b in enumerate(bells):
        # print('equiv check progress : {} / {}'.format(i, bells.shape[0]))
        equiv = False
        for c in ineq_bells:
            if equiv_check_adjacency_testing(b[:-1], c[:-1], relabels, dets, tol=tol):
                equiv = True
                break
        if not equiv:
            ineq_bells.append(b)
    return np.array(ineq_bells)


def facet_equiv_to_any(bell, bells_arr, relabels, dets, tol=1e-6):
    """ Checks if a bell expression is equivalent to any bell expression in a list. """
    # rescale bell
    divisor = bell[-1]
    if divisor < 0:
        bell = bell[:-1] / np.abs(divisor)
    else:
        assert divisor == 0
        bell = bell[:-1] + 1 / (len(outputs) ** 2)


    for i in range(bells_arr.shape[0]):
        divisor = bells_arr[i, -1]
        if divisor < 0:
            other_bell = bells_arr[i, :-1] / np.abs(divisor)
        else:
            assert divisor == 0
            other_bell = bells_arr[i, :-1] + 1 / (len(outputs) ** 2)
        if equiv_check_adjacency_testing(bell, other_bell, relabels, dets):
            return True
    return False


def get_poss_relabels(subvertices, curr_relabels):
    # do this with numpy functions to use parallel code
    # cut off the last entry
    subvertices = subvertices[:, :-1]
    new_relabels = []
    for r in curr_relabels:
        for subvertex in subvertices:
            if np.any(np.sum((subvertices - subvertex[r]) ** 2, axis=1) == 0):
                new_relabels.append(r)
                break
    return np.array(new_relabels)


# define the recursive function
def get_all_ineq_facets(vertices, facet):
    """
    Gets all inequivalent facets by adjacency decomposition
    Parameters
    ----------
    vertices: vertices of the polytope
    facet: one facet of the polytope
    Returns
    -------
    """
    # print('len vertices ', len(vertices))
    # if number of vertices is small enough -> enumerate with double description
    if vertices.shape[0] <= 24:
        write_known_vertices(vertices[:, :-1], file='knownvertices.ext')
        # run original panda to get all facets
        cmd = 'panda_org knownvertices.ext -t 1 --method=dd > out.ine'
        out = subprocess.run(cmd, shell=True)
        all_facets = read_inequalities('out.ine')
        print('Finished PANDA and found {} facets on the subpolytope'.format(all_facets.shape[0]))
        assert facet in all_facets
        return all_facets

    # list of inequivalent facets
    ineq_facets = [facet]
    # facets to process
    check_facets = [facet]
    # calculate furthest vertex
    while check_facets:
        curr_facet = check_facets.pop()
        furthestVertex = furthest_vertex(vertices, curr_facet)
        assert distance(furthestVertex, curr_facet) != 0
        # get vertices on that face
        subvertices = np.array([v for v in vertices if distance(v, curr_facet) == 0], dtype=int)
        assert subvertices.shape[0] < vertices.shape[0]
        write_known_vertices(subvertices[:, :-1], file='input.ext')
        # run panda to get one facet
        run_panda_vertices_on_facet('input.ext', outfile='out.ine')
        subfacet = read_inequalities('out.ine')[0]
        # get the ridges
        ridges = get_all_ineq_facets(subvertices, subfacet)
        if vertices.shape[0] == dets.shape[0]:
            print('TOP LEVEL: Num ineq facets: {}'.format(len(ineq_facets)))
            print('TOP LEVEL: Num facets to check: {}'.format(len(check_facets)))
        # iterate over ridges
        for i in range(ridges.shape[0]):
            # apply rotation algorithm
            new_facet = rotate(vertices, furthestVertex, curr_facet, ridges[i])
            # check if we're at the top level and check for neighbouring
            if vertices.shape[0] == dets.shape[0]:
                if not facet_equiv_to_any(new_facet, np.array(ineq_facets), relabels, dets):
                    ineq_facets.append(new_facet)
                    check_facets.append(new_facet)
                    np.savetxt('../data/facets_panda/{}{}{}{}.txt'.format(ma, mb, n, n), np.array(ineq_facets))
            else:
                # if not at top level only rotate and add as facets
                ineq_facets.append(new_facet)
    return np.array(ineq_facets)


# Write the current vertices
write_known_vertices(dets, file='input.ext')
# run panda to get one facet
run_panda_vertices_on_facet('input.ext', outfile='out.ine')
initial_facet = read_inequalities('out.ine')[0]
start = time.time()
facets = get_all_ineq_facets(dets_extended, initial_facet)
print('Calculation time: ', time.time() - start)
print('Number of facets: ', facets.shape[0])
