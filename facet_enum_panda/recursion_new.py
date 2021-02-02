import argparse
from linearbell.utils import get_deterministic_behaviors,get_all_inequiv_bell_non_rescale
from linearbell.panda_helper import write_known_vertices, read_inequalities, run_panda_vertices_on_facet
from linearbell.adjacency_decomposition import rotate, furthest_vertex, distance
import numpy as np
import subprocess

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
    global recursion_depth
    recursion_depth += 1
    print('len vertices ', len(vertices))
    # if number of vertices is small enough -> enumerate with double description
    if vertices.shape[0] <= 4:
        write_known_vertices(vertices[:, :-1], file='knownvertices.ext')
        # run original panda to get all facets
        cmd = 'panda_org knownvertices.ext -t 1 --method=dd > out.ine'
        out = subprocess.run(cmd, shell=True)
        all_facets = read_inequalities('out.ine')
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
        subvertices = np.array([v for v in vertices if distance(v, curr_facet) == 0])
        assert subvertices.shape[0] < vertices.shape[0]
        write_known_vertices(subvertices[:, :-1], file='input.ext')
        # run panda to get one facet
        run_panda_vertices_on_facet('input.ext', outfile='out.ine')
        subfacet = read_inequalities('out.ine')[0]
        # get the ridges
        ridges = get_all_ineq_facets(subvertices, subfacet)
        # iterate over ridges
        for ridge in ridges:
            # apply rotation algorithm
            new_facet = rotate(vertices, furthestVertex, curr_facet, ridge)
            # check that new face is not 0
            if np.all(new_facet == 0):
                continue
            # if the new facet is not yet there, we also check that one
            if np.sum((np.array(ineq_facets) - new_facet)**2, axis=1).all() > 0:
                ineq_facets.append(new_facet)
                check_facets.append(new_facet)

    recursion_depth -= 1
    return np.array(ineq_facets)

# Write the current vertices
write_known_vertices(test_vertices, file='input.ext')
# run panda to get one facet
run_panda_vertices_on_facet('input.ext', outfile='out.ine')
initial_facet = read_inequalities('out.ine')[0]

facets = get_all_ineq_facets(test_vertices_extended, initial_facet)
print('Number of facets: ', facets.shape[0])
