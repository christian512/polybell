import argparse
from linearbell.utils import get_deterministic_behaviors, find_local_weight_dual, get_configs, general_pr_box, \
    check_equiv_bell_vertex_enum_non_rescale
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
    global recursion_depth
    recursion_depth += 1
    print('recursion_depth ', recursion_depth)
    print('len vertices ', len(vertices))
    # TODO: IS the rotation algorithm correct?
    # TODO: IF the number here is to small it does not work -> Why?
    if vertices.shape[0] <= 2:
        if vertices.shape[0] == 2:
            return np.array([facet])
        write_known_vertices(vertices[:, :-1], file='../vertex_enum_panda/knownvertices.ext')
        # run original panda to get all facets
        cmd = 'panda_org knownvertices.ext -t 1 --method=dd > out.ine'
        out = subprocess.run(cmd, shell=True)
        all_facets = read_inequalities('../vertex_enum_panda/out.ine')
        assert facet in all_facets
        ineq_facets = [facet]
        for f in all_facets:
            equiv = False
            for c in ineq_facets:
                if recursion_depth <= 0:
                    if check_equiv_bell_vertex_enum_non_rescale(f[:-1], c[:-1], relabels, dets):
                        equiv = True
                else:
                    if np.all(c == f):
                        equiv = True
            if not equiv:
                ineq_facets.append(f)
        return np.array(ineq_facets)

    # list of inequivalent facets
    ineq_facets = [facet]
    # facets to process
    check_facets = [facet]
    # calculate furthest vertex
    while check_facets:
        curr_facet = check_facets.pop()
        furthestVertex = furthest_vertex(vertices, curr_facet)
        # get vertices on that face
        subvertices = np.array([v for v in vertices if distance(v, curr_facet) == 0])
        write_known_vertices(subvertices[:, :-1], file='../vertex_enum_panda/input.ext')
        # run panda to get one facet
        run_panda_vertices_on_facet('../vertex_enum_panda/input.ext', outfile='out.ine')
        subfacet = read_inequalities('../vertex_enum_panda/out.ine')[0]
        # get the ridges
        ridges = get_all_ineq_facets(subvertices, subfacet)
        # iterate over ridges
        for ridge in ridges:
            # TODO: I guess the rotation algorithm is not correct! If we would get only valid faces, we would always get 2
            new_facet = rotate(vertices, furthestVertex, curr_facet, ridge)
            if np.all(new_facet == 0):
                continue
            equiv = False
            for c in ineq_facets:
                assert c.shape[0] == new_facet.shape[0]
                assert relabels.shape[1] == c.shape[0] - 1
                # Only do equivalence check above some recursion level.
                if recursion_depth <= 2:
                    if check_equiv_bell_vertex_enum_non_rescale(new_facet[:-1], c[:-1], relabels, dets):
                        equiv = True
                else:
                    if np.all(c == new_facet):
                        equiv = True
            if not equiv:
                ineq_facets.append(new_facet)
                check_facets.append(new_facet)

    recursion_depth -= 1
    return np.array(ineq_facets)

# Write the current vertices
write_known_vertices(dets, file='../vertex_enum_panda/input.ext')
# run panda to get one facet
run_panda_vertices_on_facet('../vertex_enum_panda/input.ext', outfile='out.ine')
initial_facet = read_inequalities('../vertex_enum_panda/out.ine')[0]

facets = get_all_ineq_facets(dets_extended, initial_facet)
# TODO: HERE YOU HAVE TO START AGAIN, WITH EACH FACET THAT WAS NOT USED
print('Number of facets: ', facets.shape[0])
