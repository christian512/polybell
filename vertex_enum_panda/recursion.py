"""
Here we want to test the idea of a recursive Adjacency Decomposition method.

I setup PANDA to finish, when it returns the points on the current facet.
We will run PANDA as facet enumerator and give one facet as input to PANDA.
This one facet can be found by GurobiPy linear programming.
"""
import argparse
from linearbell.utils import get_deterministic_behaviors, find_local_weight_dual, get_configs, general_pr_box
from linearbell.panda_helper import write_known_vertices, read_inequalities, run_panda_vertices_on_facet
import numpy as np

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
p_origin = np.sum(dets, axis=0) / dets.shape[0]
dets = dets - p_origin
dets = dets[np.lexsort(np.rot90(dets))]
print('Number of dets: {}'.format(dets.shape[0]))
#relabels = np.loadtxt('../data/relabels/{}{}{}{}.gz'.format(ma, mb, n, n)).astype(int)

# setup pr box
configs = get_configs(inputs_a, inputs_b, outputs, outputs)
pr_box = np.array([general_pr_box(*c) for c in configs])

# storage for all facets
facets = []
all_subvertices = []

# WRITE OUT VERTICES
write_known_vertices(dets, file='input.ext')
# Run Panda and store facets
subvertices = run_panda_vertices_on_facet('input.ext', outfile='out.ine')
all_subvertices.append(subvertices)


print('Num subvertices: {}'.format(subvertices.shape[0]))

while subvertices.shape[0] > 1:
    # TODO: SOMEHOW ALWAYS ALL RELABELS ARE POSSIBLE
    # STORE FACETS
    curr_facets = read_inequalities('out.ine')
    assert curr_facets.shape[0] == 1
    print(curr_facets[0])
    facets.append(curr_facets[0])
    # WRITE KNOWN VERTICES
    write_known_vertices(subvertices, file='input.ext')
    subvertices = run_panda_vertices_on_facet('input.ext', outfile='out.ine')
    print('Num subvertices: {}'.format(subvertices.shape[0]))
    all_subvertices.append(subvertices)

assert len(facets) == len(all_subvertices) - 1

# As the last subfacet is a vertex we rotate the last facet around this vertex
facet = facets.pop()
subfacet = all_subvertices[-1]
# TODO: Rotate facet around this
new_facets = rotate(facet, [subfacet])
# TODO: Do equivalence check, to only take inequivalent faces
new_facets = equivalence_check(new_facets, relabels)


while facets:
    print('num new facets: {}'.format(len(new_facets)))
    print('remaining steps: {}'.format(len(facets)))
    facet = facets.pop()
    # TODO: Rotate facet around ridges
    new_facets = rotate(facet, new_facets)
    # TODO: Do equivalence check
    new_facets = equivalence_check(new_facets, relabels)