import argparse
import numpy as np
import subprocess
from polybell.utils import get_deterministic_behaviors, check_equiv_bell_vertex_enum_non_rescale
from polybell.panda_helper import write_known_vertices, read_inequalities

parser = argparse.ArgumentParser()
parser.add_argument(dest='ma', help="number of inputs for ALICE")
parser.add_argument(dest='mb', help='number of inputs for BOB')
parser.add_argument(dest='n', help='number of outputs')
args = parser.parse_args()
ma = int(args.ma)
mb = int(args.mb)
n = int(args.n)

# set inputs / outputs
inputs_a = range(ma)
inputs_b = range(mb)
outputs = range(n)
outputs_wo_failure = range(n - 1)

# setup output file
outfile = '../data/vertex_enum/{}{}{}{}.ext'.format(ma, mb, n, n)

# get deterministic points
dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs)
# sort deterministics
dets = dets[np.lexsort(np.rot90(dets))]
relabels = np.loadtxt('../data/relabels/{}{}{}{}.gz'.format(ma, mb, n, n)).astype(int)

test_vertices = np.array(
    [[1, 1, 1], [-1, 1, 1], [1, -1, 1], [-1, -1, 1], [1, 1, -1], [-1, 1, -1], [1, -1, -1], [-1, -1, -1]])


# write known vertices
write_known_vertices(test_vertices, file='knownvertices.ext')

# run panda
cmd = 'panda_org knownvertices.ext -t 1 > out.ine'
out = subprocess.run(cmd, shell=True)
print('done processing in PANDA')

# read inequalities
ineq = read_inequalities('out.ine')

# check for equivalence
facets = ineq[:, :-1]
print('number facets: ', facets.shape[0])

# classes = [facets[0]]
#
# for i, f in enumerate(facets):
#     print('{} / {}'.format(i, facets.shape[0]))
#     equiv = False
#     for c in classes:
#         if check_equiv_bell_vertex_enum_non_rescale(f, c, relabels, dets):
#             equiv = True
#             break
#     if not equiv:
#         classes.append(f)
# classes = np.array(classes)
# print('num classes ', classes.shape[0])
