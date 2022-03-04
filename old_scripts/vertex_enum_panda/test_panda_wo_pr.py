import argparse
from polybell.utils import get_deterministic_behaviors, get_configs, find_local_weight_dual, general_pr_box, \
    check_equiv_bell_vertex_enum_non_rescale
from polybell.panda_helper import write_panda_input_inequalities, run_panda, write_known_vertices, read_vertices_rays
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument(dest='ma', help="number of inputs for ALICE")
parser.add_argument(dest='mb', help='number of inputs for BOB')
parser.add_argument(dest='n', help='number of outputs')
parser.add_argument(dest='threads', help='number of threads to use')
args = parser.parse_args()
ma = int(args.ma)
mb = int(args.mb)
n = int(args.n)
threads = int(args.threads)

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

dets_unshifted = np.copy(dets)

relabels = np.loadtxt('../data/relabels/{}{}{}{}.gz'.format(ma, mb, n, n)).astype(int)

# Define new origin, such that the 0 origin is inside of the polytope (for dual representation)
p_origin = np.sum(dets, axis=0) / dets.shape[0]
# shift the origin of the deterministic points
dets = dets - p_origin

# get configurations
configs = get_configs(inputs_a, inputs_b, outputs, outputs)

# setup the constraints: dets @ bell <= 1
lhs = dets
rhs = np.ones(dets.shape[0])

# write the file
hrepr = write_panda_input_inequalities(lhs, rhs, symmetries=relabels, dets=dets_unshifted, file='input.ine')

# run the file
run_panda('input.ine', outfile='out.ext', threads=threads)

vertices, rays = read_vertices_rays('out.ext')

# check how many classes
classes = [vertices[0]]
for b in vertices:
    equiv = False
    for c in classes:
        if check_equiv_bell_vertex_enum_non_rescale(c, b, relabels, dets_unshifted):
            equiv = True
            continue
    if not equiv:
        classes.append(b)

classes = np.array(classes)
np.savetxt(outfile, classes)
print('stored {} classes'.format(classes.shape[0]))
