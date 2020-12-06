from linearbell.utils import get_deterministic_behaviors, find_local_weight_primal, get_configs, general_pr_box_extended
from linearbell.lrs_helper import polyhedra_h_representation, run_lrs_h_repr, run_redund
import numpy as np

# number of inputs / outputs
ma = 2
mb = 2
n = 3

# set inputs / outputs
inputs_a = range(ma)
inputs_b = range(mb)
outputs = range(n)
outputs_wo_failure = range(n - 1)

# get deterministic points
dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs)

# setup pr box
configs = get_configs(inputs_a, inputs_b, outputs, outputs)
eta = 4 / (4 + len(inputs_a))
pr_box = [general_pr_box_extended(a, b, x, y, eta, outputs_wo_failure) for (a, b, x, y) in configs]
pr_box = np.array(pr_box)


# load the possible combinations to build the PR box from file
file = '../data/pr_box_contrib_det/{}{}{}{}.gz'.format(ma, mb, n, n)
combinations = np.loadtxt(file)
if len(combinations.shape) == 1:
    # reshape to 2D array
    combinations = combinations.reshape(1, combinations.shape[0])


# check that the combination of deterministics rebuilds the pr box
all_vertices = []
for i in range(combinations.shape[0]):
    c = combinations[i]
    pr_test = np.sum(dets * c[:, np.newaxis], axis=0)
    assert np.allclose(pr_test, pr_box)

    # write file for running the vertex enumeration
    lhs = np.copy(dets)
    lhs = np.r_[lhs, [pr_box]]
    rhs = np.ones(lhs.shape[0])
    lins = list(np.arange(c.shape[0])[c > 1e-6])
    lins.append(lhs.shape[0] - 1)

    repr = polyhedra_h_representation(lhs, rhs, lins, file='input.ine')

    # run redund
    run_redund('input.ine', 'input_redund.ine')

    # run the enumeration
    outfile = 'out.ext'
    vertices, rays = run_lrs_h_repr('input_redund.ine', output_file=outfile, nproc=1)
    for v in vertices:
        all_vertices.append(v)

# store all vertices (bell inequalities we found)
all_vertices = np.array(all_vertices)
file = '../data/vertex_enum_pr_box_det_decomp/{}{}{}{}.gz'.format(ma, mb, n, n)
np.savetxt(file, all_vertices)
