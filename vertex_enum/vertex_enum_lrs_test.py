import argparse
from linearbell.utils import get_deterministic_behaviors, get_configs, general_pr_box_extended, \
    get_parametrisation_configs, parametrise_behavior, get_allowed_relabellings, check_equiv_bell
import numpy as np
from lrs_helper import polyhedra_h_representation, run_lrs_polytope

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

# get deterministic points
dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs)

# get origin
p_origin = np.sum(dets, axis=0) / (dets.shape[0])
# rescale the deterministic points
dets = dets - p_origin

# get allowed relabellings
allowed_relabellings = get_allowed_relabellings(inputs_a, inputs_b, outputs, outputs)
# get relabellings for deterministic points
file_relabels = '../data/relabels_dets/{}{}{}{}.gz'.format(ma, mb, n, n)
relabels_dets = np.loadtxt(file_relabels, dtype=float).astype(int)

# get the configurations
configs = get_configs(inputs_a, inputs_b, outputs, outputs)

# Parametrise to get lower dimensional representations
configs_param = get_parametrisation_configs(inputs_a, inputs_b, outputs, outputs)
# get parametrisations of deterministic points
dets_param = []
for d in dets:
    d_param = parametrise_behavior(d, configs, configs_param, inputs_a, inputs_b, outputs, outputs)
    dets_param.append(d_param)
dets_param = np.array(dets_param)
#dets = np.copy(dets_param)

# setup the constraints

# dets @ b <= 1
lhs_ineq = np.copy(dets)
rhs_ineq = np.ones(dets.shape[0])

# -1.0 * b <= 0 equiv to b >= 0
lhs_ineq = np.r_[lhs_ineq, -1.0 * np.eye(dets.shape[1])]
rhs_ineq = np.r_[rhs_ineq, np.zeros(dets.shape[1])]

# compute vertices of the polytope
print('Start computation of vertices')
polyhedra_h_representation(lhs_ineq, rhs_ineq, file='lrs.ine')
vertices = run_lrs_polytope('lrs.ine', 'out.ext')

#vertices = compute_polytope_vertices(lhs_ineq, rhs_ineq)

# check that vertices are facets
vertices = np.array(vertices)
# check if all points are facets
facets = []
# check facets
facets = []
for v in vertices:
    eq_dets = []
    for d in dets:
        if np.abs(d @ v - 1) < 1e-2:
            eq_dets.append(d)
    eq_dets = np.array(eq_dets)
    if eq_dets.shape[0] > 0:
        eq_dets_new = eq_dets - eq_dets[0]
        # calculate the rank of the matrix
        rank = np.linalg.matrix_rank(eq_dets_new)
        # check if it's a facet by rank check
        if rank == (ma * (n - 1) + 1) * (mb * (n - 1) + 1) - 2:
            facets.append(v)
facets = np.array(facets)

# store the vertices
file = '../data/vertex_enum/{}{}{}{}.txt'.format(ma, mb, n, n)
np.savetxt(file, vertices)
print('done')
print('number of vertices: {}'.format(vertices.shape[0]))
print('number of facets: {}'.format(facets.shape[0]))

del_facets = []
# iterate through the facets
for i in range(facets.shape[0]):
    # if this facet can already be deleted -> continue
    # print('facet: {} / {} || len deletion list: {}'.format(i, facets.shape[0], len(del_facets)))
    if i in del_facets: continue
    for j in range(i + 1, facets.shape[0]):
        # if this facet can already be deleted -> continue
        if j in del_facets: continue
        # check if the two facets are equivalent
        if check_equiv_bell(facets[i], facets[j], relabels_dets, dets, tol=1e-4):
            del_facets.append(j)
# store new facets
new_facets = np.delete(facets, del_facets, axis=0)
print('Number of unique facets: {}'.format(new_facets.shape[0]))