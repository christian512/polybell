from linearbell.utils import *
from linearbell.lrs_helper import read_v_file

# setup parameters
ma = 2
mb = 2
n = 3

# file where all the facets stored
file = '../data/vertex_enum_pr_box/2233_local_weight_lins.ext'

# set inputs / outputs
inputs_a = range(ma)
inputs_b = range(mb)
outputs = range(n)
outputs_wo_failure = range(n - 1)

# setup configurations and deterministic behaviors
configs = get_configs(inputs_a, inputs_b, outputs, outputs)
configs_param = get_parametrisation_configs(inputs_a, inputs_b, outputs, outputs)
# get deterministics
dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs)
# define the origin
p_origin = np.sum(dets, axis=0) / dets.shape[0]

# get allowed relabellings
allowed_relabellings = get_allowed_relabellings(inputs_a, inputs_b, outputs, outputs)
# get relabellings for deterministic points
file_relabels = '../data/relabels_dets/{}{}{}{}.gz'.format(ma, mb, n, n)
try:
    relabels_dets = np.loadtxt(file_relabels, dtype=float).astype(int)
except IOError:
    print('Have to calculate the possible relabels before actual start')
    relabels_dets = get_relabels_dets(dets, allowed_relabellings)
    np.savetxt(file_relabels, relabels_dets)

# load bell expressions from file
vertices, rays = read_v_file(file)

# define PR box
eta = 4 / (ma + 4)
pr = np.array(
    [general_pr_box_extended(a, b, x, y, eta, outputs_wo_failure) for a, b, x, y in configs])
pr = pr - p_origin
pr_param = parametrise_behavior(pr, configs, configs_param, inputs_a, inputs_b, outputs, outputs)

# selected vertices
sel_vertices = []
# check which vertices fulfill the condition that PR_eta * B == 1 for eta = 4 / (m+4)
for v_param in vertices:
    bell = deparametrise_bell_expression(v_param, configs, configs_param, inputs_a, inputs_b, outputs, outputs)
    if bell @ pr == 1.0:
        sel_vertices.append(v_param)
print('number of vertices v with v @ pr_box == 1: {}'.format(len(sel_vertices)))

# rescale the bell expressions
bell_expressions = []
sum_ones = np.sum(dets, axis=1)[0]
for v_param in sel_vertices:
    bell = deparametrise_bell_expression(v_param, configs, configs_param, inputs_a, inputs_b, outputs, outputs)
    v = np.array([bell @ d for d in dets])
    min_val = np.min(v)
    v = v + np.abs(min_val)
    sec_min_pos_val = np.min(v[v > 1e-4])
    # shift the bell expression
    bell = (bell + np.abs(min_val) / sum_ones) / sec_min_pos_val
    bell_expressions.append(bell)
    # new v vector
    v = np.array([bell @ d for d in dets])
facets = np.array(bell_expressions)

# find the number of classes
del_facets = []
for i in range(facets.shape[0]):
    # if this facet can already be deleted -> continue
    if i in del_facets: continue
    for j in range(i + 1, facets.shape[0]):
        # if this facet can already be deleted -> continue
        if j in del_facets: continue
        # check if the two facets are equivalent
        if check_equiv_bell(facets[i], facets[j], relabels_dets, dets, tol=1e-4):
            del_facets.append(j)
# store new facets
new_facets = np.delete(facets, del_facets, axis=0)

# print the results
print('Number of classes: {}'.format(new_facets.shape[0]))
print(new_facets)
