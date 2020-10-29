from linearbell.utils import get_deterministic_behaviors, get_allowed_relabellings, check_equiv_bell
import numpy as np

# set inputs / outputs
inputs_a = range(4)
inputs_b = range(4)
outputs = range(2)

# get deterministic points
dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs)

# get allowed relabellings
allowed_relabellings = get_allowed_relabellings(inputs_a, inputs_b, outputs, outputs)

# load facets from file
facets_file = '../data/facets/4422.txt'
facets = np.loadtxt(facets_file)

# indices of facets to delete
del_facets = []
# iterate through the facets
for i in range(facets.shape[0]):
    # if this facet can already be deleted -> continue
    if i in del_facets: continue
    for j in range(i+1, facets.shape[0]):
        # if this facet can already be deleted -> continue
        if j in del_facets: continue
        # check if the two facets are equivalent
        if check_equiv_bell(facets[i], facets[j], allowed_relabellings, dets):
            del_facets.append(j)
# store new facets
new_facets = np.delete(facets, del_facets, axis=0)
np.savetxt(facets_file, new_facets)
