""" Writing the vertices to a file and printing the generators of the automorphism group for a bell scenario """

from linearbell.utils import get_relabelling_generators, get_deterministic_behaviors_two_party, get_relabels_dets
from linearbell.gap_helper import relabels_dets_to_disjoint_cycles
from linearbell.panda_helper import write_known_vertices

# set the scenario here
ma, mb, na, nb = 2, 2, 4, 4

inputs_a, inputs_b, outputs_a, outputs_b = range(ma), range(mb), range(na), range(nb)

generators = get_relabelling_generators(inputs_a, inputs_b, outputs_a, outputs_b)
vertices = get_deterministic_behaviors_two_party(inputs_a, inputs_b, outputs_a, outputs_b)
automorphisms = get_relabels_dets(vertices, generators, show_progress=1)
disjoint_cycles = relabels_dets_to_disjoint_cycles(automorphisms)
write_known_vertices(vertices, 'test2244.ext')
print("Number of relabelling generators: ", len(generators))
print(disjoint_cycles)