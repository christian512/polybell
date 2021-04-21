from linearbell.utils import get_relabelling_generators, get_deterministic_behaviors, get_relabels_dets
from linearbell.gap_helper import relabels_dets_to_disjoint_cycles
ma, mb, na, nb = 4, 4, 2, 2

inputs_a, inputs_b, outputs_a, outputs_b = range(ma), range(mb), range(na), range(nb)

generators = get_relabelling_generators(inputs_a, inputs_b, outputs_a, outputs_b)
vertices = get_deterministic_behaviors(inputs_a, inputs_b, outputs_a)
automorphisms = get_relabels_dets(vertices, generators, show_progress=1)
disjoint_cycles = relabels_dets_to_disjoint_cycles(automorphisms)
print(disjoint_cycles)