from polybell.utils import get_deterministic_behaviors, get_relabelling_generators, get_relabels_dets, get_allowed_relabellings
from polybell.panda_helper import write_known_vertices
from polybell.gap_helper import relabels_dets_to_disjoint_cycles
import numpy as np

vertices = get_deterministic_behaviors(range(4), range(4), range(2))
write_known_vertices(vertices, 'test.ext')
