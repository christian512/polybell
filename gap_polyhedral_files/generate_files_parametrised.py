""" Writing the vertices to a file and printing the generators of the automorphism group for a bell scenario """

from linearbell.utils import get_relabelling_generators, get_deterministic_behaviors_two_party, get_relabels_dets, get_configs, get_parametrisation_configs, parametrise_behavior
from linearbell.gap_helper import relabels_dets_to_disjoint_cycles
from linearbell.panda_helper import write_known_vertices
import argparse
import numpy as np

# parse the inputs
parser = argparse.ArgumentParser()
parser.add_argument(dest='ma', help="number of inputs for ALICE")
parser.add_argument(dest='mb', help='number of inputs for BOB')
parser.add_argument(dest='na', help='number of outputs for ALICE')
parser.add_argument(dest='nb', help='number of outputs for BOB')
args = parser.parse_args()

# set the scenario
ma, mb, na, nb = int(args.ma), int(args.mb), int(args.na), int(args.nb)
inputs_a, inputs_b, outputs_a, outputs_b = range(ma), range(mb), range(na), range(nb)

# setup configurations for extended and parametrised behavior
configs = get_configs(inputs_a, inputs_b, outputs_a, outputs_b)
configs_param = get_parametrisation_configs(inputs_a, inputs_b, outputs_a, outputs_b)

# Generate parametrised vertices
vertices = get_deterministic_behaviors_two_party(inputs_a, inputs_b, outputs_a, outputs_b)

# Calculate the generators
generators = get_relabelling_generators(inputs_a, inputs_b, outputs_a, outputs_b)
automorphisms = get_relabels_dets(vertices, generators, show_progress=1)
disjoint_cycles = relabels_dets_to_disjoint_cycles(automorphisms)

# parametrise vertices
vertices_param = []
for v in vertices:
    v_param = parametrise_behavior(v, configs, configs_param, inputs_a, inputs_b, outputs_a, outputs_b)
    vertices_param.append(v_param)
vertices = np.array(vertices_param)

# Vertices
gap_script = "EXT := ["
for v in vertices:
    gap_script += str(v).replace('.', ',')[:-2] + '],\n'
gap_script = gap_script[:-2] + "]; \n"
# Generators
gap_script += "GRP := Group(" + disjoint_cycles + "); \n"
# Function Execution
gap_script += "cur_time := NanosecondsSinceEpoch();\n"
gap_script += "DualDescriptionStandard(EXT,GRP);\n"
gap_script += "calc_time := NanosecondsSinceEpoch() - cur_time;\n"
gap_script += "Print(calc_time * 10^(-9));"

f = open('scripts/{}{}{}{}_polyhedral_parametrised.g'.format(ma, mb, na, nb), 'w+')
f.write(gap_script)
f.close()
