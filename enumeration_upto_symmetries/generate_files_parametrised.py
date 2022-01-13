""" Writing the vertices to a file and printing the generators of the automorphism group for a bell scenario """

from linearbell.utils import get_relabelling_generators, get_deterministic_behaviors_two_party, get_relabels_dets, \
    get_configs, get_parametrisation_configs, parametrise_behavior
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

# Write Vertices to a File
vertices = get_deterministic_behaviors_two_party(inputs_a, inputs_b, outputs_a, outputs_b)
vertices_param = []
for v in vertices:
    v_param = parametrise_behavior(v, configs, configs_param, inputs_a, inputs_b, outputs_a, outputs_b)
    vertices_param.append(v_param)
vertices_param = np.array(vertices_param)
write_known_vertices(vertices_param, 'randa_files/{}{}{}{}.ext'.format(ma, mb, na, nb))

# Calculate the generators
generators = get_relabelling_generators(inputs_a, inputs_b, outputs_a, outputs_b)
automorphisms = get_relabels_dets(vertices, generators, show_progress=1)
disjoint_cycles = relabels_dets_to_disjoint_cycles(automorphisms)

# write gap script
gap_script_start = """LoadPackage("json");

#Define the Symmetry Group
GRP_RED := Group("""

gap_script_end = """);
# Storage for all polytopes found
all_polys := [];
max_recursion := 30;
for i in [1..max_recursion] do
    Add(all_polys, []);
od;

# setup files

outfile := IO_File(Concatenation(IO_getcwd(), "/fromgap.pipe"), "w");
infile := IO_File(Concatenation(IO_getcwd(), "/togap.pipe"), "r");

while true do
        # read command from input
        str := IO_ReadLine(infile);
        if str <> "" then
            # Print("GAP READ: ", str);
            # Convert to GAP Object
            arr := JsonStringToGap(str);
            # Extract Level
            level := arr[1][1] + 1;

            Remove(arr,1);
            if level > max_recursion then
                Print("WATCH OUT: MAX RECURSION LEVEL IN GAP REACHED");
            fi;

            response := [];
            for i in [1..Length(arr)] do
                cpoly := arr[i];
                equiv := 0;
                for tpoly in all_polys[level] do
                    res := RepresentativeAction(GRP_RED, tpoly, cpoly, OnSets);
                    if res <> fail then
                        equiv := 1;
                        break;
                    fi;
                od;
                if equiv = 0 then
                    Add(all_polys[level], cpoly);
                    Add(response, i-1);
                    Print("Level ", level, ": ",  Length(all_polys[level]), "\\n");
                fi;
            od;
            # write false in case there is no equivalent polytope
            if Length(response) = 0 then
                IO_WriteLine(outfile, "false");
            else
                IO_WriteLine(outfile, response);
            fi;
        fi;
od;
"""

gap_script = gap_script_start + disjoint_cycles + gap_script_end
f = open('gap_scripts/{}{}{}{}_parametrised.g'.format(ma, mb, na, nb), 'w+')
f.write(gap_script)
f.close()
