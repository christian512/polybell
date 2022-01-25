""" Writing the vertices to a file and printing the generators of the automorphism group for a bell scenario """

from linearbell.utils import get_relabelling_generators, get_deterministic_behaviors_two_party, get_relabels_dets
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
ma, mb, na, nb = 5,3,2,2
inputs_a, inputs_b, outputs_a, outputs_b = range(ma), range(mb), range(na), range(nb)

# Write Vertices to a File
vertices = get_deterministic_behaviors_two_party(inputs_a, inputs_b, outputs_a, outputs_b)
vertices = vertices[np.lexsort(np.rot90(vertices))]
write_known_vertices(vertices, 'randa_files/{}{}{}{}.ext'.format(ma, mb, na, nb))

# Calculate the generators
generators = get_relabelling_generators(inputs_a, inputs_b, outputs_a, outputs_b)
automorphisms = get_relabels_dets(vertices, generators, show_progress=1)
disjoint_cycles = relabels_dets_to_disjoint_cycles(automorphisms)

# write gap script
gap_script_start = """LoadPackage("json");

#Define the Symmetry Group
GRP_RED := Group("""

gap_script_end = """);
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
            # First is the recursion level
            recursion_level := arr[1][1] + 1;
            # Print("Recursion Level: ", recursion_level, "\\n");
            # Second is the list of Facets to test
            tarr := arr[2];
            # Third is the list of known facets
            karr := arr[3];
            # Fourth is the list of vertices on the polytope
            vertices := arr[4];
            # Caluclate the stabilizer
            GRP := GRP_RED;
            if recursion_level > 1 then
                GRP := Stabilizer(GRP_RED, vertices, OnSets);
            fi;
            
            # Response to return to RANDA
            response := [];
            # Iterate through all arrays to test
            for i in [1..Length(tarr)] do
                tpoly := tarr[i];
                equiv := 0;
                for kpoly in karr do
                    res := RepresentativeAction(GRP, tpoly, kpoly, OnSets);
                    if res <> fail then
                        equiv := 1;
                        break;
                    fi;
                od;
                if equiv = 0 then
                    Add(karr, tpoly);
                    Add(response, i-1);
                fi;
            od;
            
            # Respond
            if Length(response) = 0 then
                IO_WriteLine(outfile, "false");
            else
                IO_WriteLine(outfile, response);
            fi;
        fi;
od;

"""

gap_script = gap_script_start + disjoint_cycles + gap_script_end
f = open('gap_scripts/{}{}{}{}_local_stabilizer.g'.format(ma, mb, na, nb), 'w+')
f.write(gap_script)
f.close()
