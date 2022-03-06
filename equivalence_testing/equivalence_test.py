"""
This script uses GAP to check the equivalence of facets given in the PANDA file format.
"""
import argparse
from polybell.panda_helper import read_inequalities, write_known_inequalities
from polybell.utils import get_deterministic_behaviors_two_party, get_relabelling_generators, get_relabels_dets
from polybell.gap_helper import relabels_dets_to_disjoint_cycles
import numpy as np
import os

# parse the inputs
parser = argparse.ArgumentParser()
parser.add_argument(dest='ma', help="number of inputs for ALICE")
parser.add_argument(dest='mb', help='number of inputs for BOB')
parser.add_argument(dest='na', help='number of outputs for ALICE')
parser.add_argument(dest='nb', help='number of outputs for BOB')
parser.add_argument(dest='input_file', type=str,
                    help='Path to input file (relative to this script). Input needs to be given in the standard PANDA/PORTA format.')
parser.add_argument('-o', type=str,default='', dest='output_file', help='Output file (default: no output stored).')
args = parser.parse_args()

# set scenario
ma, mb, na, nb = int(args.ma), int(args.mb), int(args.na), int(args.nb)

# set filename with inequalities to check
try:
    ineqs = read_inequalities(args.input_file)
except:
    print('Error opening file -> Did you give the right path?')
print('Obtained inequalities from file')


# get vertices and symmetry generators
vertices = get_deterministic_behaviors_two_party(range(ma), range(mb), range(na), range(nb))
generators = get_relabelling_generators(range(ma), range(mb), range(na), range(nb))
automorphisms = get_relabels_dets(vertices, generators, show_progress=0)
disjoint_cycles = relabels_dets_to_disjoint_cycles(automorphisms)

# read the inequalities and generate list of vertex indices on each facet
vertex_indices_on_facets = []
for ineq in ineqs:
    indices = []
    for i in range(vertices.shape[0]):
        if vertices[i] @ ineq[:-1] == -1.0 * ineq[-1]:
            # add one since GAP deals with the vertices starting by one
            indices.append(i + 1)
    vertex_indices_on_facets.append(indices)
print('calculated the indices of vertices on each facet')

# Add vertex indices to the GAP script
gap_script = "VERTLIST := ["
for indices in vertex_indices_on_facets:
    gap_script += str(indices).replace('.', ',') + ',\n'
gap_script = gap_script[:-2] + ']; \n'

# add symmetry generators to gap script
gap_script += "GRP := Group(" + disjoint_cycles + "); \n"

# Write functions to compute all inequivalent facets
gap_script += """
INEQUIV_FACETS := [];
INDICES_INEQUIV_FACETS := [];
for i in [1..Length(VERTLIST)] do
    facet := VERTLIST[i];
    equiv := 0;
    for test_facet in INEQUIV_FACETS do
        res := RepresentativeAction(GRP, facet, test_facet, OnSets);
        if res <> fail then
            equiv := 1;
            break;
        fi;
    od;
    if equiv = 0 then
        Add(INEQUIV_FACETS, facet);
        Add(INDICES_INEQUIV_FACETS, i-1);
    fi;
od;
outfile := IO_File(Concatenation(IO_getcwd(), "/indices_inequiv_facets.txt"), "w");
IO_WriteLine(outfile, INDICES_INEQUIV_FACETS);
IO_Close(outfile);"""
f = open('equivalence_script.g', 'w+')
f.write(gap_script)
f.close()

# Execute Gap script
os.system("gap -q --nointeract equivalence_script.g")
os.system("rm equivalence_script.g")

# read output from gap script
f = open('indices_inequiv_facets.txt')
lines = f.read().replace(' ', '').replace('[', '').replace(']', '')
f.close()
indices_inequiv_facets = np.fromstring(lines, sep=',', dtype=int)

print('Found {} inequivalent facets'.format(indices_inequiv_facets.shape[0]))

# write inequivalent facets to file
if args.output_file:
    inequiv_facets = ineqs[indices_inequiv_facets]
    lhs = inequiv_facets[:, :-1]
    rhs = -1.0 * inequiv_facets[:, -1]
    write_known_inequalities(lhs, rhs, args.output_file)
    print('wrote inequivalent facets to {}'.format(args.output_file))
os.system("rm indices_inequiv_facets.txt")