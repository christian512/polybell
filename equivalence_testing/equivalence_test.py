"""
This script uses GAP to check the equivalence of facets given in the PANDA file format.
"""

from linearbell.panda_helper import read_inequalities, write_known_inequalities
from linearbell.utils import get_deterministic_behaviors_two_party, get_relabelling_generators, get_relabels_dets
from linearbell.gap_helper import relabels_dets_to_disjoint_cycles
import numpy as np
import os

# set scenario
ma,mb,na,nb = 3,3,3,3

# set filename with inequalities to check
filename = '../facet_classes/3333_partial_tom.ine'
ineqs = read_inequalities(filename)
print('read inequalities from file')


# +++ PROGRAM BEGINS HERE +++

# get vertices and symmetry generators
vertices = get_deterministic_behaviors_two_party(range(ma),range(mb),range(na),range(nb))
generators = get_relabelling_generators(range(ma),range(mb),range(na),range(nb))
automorphisms = get_relabels_dets(vertices, generators, show_progress=0)
disjoint_cycles = relabels_dets_to_disjoint_cycles(automorphisms)

# read the inequalities and generate list of vertex indices on each facet
vertex_indices_on_facets = []
for ineq in ineqs:
    indices = []
    for i in range(vertices.shape[0]):
        if vertices[i] @ ineq[:-1] == -1.0 * ineq[-1]:
            # add one since GAP deals with the vertices starting by one
            indices.append(i+1)
    vertex_indices_on_facets.append(indices)
print('calculated the indices of vertices on each facet')

# Add vertex indices to the GAP script
gap_script = "VERTLIST := ["
for indices in vertex_indices_on_facets:
    gap_script += str(indices).replace('.',',') + ',\n'
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
lines = f.read().replace(' ','').replace('[','').replace(']','')
f.close()
indices_inequiv_facets = np.fromstring(lines,sep=',',dtype=int)

print('Found {} inequivalent facets'.format(indices_inequiv_facets.shape[0]))

# write inequivalent facets to file
inequiv_facets = ineqs[indices_inequiv_facets]
lhs = inequiv_facets[:,:-1]
rhs = -1.0 * inequiv_facets[:,-1]
fname = '{}{}{}{}_inequiv.ine'.format(ma,mb,na,nb)
write_known_inequalities(lhs, rhs, fname)
print('wrote inequivalent facets to {}'.format(fname))




