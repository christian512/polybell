""" Writing the vertices to a file and printing the generators of the automorphism group for a bell scenario """

from linearbell.utils import get_relabelling_generators, get_deterministic_behaviors_two_party, get_relabels_dets
from linearbell.gap_helper import relabels_dets_to_disjoint_cycles
from linearbell.panda_helper import write_known_vertices
import argparse

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

# Write Vertices to a File
vertices = get_deterministic_behaviors_two_party(inputs_a, inputs_b, outputs_a, outputs_b)
# Calculate the generators
generators = get_relabelling_generators(inputs_a, inputs_b, outputs_a, outputs_b)
automorphisms = get_relabels_dets(vertices, generators, show_progress=1)
disjoint_cycles = relabels_dets_to_disjoint_cycles(automorphisms)

# Vertices
gap_script = "EXT := ["
for v in vertices:
    gap_script += str(v).replace('.', ',') + ',\n'
gap_script = gap_script[:-2] + "]; \n"
# Generators
gap_script += "GRP := Group(" + disjoint_cycles + "); \n"
# Function Execution
gap_script += """MaxRecursionLvl:= 7;

eRec := GetIAI_FromEXT_GRP(EXT, GRP, 500);
FuncStabilizer:=eRec.FuncStabilizer;
FuncIsomorphy:=eRec.FuncIsomorphy;
FuncInvariant:=eRec.FuncInvariant;

IsRespawn:=function(OrdGRP, EXT, TheDepth)
  if TheDepth > MaxRecursionLvl then
    return false;
  fi;
  if Length(EXT) < 100 then
    return false;
  fi;
  return true;
end;

IsBankSave:=function(EllapsedTime, OrdGRP, EXT, TheDepth)
	return false;
end;

ThePathWork:="./TheWork/";
Exec("mkdir -p ", ThePathWork);
ThePath:=Concatenation(ThePathWork, "tmp/");
Exec("mkdir -p ", ThePath);
PathSave:=Concatenation(ThePathWork, "PathSAVE/");
Exec("mkdir -p ", PathSave);

TheFunc:=__DualDescriptionCDD_Reduction;

GetInitialRays_SamplingFramework:=function(EXT,nb)
  return GetInitialRays_LinProg(EXT,nb);
end;

TestNeedMoreSymmetry:=function(EXT)
	return false;
end;


Data:=rec(TheDepth:=0,
          ThePath:=ThePath,
          GetInitialRays:=GetInitialRays_SamplingFramework,
          IsBankSave:=IsBankSave,
          GroupFormalism:=OnSetsGroupFormalism(500),
          DualDescriptionFunction:=TheFunc,
          TestNeedMoreSymmetry:=TestNeedMoreSymmetry,
          IsRespawn:=IsRespawn,
          Saving:=true,
          ThePathSave:=PathSave);

BankPath:=Concatenation(ThePathWork, "TheBank/");
Exec("mkdir -p ", BankPath);
DataBank:=rec(BankPath:=BankPath, Saving:=false);
BF:=BankRecording(DataBank, FuncStabilizer, FuncIsomorphy, FuncInvariant, OnSetsGroupFormalism(500));


cur_time := NanosecondsSinceEpoch();
LORB:=__ListFacetByAdjacencyDecompositionMethod(EXT,
                                                GRP,
                                                Data,
                                                BF);

calc_time := Float((NanosecondsSinceEpoch() - cur_time)*10^(-6));
SaveDataToFile("ListOrbitEXT", LORB);
"""


f = open('scripts/{}{}{}{}_polyhedral_adm_dd.g'.format(ma, mb, na, nb), 'w+')
f.write(gap_script)
f.close()
