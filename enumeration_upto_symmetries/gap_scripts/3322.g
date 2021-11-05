LoadPackage("json");

#Define the Symmetry Group
GRP_RED := Group((17,33)(18,34)(19,35)(20,36)(21,37)(22,38)(23,39)(24,40)(25,41)(26,42)(27,43)(28,44)(29,45)(30,46)(31,47)(32,48),
(9,33)(10,34)(11,35)(12,36)(13,37)(14,38)(15,39)(16,40)(25,49)(26,50)(27,51)(28,52)(29,53)(30,54)(31,55)(32,56),
(3,5)(4,6)(11,13)(12,14)(19,21)(20,22)(27,29)(28,30)(35,37)(36,38)(43,45)(44,46)(51,53)(52,54)(59,61)(60,62),
(2,5)(4,7)(10,13)(12,15)(18,21)(20,23)(26,29)(28,31)(34,37)(36,39)(42,45)(44,47)(50,53)(52,55)(58,61)(60,63),
(1,33)(2,34)(3,35)(4,36)(5,37)(6,38)(7,39)(8,40)(9,41)(10,42)(11,43)(12,44)(13,45)(14,46)(15,47)(16,48)(17,49)(18,50)(19,51)(20,52)(21,53)(22,54)(23,55)(24,56)(25,57)(26,58)(27,59)(28,60)(29,61)(30,62)(31,63)(32,64),
(1,5)(2,6)(3,7)(4,8)(9,13)(10,14)(11,15)(12,16)(17,21)(18,22)(19,23)(20,24)(25,29)(26,30)(27,31)(28,32)(33,37)(34,38)(35,39)(36,40)(41,45)(42,46)(43,47)(44,48)(49,53)(50,54)(51,55)(52,56)(57,61)(58,62)(59,63)(60,64),
(2,9)(3,17)(4,25)(5,33)(6,41)(7,49)(8,57)(11,18)(12,26)(13,34)(14,42)(15,50)(16,58)(20,27)(21,35)(22,43)(23,51)(24,59)(29,36)(30,44)(31,52)(32,60)(38,45)(39,53)(40,61)(47,54)(48,62)(56,63));
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
                    Print("Level ", level, ": ",  Length(all_polys[level]), "
");
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
