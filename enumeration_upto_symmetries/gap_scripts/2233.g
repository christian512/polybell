LoadPackage("json");

#Define the Symmetry Group
GRP_RED := Group((10,28)(11,29)(12,30)(13,31)(14,32)(15,33)(16,34)(17,35)(18,36)(19,55)(20,56)(21,57)(22,58)(23,59)(24,60)(25,61)(26,62)(27,63)(46,64)(47,65)(48,66)(49,67)(50,68)(51,69)(52,70)(53,71)(54,72),
(2,4)(3,7)(6,8)(11,13)(12,16)(15,17)(20,22)(21,25)(24,26)(29,31)(30,34)(33,35)(38,40)(39,43)(42,44)(47,49)(48,52)(51,53)(56,58)(57,61)(60,62)(65,67)(66,70)(69,71)(74,76)(75,79)(78,80),
(1,28)(2,29)(3,30)(4,31)(5,32)(6,33)(7,34)(8,35)(9,36)(10,37)(11,38)(12,39)(13,40)(14,41)(15,42)(16,43)(17,44)(18,45)(19,46)(20,47)(21,48)(22,49)(23,50)(24,51)(25,52)(26,53)(27,54),
(1,55)(2,56)(3,57)(4,58)(5,59)(6,60)(7,61)(8,62)(9,63)(10,64)(11,65)(12,66)(13,67)(14,68)(15,69)(16,70)(17,71)(18,72)(19,73)(20,74)(21,75)(22,76)(23,77)(24,78)(25,79)(26,80)(27,81),
(1,4)(2,5)(3,6)(10,13)(11,14)(12,15)(19,22)(20,23)(21,24)(28,31)(29,32)(30,33)(37,40)(38,41)(39,42)(46,49)(47,50)(48,51)(55,58)(56,59)(57,60)(64,67)(65,68)(66,69)(73,76)(74,77)(75,78),
(1,7)(2,8)(3,9)(10,16)(11,17)(12,18)(19,25)(20,26)(21,27)(28,34)(29,35)(30,36)(37,43)(38,44)(39,45)(46,52)(47,53)(48,54)(55,61)(56,62)(57,63)(64,70)(65,71)(66,72)(73,79)(74,80)(75,81),
(2,10)(3,19)(4,28)(5,37)(6,46)(7,55)(8,64)(9,73)(12,20)(13,29)(14,38)(15,47)(16,56)(17,65)(18,74)(22,30)(23,39)(24,48)(25,57)(26,66)(27,75)(32,40)(33,49)(34,58)(35,67)(36,76)(42,50)(43,59)(44,68)(45,77)(52,60)(53,69)(54,78)(62,70)(63,79)(72,80));
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
                    Print("Level ", level, ": ",  Length(all_polys[level]), "\n");
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
