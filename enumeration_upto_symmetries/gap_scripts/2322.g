LoadPackage("json");

#Define the Symmetry Group
GRP_RED := Group((5,9)(6,10)(7,11)(8,12),
(2,3)(6,7)(10,11)(14,15),
(1,9)(2,10)(3,11)(4,12)(5,13)(6,14)(7,15)(8,16),
(1,3)(2,4)(5,7)(6,8)(9,11)(10,12)(13,15)(14,16),
(2,5)(3,9)(4,13)(7,10)(8,14)(12,15));
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
