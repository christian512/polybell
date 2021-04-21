LoadPackage("json");

#Define the Symmetry Group
GRP := Group((1,13)(2,14)(3,15)(4,16)(5,9)(6,10)(7,11)(8,12),(1,4)(2,3)(5,8)(6,7)(9,12)(10,11)(13,16)(14,15),(5,9)(6,10)(7,11)(8,12),(2,3)(6,7)(10,11)(14,15));
GRP_RED := Group(MinimalGeneratingSet(GRP));

# Storage for all polytopes found
all_polys := [];

while true do
        # setup files
        outfile := IO_File("/home/chris/fromgap.pipe", "w");
        infile := IO_File("/home/chris/togap.pipe", "r");
        # read command from input
        str := IO_ReadLine(infile);
        if str <> "" then
            Print("GAP read: ", str);
            # Conver to Gap object
            arr := JsonStringToGap(str);
            level := arr[1];
            Remove(arr, 1);
            Print("Level: ", level, "\n");
            # iterate over all polytopes
            equiv := 0;
            for poly in all_polys do
                res := RepresentativeAction(GRP_RED, poly, arr, OnSets);
                if res <> fail then
                    equiv := 1;
                    IO_WriteLine(outfile, "true");
                    break;
                fi;
            od;
            if equiv = 0 then
                Add(all_polys, arr);
                Print("All_polys: ", all_polys, "\n");
                IO_WriteLine(outfile, "false");
            fi;
        fi;
od;

