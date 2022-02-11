""" Templates for the three different equivalence checks provided """

template_global_equiv_test = """
LoadPackage("json");

#Define the Symmetry Group
GRP_RED := Group(%s);

# Storage for all polytopes found
all_polys := [];
max_recursion := 30;
for i in [1..max_recursion] do
    Add(all_polys, []);
od;

# setup files
outfile := IO_File(Concatenation(IO_getcwd(), "/fromgap.pipe"), "w");
infile := IO_File(Concatenation(IO_getcwd(), "/togap.pipe"), "r");

Print("Running GAP Equivalence check server with global storage option \\n");

while true do
        # read command from input
        str := IO_ReadLine(infile);
        if str <> "" then
            # Print("GAP READ: ", str);
            # Convert to GAP Object
            arr := JsonStringToGap(str);
            
            # First is recrusion level
            level := arr[1][1] + 1;
            if level < 0 then
                break;
            fi;
            # Seond is the list of facets to test
            tarr := arr[2];

            if level > max_recursion then
                Print("WATCH OUT: MAX RECURSION LEVEL IN GAP REACHED");
            fi;

            response := [];
            for i in [1..Length(tarr)] do
                cpoly := tarr[i];
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

template_local_equiv_test = """
LoadPackage("json");

#Define the Symmetry Group
GRP_RED := Group(%s);

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
            if recursion_level < 0 then
                break;
            fi;
            # Print("Recursion Level: ", recursion_level, "\\n");
            # Second is the list of Facets to test
            tarr := arr[2];
            # Third is the list of known facets
            karr := arr[3];
            
            # Response to return to RANDA
            response := [];
            # Iterate through all arrays to test
            for i in [1..Length(tarr)] do
                tpoly := tarr[i];
                equiv := 0;
                for kpoly in karr do
                    res := RepresentativeAction(GRP_RED, tpoly, kpoly, OnSets);
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

template_stabilizer_equiv_test = """
LoadPackage("json");

#Define the Symmetry Group
GRP_RED := Group(%s);

# setup files
outfile := IO_File(Concatenation(IO_getcwd(), "/fromgap.pipe"), "w");
infile := IO_File(Concatenation(IO_getcwd(), "/togap.pipe"), "r");

str := IO_ReadLine(infile);
while str <> "break" do
        # read command from input
        if str <> "" then
            #Print("GAP READ: ", str);
            # Convert to GAP Object
            arr := JsonStringToGap(str);
            # First is the recursion level
            recursion_level := arr[1][1] + 1;
            if recursion_level < 0 then
                break;
            fi;
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
        str := IO_ReadLine(infile); 
od;
"""