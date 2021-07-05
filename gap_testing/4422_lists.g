LoadPackage("json");

#Define the Symmetry Group
GRP_RED := Group((65,129)(66,130)(67,131)(68,132)(69,133)(70,134)(71,135)(72,136)(73,137)(74,138)(75,139)(76,140)(77,141)(78,142)(79,143)(80,144)(81,145)(82,146)(83,147)(84,148)(85,149)(86,150)(87,151)(88,152)(89,153)(90,154)(91,155)(92,156)(93,157)(94,158)(95,159)(96,160)(97,161)(98,162)(99,163)(100,164)(101,165)(102,166)(103,167)(104,168)(105,169)(106,170)(107,171)(108,172)(109,173)(110,174)(111,175)(112,176)(113,177)(114,178)(115,179)(116,180)(117,181)(118,182)(119,183)(120,184)(121,185)(122,186)(123,187)(124,188)(125,189)(126,190)(127,191)(128,192),
(33,129)(34,130)(35,131)(36,132)(37,133)(38,134)(39,135)(40,136)(41,137)(42,138)(43,139)(44,140)(45,141)(46,142)(47,143)(48,144)(49,145)(50,146)(51,147)(52,148)(53,149)(54,150)(55,151)(56,152)(57,153)(58,154)(59,155)(60,156)(61,157)(62,158)(63,159)(64,160)(97,193)(98,194)(99,195)(100,196)(101,197)(102,198)(103,199)(104,200)(105,201)(106,202)(107,203)(108,204)(109,205)(110,206)(111,207)(112,208)(113,209)(114,210)(115,211)(116,212)(117,213)(118,214)(119,215)(120,216)(121,217)(122,218)(123,219)(124,220)(125,221)(126,222)(127,223)(128,224),
(17,129)(18,130)(19,131)(20,132)(21,133)(22,134)(23,135)(24,136)(25,137)(26,138)(27,139)(28,140)(29,141)(30,142)(31,143)(32,144)(49,161)(50,162)(51,163)(52,164)(53,165)(54,166)(55,167)(56,168)(57,169)(58,170)(59,171)(60,172)(61,173)(62,174)(63,175)(64,176)(81,193)(82,194)(83,195)(84,196)(85,197)(86,198)(87,199)(88,200)(89,201)(90,202)(91,203)(92,204)(93,205)(94,206)(95,207)(96,208)(113,225)(114,226)(115,227)(116,228)(117,229)(118,230)(119,231)(120,232)(121,233)(122,234)(123,235)(124,236)(125,237)(126,238)(127,239)(128,240),
(5,9)(6,10)(7,11)(8,12)(21,25)(22,26)(23,27)(24,28)(37,41)(38,42)(39,43)(40,44)(53,57)(54,58)(55,59)(56,60)(69,73)(70,74)(71,75)(72,76)(85,89)(86,90)(87,91)(88,92)(101,105)(102,106)(103,107)(104,108)(117,121)(118,122)(119,123)(120,124)(133,137)(134,138)(135,139)(136,140)(149,153)(150,154)(151,155)(152,156)(165,169)(166,170)(167,171)(168,172)(181,185)(182,186)(183,187)(184,188)(197,201)(198,202)(199,203)(200,204)(213,217)(214,218)(215,219)(216,220)(229,233)(230,234)(231,235)(232,236)(245,249)(246,250)(247,251)(248,252),
(3,9)(4,10)(7,13)(8,14)(19,25)(20,26)(23,29)(24,30)(35,41)(36,42)(39,45)(40,46)(51,57)(52,58)(55,61)(56,62)(67,73)(68,74)(71,77)(72,78)(83,89)(84,90)(87,93)(88,94)(99,105)(100,106)(103,109)(104,110)(115,121)(116,122)(119,125)(120,126)(131,137)(132,138)(135,141)(136,142)(147,153)(148,154)(151,157)(152,158)(163,169)(164,170)(167,173)(168,174)(179,185)(180,186)(183,189)(184,190)(195,201)(196,202)(199,205)(200,206)(211,217)(212,218)(215,221)(216,222)(227,233)(228,234)(231,237)(232,238)(243,249)(244,250)(247,253)(248,254),
(2,9)(4,11)(6,13)(8,15)(18,25)(20,27)(22,29)(24,31)(34,41)(36,43)(38,45)(40,47)(50,57)(52,59)(54,61)(56,63)(66,73)(68,75)(70,77)(72,79)(82,89)(84,91)(86,93)(88,95)(98,105)(100,107)(102,109)(104,111)(114,121)(116,123)(118,125)(120,127)(130,137)(132,139)(134,141)(136,143)(146,153)(148,155)(150,157)(152,159)(162,169)(164,171)(166,173)(168,175)(178,185)(180,187)(182,189)(184,191)(194,201)(196,203)(198,205)(200,207)(210,217)(212,219)(214,221)(216,223)(226,233)(228,235)(230,237)(232,239)(242,249)(244,251)(246,253)(248,255),
(1,129)(2,130)(3,131)(4,132)(5,133)(6,134)(7,135)(8,136)(9,137)(10,138)(11,139)(12,140)(13,141)(14,142)(15,143)(16,144)(17,145)(18,146)(19,147)(20,148)(21,149)(22,150)(23,151)(24,152)(25,153)(26,154)(27,155)(28,156)(29,157)(30,158)(31,159)(32,160)(33,161)(34,162)(35,163)(36,164)(37,165)(38,166)(39,167)(40,168)(41,169)(42,170)(43,171)(44,172)(45,173)(46,174)(47,175)(48,176)(49,177)(50,178)(51,179)(52,180)(53,181)(54,182)(55,183)(56,184)(57,185)(58,186)(59,187)(60,188)(61,189)(62,190)(63,191)(64,192)(65,193)(66,194)(67,195)(68,196)(69,197)(70,198)(71,199)(72,200)(73,201)(74,202)(75,203)(76,204)(77,205)(78,206)(79,207)(80,208)(81,209)(82,210)(83,211)(84,212)(85,213)(86,214)(87,215)(88,216)(89,217)(90,218)(91,219)(92,220)(93,221)(94,222)(95,223)(96,224)(97,225)(98,226)(99,227)(100,228)(101,229)(102,230)(103,231)(104,232)(105,233)(106,234)(107,235)(108,236)(109,237)(110,238)(111,239)(112,240)(113,241)(114,242)(115,243)(116,244)(117,245)(118,246)(119,247)(120,248)(121,249)(122,250)(123,251)(124,252)(125,253)(126,254)(127,255)(128,256),
(1,9)(2,10)(3,11)(4,12)(5,13)(6,14)(7,15)(8,16)(17,25)(18,26)(19,27)(20,28)(21,29)(22,30)(23,31)(24,32)(33,41)(34,42)(35,43)(36,44)(37,45)(38,46)(39,47)(40,48)(49,57)(50,58)(51,59)(52,60)(53,61)(54,62)(55,63)(56,64)(65,73)(66,74)(67,75)(68,76)(69,77)(70,78)(71,79)(72,80)(81,89)(82,90)(83,91)(84,92)(85,93)(86,94)(87,95)(88,96)(97,105)(98,106)(99,107)(100,108)(101,109)(102,110)(103,111)(104,112)(113,121)(114,122)(115,123)(116,124)(117,125)(118,126)(119,127)(120,128)(129,137)(130,138)(131,139)(132,140)(133,141)(134,142)(135,143)(136,144)(145,153)(146,154)(147,155)(148,156)(149,157)(150,158)(151,159)(152,160)(161,169)(162,170)(163,171)(164,172)(165,173)(166,174)(167,175)(168,176)(177,185)(178,186)(179,187)(180,188)(181,189)(182,190)(183,191)(184,192)(193,201)(194,202)(195,203)(196,204)(197,205)(198,206)(199,207)(200,208)(209,217)(210,218)(211,219)(212,220)(213,221)(214,222)(215,223)(216,224)(225,233)(226,234)(227,235)(228,236)(229,237)(230,238)(231,239)(232,240)(241,249)(242,250)(243,251)(244,252)(245,253)(246,254)(247,255)(248,256),
(2,17)(3,33)(4,49)(5,65)(6,81)(7,97)(8,113)(9,129)(10,145)(11,161)(12,177)(13,193)(14,209)(15,225)(16,241)(19,34)(20,50)(21,66)(22,82)(23,98)(24,114)(25,130)(26,146)(27,162)(28,178)(29,194)(30,210)(31,226)(32,242)(36,51)(37,67)(38,83)(39,99)(40,115)(41,131)(42,147)(43,163)(44,179)(45,195)(46,211)(47,227)(48,243)(53,68)(54,84)(55,100)(56,116)(57,132)(58,148)(59,164)(60,180)(61,196)(62,212)(63,228)(64,244)(70,85)(71,101)(72,117)(73,133)(74,149)(75,165)(76,181)(77,197)(78,213)(79,229)(80,245)(87,102)(88,118)(89,134)(90,150)(91,166)(92,182)(93,198)(94,214)(95,230)(96,246)(104,119)(105,135)(106,151)(107,167)(108,183)(109,199)(110,215)(111,231)(112,247)(121,136)(122,152)(123,168)(124,184)(125,200)(126,216)(127,232)(128,248)(138,153)(139,169)(140,185)(141,201)(142,217)(143,233)(144,249)(155,170)(156,186)(157,202)(158,218)(159,234)(160,250)(172,187)(173,203)(174,219)(175,235)(176,251)(189,204)(190,220)(191,236)(192,252)(206,221)(207,237)(208,253)(223,238)(224,254)(240,255));

# Storage for all polytopes found
all_polys := [];
max_recursion := 20;
for i in [1..max_recursion] do
    Add(all_polys, []);
od;

# setup files
outfile := IO_File("/home/chris/fromgap.pipe", "w");
infile := IO_File("/home/chris/togap.pipe", "r");

while true do
        # read command from input
        str := IO_ReadLine(infile);
        if str <> "" then
            # Print("GAP READ: ", str);
            # Convert to GAP Object
            arr := JsonStringToGap(str);
            # Extract Level
            level := arr[1][1] + 1;
            Print("Level ", level, "\n");
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
