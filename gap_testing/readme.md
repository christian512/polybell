# GAP Testing
This directory contains files for running GAP codes. Especially for checking equivalences when running PANDA to enumerate 
facets of a polytope. Hence there are also "_*.g_" files, which are scripts that can be run by GAP using "_gap file.g_".

#### GAP scripts
This folder contains the code to run the equivalence check server in GAP for different cases. Each of the files can be run by 
"_gap file.g_". The name of the file gives the scenario to which the code corresponds. The only difference in the files are the 
generators for the automorphism groups. This can be generated with a python script (see below).
In the files you might change the path of the **FIFO-pipe** files to work with your setup.

#### generators_perm_group.py
With this script you can generate the vertices of Bell Scenario and the generators of the corresponding automorphism group.
This data can then be used to run the GAP server for equivalence tests and to start a PANDA calculation for face enumeration.

#### interface_test.py
Just a test file for using the **FIFO-pipe** with GAP. **Can be deleted, I think.**

#### vertices_on_facet.py
Calculates the vertices on facets, where the FACETS are given by Toms representation. **THIS CAN BE MOVED TO ELSEWHERE, e.g. examples**.


