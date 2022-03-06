# Degeneracy
Scripts for calculating if a polyhedron is non-simplical, i.e., if the facet-enumeration is degenerate.

## check_degeneracy.py
This script can check if the facet enumeration of a Bell polytope was degenerate. To check this, the polytope needs to be fully enumerated.
You have to provide the scenario *(#inputs Alice, #inputs Bob, #outputs Alice, #outputs Bob)* and the path to a file, that contains the facets 
of the particular Bell polytope. The file needs to be in the PANDA format. 

For example, you can apply this script to a polytope, whose facets are provided in the repository by calling:

```shell
python check_degeneracy.py 2 2 2 2 ../facet_classes/2222.ine
```
