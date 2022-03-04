# Polyhedral Algorithms for Enumerating Bell Inequality Classes

This repository includes functionality to generate facet-classes (Bell Inequality classes) for the two-party Bell
polytope, which I developed for my master thesis. To generate a single Bell inequality, I utilise the technique of
Linear Programming with the fast solvers provided by [*Gurobi*](https://www.gurobi.com/). To enumerate all Bell Inequalities
or all Bell Inequality classes I provide codes to apply dual description algorithms for convex polyhedra. Especially, I make use 
of the symmetric Adjacency Decomposition Method. The Adjacency Decomposition method is implemented in a [different repository](https://github.com/christian512/randa)
and the functionality regarding computational group theory uses functions from [*GAP*](https://www.gap-system.org/).


## Overview

The repository is structured using different directories. Here, I give an overview over the provided programs. For more detailed information there
is a *README*-file in each directory.


* [**degeneracy**](./degeneracy/readme.md) : Identifying degeneracy of the Bell polytopes.
* [**enumeration_upto_symmetries**](./enumeration_upto_symmetries/readme.md) : Enumerate (Bell) polytopes using different dual description algorithms.
* [**equivalence_testing**](./equivalence_testing/readme.md) : Test the equivalences of facets 
* **facet_classes** : Facet-classes of enumerated  Bell polytopes.
* [**gap_polyhedral_scripts**](./gap_polyhedral_scripts/readme.md) : Scripts for enumerating polyhedra with the [*polyhedral-package*](http://mathieudutour.altervista.org/Polyhedral/index.html) for GAP.

## Installation




## Manual Installation



### Python environment
We provide a *requirements.txt* file for the needed pip-packages in the virtual environment.
**Make sure you activated your new virtual environment (based on python 3)**. You can install these by:

```bash
pip install -r requirements.txt
```

As we are providing the functionalities of this repository as a package itself, we have to install this package directly
from the local repository. To do so, run the following command from the main directory, where the *setup.py* is located:

```bash
python -m pip install -e .
```

### RANDA installation

### GAP installation


### Gurobi Setup

Some of the codes use the commercial Gurobi solver. You need a license to run the solver, but you can get a free
academic at their website.

To install Gurobi follow this guide to install the software and set the environment variables:

https://www.gurobi.com/documentation/9.1/quickstart_linux/software_installation_guid.html

For installing *gurobipy*, the python library for Gurobi, you can simply run:

*python -m pip install -i https://pypi.gurobi.com gurobipy*

More information for the pyhton installation here:

https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python-







