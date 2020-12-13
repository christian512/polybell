# Bell Inequalities, Linear Programming and Vertex Enumeration

**This is still work in progress (code can always run faster...)**

This repo is my collection of codes for finding Bell Inequalities for no-signalling probability distributions. I use different 
techniques of Linear Programming and Vertex Enumeration for finding those inequalities in different cases. One main topic 
is  finding Inequalities for the finite efficiency PR Box, which is a probability distribution relevant to finite efficient detectors.

I try to make 
the codes as reusable as possible and collect general codes in a package that might later lead to it's own repository.


## Overview

The repository is structured using different directories. Here is a small overview. For more detailed information there
is a Readme in each directory.

* **data** : contains all data generated using the codes in this repository
* **facet_enum** : codes for enumerating facets (without using vertex enumeration codes)
* **lifted_bell** : small tests for lifting bell inequalities (but mostly abandoned)
* **linearbell** : this is the package directory, where all functionality is stored that is used throughout the repository
* **notebooks** : some explanatory notebooks
* **vertex_enum** : vertex enumeration to find bell inequalities (using MPLRS)

## Software setup 
The external software I use in this repository:
* Gurobi as a fast solver for linear programs
* LRS for running vertex enumeration 
* MPLRS for running LRS in parallel

### Gurobi Setup
Some of the codes use the commercial Gurobi solver. You need a license to run the solver, but you can get a free academic
at their website.

To install Gurobi follow this guide to install the software and set the environment variables:

https://www.gurobi.com/documentation/9.1/quickstart_linux/software_installation_guid.html

For installing *gurobipy*, the python library for Gurobi, you can simply run:

*python -m pip install -i https://pypi.gurobi.com gurobipy*

More information for the pyhton installation here:

https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python-

### (MP)LRS setup
LRS is a code for vertex enumeration using an implementation of the reverse search algorithm. Additionally one can also 
use the MPLRS implementation, which runs the LRS vertex enumeration in parallel. 

Since there is no python interface for the library, I provide some simple helper functions in *linearbell/lrs_helper.py*.
This might be worth it's own repository...

The documentation and installation notes are maintained here:

http://cgm.cs.mcgill.ca/~avis/C/lrs.html
# Sources

* [*Bell nonlocality*, *Brunner et. al. (2014)*](https://arxiv.org/abs/1303.2849)
* [*Bell inequalities from no-signalling distribution*, *Cope & Colbeck (2019)*](https://arxiv.org/abs/1812.10017)
* add paper on MPLRS and others from computer






