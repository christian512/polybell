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
* [**detection_threshold_inequalities**](./detection_threshold_inequalities/readme.md) : Calculating the detection threshold efficiency using linear programming.
* [**enumeration_upto_symmetries**](./enumeration_upto_symmetries/readme.md) : Enumerate (Bell) polytopes using different dual description algorithms.
* [**equivalence_testing**](./equivalence_testing/readme.md) : Test the equivalences of facets.
* **facet_classes** : Facet-classes of enumerated  Bell polytopes.
* [**gap_polyhedral_scripts**](./gap_polyhedral_scripts/readme.md) : Scripts for enumerating polyhedra with the [*polyhedral-package*](http://mathieudutour.altervista.org/Polyhedral/index.html) for GAP.
* [**method_comparison**](./method_comparison/readme.md) : Run-time comparison of different dual description methods applied on the Bell polytope.
* [**plots**](.plots/readme.md) : Codes for generating the plots in my thesis.
* **polybell** : All functions that are used in the other scripts to handle Bell inequalities, which can be used to write new programs.
* [**small_scripts**](./small_scripts/readme.md) : Collection of some small scripts. 

## Installation
For a simple installation process I provide a [Dockerfile](Dockerfile) that allows a quick installation using the Docker engine. 
A guide for installing Docker engine on your computer can be found [here](https://docs.docker.com/engine/install/).
This guide assumes that you use Ubuntu, but the instructions should work similarly on other operating systems.

To install all software required to run the codes in this repository, you need to clone this repository and build a docker image from the Dockerfile:

```bash
git clone https://github.com/christian512/polybell.git
cd polybell
sudo docker build -t polybell .
```
The installation process might take a while. You can now spawn a container from the generated image, which allows you to run any code from this repository. Start the container by:

```bash
sudo docker run -it polybell
```
You now have access to a Bash-Terminal within the container. You can change the directory (using the `cd` command) and run the Python-codes by:
```bash
python <script_name> <arguments>
```
To get a full list of options for all scripts, you can use the `-h` flag:

```bash
python <script_name> -h
```

Please refer to the README-files listed before, for a description of the codes.

## Manual Installation
If you want to make changes to the code, it might be useful to install all tools manually. Here we give a quick introduction 
to do so on Ubuntu.

### Python environment
I recommend the usage of a [virtual environment](https://docs.python.org/3/library/venv.html), for setting up Python. 
However, you can also proceed without it, if you have Python 3 installed. 

First install the packages given as requirements by:
```bash
pip install -r requirements.txt
```

As we are providing the functionalities of this repository as a package itself, we have to install this package directly
from the local repository. To do so, run the following command:

```bash
python -m pip install -e .
```

### RANDA and GAP installation
The Python codes use the two programs [RANDA](https://github.com/christian512/randa) and [GAP](https://www.gap-system.org/) to perform dual description algorithms for the Bell polytope. 
The installation of both can be found in the [RANDA repository](https://github.com/christian512/randa).

### Gurobi Setup

Some of the codes use the commercial Gurobi solver. You need a license to run the solver, but you can get a free
academic license at their [website](https://www.gurobi.com/).

To install Gurobi follow this guide:

https://www.gurobi.com/documentation/9.1/quickstart_linux/software_installation_guid.html

For installing *gurobipy*, the python library for Gurobi, you can simply run:

*python -m pip install -i https://pypi.gurobi.com gurobipy*

More information for the Python installation here:

https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python-







