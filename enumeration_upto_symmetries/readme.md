# Enumeration of Bell Inequalities up to Symmetries
This directory contains files for running  *RANDA* to enumerate all faces for a two party setting of a Bell setting.
Additionally, it generates the script for GAP to check for equivalences during the enumeration.

## Basic Usage
The script *generate_files.py* creates all needed files for an enumeration of Bell Inequalities up to symmetries. 
The file is run by the command "_python generate_files.py ma mb na nb", where *ma,mb* are the number of inputs and 
*na,nb* the number of outputs. This creates files in *gap_scripts* and *randa_files*, which are used as input for GAP and RANDA.

The two programs, GAP and RANDA, communicate using so called *unnamed pipes*. These need to be created once manually by running:
```shell
cd linearbell/enumeration_upto_symmetries
mkfifo fromgap.pipe
mkfifo togap.pipe
```

Assume we want to enumerate all Bell Inequality Classes for the setting of 3 inputs and 2 outputs for both parties. 
To create the input files we run the Python script from within this directory: 

```shell
python generate_files.py 3 3 2 2
```
Now we are set for starting the enumeration. First, we start GAP for checking equivalences:
```shell
gap gap_scripts/3322.g
```
To start the enumeration by running RANDA:
```shell
randa randa_files/3322.ext > out.ine
```
Here *out.ine* will contain the enumerated Inequalities. There are multiple optional parameters that can be given to *randa*, 
which you can find [here](https://github.com/christian512/randa)

#### vertices_on_facet.py
Calculates the vertices on facets, where the FACETS are given by Toms representation. **THIS CAN BE MOVED TO ELSEWHERE, e.g. examples**.


