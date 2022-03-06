# Method Comparisons
Comparison of run-times between different methods.

## *gap_program_templates.py**
Templates for the GAP files used by the enumeration with RANDA. You can ignore this file.

## *time_ad_panda.py* 
Obtain the run-time for running adjacency decomposition with PANDA for different Bell polytopes.
Note that you need to have PANDA installed, to use this script.

## *time_dd_panda.py* 
Obtain the run-time for running the double description method with PANDA for different Bell polytopes.
Note that you need to have PANDA installed, to use this script.

## *time_enumerate_bell_polytope.py*
Obtain run-times for the algorithms that are implemented by the combination of RANDA and GAP. 
You can give different arguments to the script, which you can see by executing

```shell
python time_enumerate_bell_polytope.py -h
```
