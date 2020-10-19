# Bell Inequalities and Linear Programming

*This is still work in progress*

The connection between Bell inequalities and linear programming interested me, so I though a short python demo would be
good to learn about the topic by myself. This repo is just what I've done to understand it better. Maybe it helps someone.

# Notebooks
Here I provide an overview over the *Jupyter notebooks* that are provided within this repository. 

* **bell_inequalities.ipynb** - An introduction on how to find Bell inequalities by linear programming
* **binary_outcome.ipynb** - Finding Bell inequalities with arbitrary number of inputs and binary outcome. And extension to a 3-outcome model, for inefficient detectors.
* **check_facet.ipynb** - Demonstration on how to find a facet bell inequality in a 2 inputs / 2 outputs case.

# Sources

* [*Bell nonlocality*, *Brunner et. al. (2014)*](https://arxiv.org/abs/1303.2849)
* [*Bell inequalities from no-signalling distribution*, *Cope & Colbeck (2019)*](https://arxiv.org/abs/1812.10017)

### Note for development
Just some notes I don't want to forget about. Maybe I have time for that at some point..


* Currently the local deterministic behaviors are created just from the list of inputs and outputs. However we have to 
define the local deterministic behaviors in the same way as the probability distribution. I.e. the entries at the same
indices must correspond to the same setting of inputs and outputs. Maybe you can find a nicer way to set this
