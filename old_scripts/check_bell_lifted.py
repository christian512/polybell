"""
Checks whether a bell inequality is of lifted form.
For this we have to check if there is another input added to a Bell Inequality for a lower dimensional case.

Method:
Remove one input from the Bell Inequality and check whether the result is a Bell Inequality in the lower dimensional case.

Result:
It seems that the Bell Inequalities that the Heuristic ADM algorithm missed out on are not of lifted form.

"""
from polybell.utils import get_configs, get_deterministic_behaviors
import numpy as np

inputs = range(4)
inputs_lower = range(3)
outputs = range(2)

# get configurations
configs = get_configs(inputs, inputs, outputs, outputs)
dets = get_deterministic_behaviors(inputs, inputs, outputs)

# get configurations for lower dimensional bell polytope
configs_reduced = get_configs(inputs_lower, inputs, outputs, outputs)
dets_reduced = get_deterministic_behaviors(inputs_lower, inputs, outputs)

# load Bell Inequality, should be in form b * p >= 1
file = "../data/facets/heuristic_adm_4422_notfound.txt"
ineqs = np.loadtxt(file)


# iterate over inequalities
for bell in ineqs:
    new_bell = np.zeros(len(configs_reduced))
    # drop one input of second party
    for drop_input in inputs:
        for i, c in enumerate(configs):
            a, b, x, y = c
            if x == drop_input:
                continue
            if x > drop_input:
                x = x - 1
            idx = configs_reduced.index((a,b,x,y))
            new_bell[idx] = bell[i]

        lifting = True
        for d in dets_reduced:
            if d @ new_bell < 1 - 1e-3:
                lifting = False
                break
        if lifting:
            print('Found a lifting!')
