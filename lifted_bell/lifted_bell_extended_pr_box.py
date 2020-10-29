"""
Checking if there is a lifted PR box for the finite efficiency PR box with m inputs
"""
import numpy as np
from linearbell.utils import get_deterministic_behaviors, get_possible_liftings, get_configs, general_pr_box_extended, \
    reduce_extended_pr_box, find_local_weight

# set input parameters
inputs_a = range(6)
inputs_b = range(6)
outputs_failure = range(3)
outputs_wo_failure = range(2)

# set efficiency
epsilon = 0.01
eta = 4 / (len(inputs_a) + 4) + epsilon
print('epsilon: {}'.format(epsilon))

# options for local weight optimizer (bland is only needed for higher dimensions (m_a , m_b > 4)
options = {"disp": False, "maxiter": 50000, "bland": True}
method = 'simplex'

# get deterministics for the non output
dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs_wo_failure)

# get all possible liftings for both parties
poss_lifts_a = get_possible_liftings(inputs_a, outputs_wo_failure)
poss_lifts_b = get_possible_liftings(inputs_b, outputs_wo_failure)

# get configs with and without failure
configs_wo_failure = get_configs(inputs_a, inputs_b, outputs_wo_failure, outputs_wo_failure)
configs_failure = get_configs(inputs_a, inputs_b, outputs_failure, outputs_failure)

# create extended pr_box
pr_ext = [general_pr_box_extended(a, b, x, y, eta, outputs_wo_failure) for (a, b, x, y) in configs_failure]
pr_ext = np.array(pr_ext)

# iterate through all possible combinations of liftings for A and B
niter = len(poss_lifts_a) * len(poss_lifts_b)
counter = 0
for lift_a in poss_lifts_a:
    for lift_b in poss_lifts_b:
        print("{} / {}".format(counter, niter))
        counter += 1
        # get reduced pr box under this liftings
        pr_red = reduce_extended_pr_box(pr_ext, configs_failure, configs_wo_failure, lift_a, lift_b)
        # find the local weight of the reduced pr box
        opt, bell_exp = find_local_weight(pr_red, dets, method=method, options=options, retry=False)
        # check if bell expression is correct
        if not opt.success:
            print('optimizer failed')
        if bell_exp @ pr_red < 1:
            print('Found a Bell expression')