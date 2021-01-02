"""
Checking if there is a lifted PR box for the finite efficiency PR box with m inputs
"""
import numpy as np
from linearbell.utils import get_deterministic_behaviors, get_possible_liftings, get_configs, general_pr_box_extended, \
    partially_reduce_extended_pr_box, find_local_weight_dual
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(dest='ma', help="number of inputs for ALICE")
parser.add_argument(dest='mb', help='number of inputs for BOB')
args = parser.parse_args()
ma = int(args.ma)
mb = int(args.mb)

# set input parameters
inputs_a = range(ma)
inputs_b = range(mb)
outputs_failure = range(3)
outputs_wo_failure = range(2)

# for which inputs should the output NOT be a lifting
no_lift_a = [2]
no_lift_b = [2]

# set efficiency
epsilon = 0.01
eta = 4 / (len(inputs_a) + 4) + epsilon
print('epsilon: {}'.format(epsilon))

# tolerance for checking the bell expression
tol = 1e-6

# set file
file = '../data/pr_box_finite_efficiency_bell_lifted/{}{}{}{}.txt'.format(len(inputs_a), len(inputs_b),
                                                                          len(outputs_wo_failure),
                                                                          len(outputs_wo_failure))

# get deterministics for the non output
dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs_failure)

# get all possible liftings for both parties
poss_lifts_a = get_possible_liftings(inputs_a, outputs_wo_failure)
poss_lifts_b = get_possible_liftings(inputs_b, outputs_wo_failure)

# get configs with and without failure
configs_failure = get_configs(inputs_a, inputs_b, outputs_failure, outputs_failure)

# create extended pr_box
pr_ext = [general_pr_box_extended(a, b, x, y, eta, outputs_wo_failure) for (a, b, x, y) in configs_failure]
pr_ext = np.array(pr_ext)

# iterate through all possible combinations of liftings for A and B
niter = len(poss_lifts_a) * len(poss_lifts_b)
counter = 0

# list of bell expressions
bells = []
for lift_a in poss_lifts_a:
    # set the no lift variable for party A
    lift_a = list(lift_a)
    for na in no_lift_a:
        lift_a[na] = -1
    for lift_b in poss_lifts_b:
        # print("{} / {}".format(counter, niter))
        counter += 1
        # set the no lift variable for party B
        lift_b = list(lift_b)
        for nb in no_lift_b:
            lift_b[nb] = -1
        # get reduced pr box under this liftings
        pr_red = partially_reduce_extended_pr_box(pr_ext, configs_failure, lift_a, lift_b)
        # find the local weight of the reduced pr box
        bell_exp = find_local_weight_dual(pr_red, dets)
        # check if bell expression is correct
        if bell_exp @ pr_red < 1 - tol:
            bells.append(bell_exp)
            print('Found a Bell expression : {}'.format(bell_exp @ pr_red))

# store the bell expressions to file
# np.savetxt(file, np.array(bells))
