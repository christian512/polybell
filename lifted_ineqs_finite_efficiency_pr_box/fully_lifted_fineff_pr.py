"""
Checking if there is a lifted PR box for the finite efficiency PR box with m inputs
"""
import numpy as np
from linearbell.utils import get_deterministic_behaviors, get_possible_liftings, get_configs, general_pr_box_extended, \
    reduce_extended_pr_box, find_local_weight_dual
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

# set efficiency
epsilon = 0.01
eta = 4 / (len(inputs_a) + 4)
print('epsilon: {}'.format(epsilon))

# tolerance for checking the bell expression
tol = 1e-6

# set file
file = '../data/pr_box_finite_efficiency_bell_lifted/{}{}{}{}.txt'.format(len(inputs_a), len(inputs_b), len(outputs_wo_failure), len(outputs_wo_failure))

# get deterministics for the non output
dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs_wo_failure)

# get all possible liftings for both parties
poss_lifts_a = get_possible_liftings(inputs_a, outputs_wo_failure)
poss_lifts_b = get_possible_liftings(inputs_b, outputs_wo_failure)

# get configs with and without failure
configs_wo_failure = get_configs(inputs_a, inputs_b, outputs_wo_failure, outputs_wo_failure)
configs_failure = get_configs(inputs_a, inputs_b, outputs_failure, outputs_failure)

# create extended pr_box
pr_ext = [general_pr_box_extended(a, b, x, y, eta + epsilon, outputs_wo_failure) for (a, b, x, y) in configs_failure]
pr_ext = np.array(pr_ext)

pr_low_eff = [general_pr_box_extended(a, b, x, y, eta - epsilon, outputs_wo_failure) for (a, b, x, y) in configs_failure]
pr_low_eff = np.array(pr_low_eff)

# iterate through all possible combinations of liftings for A and B
niter = len(poss_lifts_a) * len(poss_lifts_b)
counter = 0

# list of bell expressions
bells = []
for lift_a in poss_lifts_a:
    for lift_b in poss_lifts_b:
        print("{} / {}".format(counter, niter))
        counter += 1
        # get reduced pr box under this liftings
        pr_red = reduce_extended_pr_box(pr_ext, configs_failure, configs_wo_failure, lift_a, lift_b)
        pr_red_low_eff = reduce_extended_pr_box(pr_low_eff, configs_failure, configs_wo_failure, lift_a, lift_b)

        # find the local weight of the reduced pr box
        bell_exp = find_local_weight_dual(pr_red, dets)
        # check if bell expression is correct
        if bell_exp @ pr_red < 1 - tol and bell_exp @ pr_red_low_eff > 1:
            bells.append(bell_exp)
            print('Found a Bell expression : {}'.format(bell_exp @ pr_red))

print('Found {} Bell Inequalities in {} lifting cases'.format(len(bells), niter))
# store the bell expressions to file
np.savetxt(file, np.array(bells))