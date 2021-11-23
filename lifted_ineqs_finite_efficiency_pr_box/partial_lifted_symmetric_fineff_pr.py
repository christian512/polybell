"""
Checking if there is a lifted PR box for the finite efficiency PR box with m inputs
"""
import numpy as np
from linearbell.utils import get_deterministic_behaviors, get_possible_liftings, get_configs, general_pr_box_extended, \
    partially_reduce_extended_pr_box, find_local_weight_dual, update_partial_lifted_bell_expression
import argparse
from linearbell.representations import get_configs_mat, transform_vec_to_mat

parser = argparse.ArgumentParser()
parser.add_argument(dest='m', help="number of inputs for ALICE and BOB")
parser.add_argument(dest='num_no_lift', help='number of inputs w/o lifted failure output for ALICE and BOB')
args = parser.parse_args()
m = int(args.m)
num_no_lift = int(args.num_no_lift)

assert num_no_lift <= m, 'Number of inputs must be larger than number of inputs w/o lifted failure output'

# set input parameters
inputs_a = range(m)
inputs_b = range(m)
outputs_failure = range(3)
outputs_wo_failure = range(2)

# for which inputs should the output NOT be a lifting
no_lift = list(range(num_no_lift))

# set efficiency
epsilon = 0.01
eta = 4 / (len(inputs_a) + 4)
print('epsilon[ % of eta]: {} %'.format(epsilon / eta * 100))

# tolerance for checking the bell expression
tol = 1e-3

# set file
file = '../data/pr_box_finite_efficiency_bell_lifted_partial/{}{}{}{}.txt'.format(len(inputs_a), len(inputs_b),
                                                                                  len(outputs_wo_failure),
                                                                                  len(outputs_wo_failure))

# get deterministics for the non output
dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs_failure)

# get all possible liftings for both parties
poss_lifts_a = get_possible_liftings(inputs_a, outputs_wo_failure)

# get configs with and without failure
configs_failure = get_configs(inputs_a, inputs_b, outputs_failure, outputs_failure)

# create extended pr_box and low efficiency pr box
pr_ext = [general_pr_box_extended(a, b, x, y, eta + epsilon, outputs_wo_failure) for (a, b, x, y) in configs_failure]
pr_low_eff = [general_pr_box_extended(a, b, x, y, eta - epsilon, outputs_wo_failure) for (a, b, x, y) in
              configs_failure]
pr_ext = np.array(pr_ext)
pr_low_eff = np.array(pr_low_eff)

# iterate through all possible combinations of liftings for A and B
niter = len(poss_lifts_a)
counter = 0

# list of bell expressions
bells = []
for lift_a in poss_lifts_a:
    # set the no lift variable for party A
    lift_a = list(lift_a)
    for na in no_lift:
        lift_a[na] = -1
    print("{} / {}".format(counter, niter))
    counter += 1

    # get reduced pr box under this liftings
    pr_red, del_idx = partially_reduce_extended_pr_box(pr_ext, configs_failure, lift_a, lift_a)

    # try to find a bell inequality
    # find the local weight of the reduced pr box
    bell_exp = find_local_weight_dual(pr_red, dets)
    bell = update_partial_lifted_bell_expression(bell_exp, configs_failure, lift_a, lift_a)
    # check if bell expression is correct
    if bell @ pr_red < 1 - tol and bell @ pr_low_eff > 1:
        bells.append(bell)
        print('Found a Bell expression : {}'.format(bell_exp @ pr_red))

print('Found {} in {} lifting cases'.format(len(bells), niter))