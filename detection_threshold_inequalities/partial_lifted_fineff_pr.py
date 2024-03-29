"""
Checking if there is a lifted PR box for the finite efficiency PR box with m inputs
"""
import numpy as np
from polybell.utils import get_deterministic_behaviors, get_possible_liftings, get_configs, general_pr_box_extended, \
    partially_reduce_extended_pr_box, find_local_weight_dual, update_partial_lifted_bell_expression
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(dest='ma', help="number of inputs for ALICE")
parser.add_argument(dest='mb', help='number of inputs for BOB')
parser.add_argument(dest='num_no_lift_a', help="number of inputs w/0 lifted failure output for ALICE")
parser.add_argument(dest='num_no_lift_b', help="number of inputs w/0 lifted failure output for BOB")
parser.add_argument('-e', dest='efficiency_tolerance', help='Numerical tolerance for the efficency', default=0.01,
                    type=float)
args = parser.parse_args()
ma = int(args.ma)
mb = int(args.mb)
num_no_lift_a = int(args.num_no_lift_a)
num_no_lift_b = int(args.num_no_lift_b)

assert num_no_lift_a <= ma, 'Number of inputs must be larger than number of inputs w/o lifted failure output'
assert num_no_lift_b <= mb, 'Number of inputs must be larger than number of inputs w/o lifted failure output'

# set input parameters
inputs_a = range(ma)
inputs_b = range(mb)
outputs_failure = range(3)
outputs_wo_failure = range(2)

# for which inputs should the output NOT be a lifting
no_lift_a = list(range(num_no_lift_a))
no_lift_b = list(range(num_no_lift_b))

# set threshold efficiency by choosing the lower number of inputs
epsilon = args.efficiency_tolerance
eta = 4 / (min(ma, mb) + 4)
print('epsilon[ % of eta]: {} %'.format(epsilon / eta * 100))

# tolerance for checking the bell expression
tol = 1e-3

# set file
file = '{}{}{}{}_partial.txt'.format(len(inputs_a), len(inputs_b),
                                     len(outputs_wo_failure),
                                     len(outputs_wo_failure))

# get deterministics for the non output
dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs_failure)

# get all possible liftings for both parties
poss_lifts_a = get_possible_liftings(inputs_a, outputs_wo_failure)
poss_lifts_b = get_possible_liftings(inputs_b, outputs_wo_failure)

# get configs with and without failure
configs_failure = get_configs(inputs_a, inputs_b, outputs_failure, outputs_failure)

# create extended pr_box and low efficiency pr box
pr_ext = [general_pr_box_extended(a, b, x, y, eta + epsilon, outputs_wo_failure) for (a, b, x, y) in configs_failure]
pr_low_eff = [general_pr_box_extended(a, b, x, y, eta - epsilon, outputs_wo_failure) for (a, b, x, y) in
              configs_failure]
pr_ext = np.array(pr_ext)
pr_low_eff = np.array(pr_low_eff)

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
        print("{} / {}".format(counter, niter))
        counter += 1
        # set the no lift variable for party B
        lift_b = list(lift_b)
        for nb in no_lift_b:
            lift_b[nb] = -1
        # get reduced pr box under this liftings
        pr_red, del_idx = partially_reduce_extended_pr_box(pr_ext, configs_failure, lift_a, lift_b)

        # find the local weight of the reduced pr box
        bell_exp = find_local_weight_dual(pr_red, dets)
        bell = update_partial_lifted_bell_expression(bell_exp, configs_failure, lift_a, lift_b)
        # check if bell expression is correct
        if bell @ pr_red < 1 - tol and bell @ pr_low_eff > 1:
            bells.append(bell)
            print('Found a Bell expression : {}'.format(bell_exp @ pr_red))

print('Used {} as tolerance in the efficiency'.format(epsilon))
print('Used PR-Box with efficiency: {}'.format(eta + epsilon))
print('Found {} Bell Inequalities in {} lifting cases'.format(len(bells), niter))
# store the bell expressions to file
np.savetxt(file, np.array(bells))
