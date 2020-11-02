"""
Checking if there is a lifted PR box for the finite efficiency PR box with m inputs
"""
import numpy as np
from linearbell.utils import get_deterministic_behaviors, get_possible_liftings_extended, get_configs, \
    general_pr_box_extended, \
    reduce_extended_pr_box_extended_lifts, find_local_weight
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
eta = 4 / (len(inputs_a) + 4) + epsilon
print('epsilon: {}'.format(epsilon))

# tolerance for checking the bell expression
tol = 1e-6

# set file
file = '../data/pr_box_finite_efficiency_bell_lifted_joint/{}{}{}{}.txt'.format(len(inputs_a), len(inputs_b),
                                                                                len(outputs_wo_failure),
                                                                                len(outputs_wo_failure))

# get deterministics for the non output
dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs_wo_failure)

# get all possible liftings for both parties
poss_lifts_a = get_possible_liftings_extended(inputs_a, outputs_wo_failure, inputs_b)
poss_lifts_b = get_possible_liftings_extended(inputs_b, outputs_wo_failure, inputs_a)

# get configs with and without failure
configs_wo_failure = get_configs(inputs_a, inputs_b, outputs_wo_failure, outputs_wo_failure)
configs_failure = get_configs(inputs_a, inputs_b, outputs_failure, outputs_failure)

# create extended pr_box
pr_ext = [general_pr_box_extended(a, b, x, y, eta, outputs_wo_failure) for (a, b, x, y) in configs_failure]
pr_ext = np.array(pr_ext)

# iterate through all possible combinations of liftings for A and B
niter = len(poss_lifts_a) * len(poss_lifts_b)
counter = 0

# TODO: open file
f = open(file, 'w+')

# identifier if a bell inequality was found
bell_found = False

# iterate through all possible liftings
for lift_a in poss_lifts_a:
    for lift_b in poss_lifts_b:
        print("{} / {}".format(counter, niter))
        counter += 1
        # get reduced pr box under this liftings
        pr_red = reduce_extended_pr_box_extended_lifts(pr_ext, configs_failure, configs_wo_failure, lift_a, lift_b)
        # find the local weight of the reduced pr box
        bell_exp = find_local_weight(pr_red, dets)
        assert np.abs(np.min(dets @ bell_exp) - 1) < tol, 'Bell * deterministic is not 1 at minimum'
        # check if bell expression is correct
        if not bell_exp @ pr_red < 1 - tol:
            # assert that the liftings is of the form, where it only depends on the input of the choosing party
            # This means that there is only one unique entry per row in the lifting arrays
            if not np.unique(lift_a, axis=1).shape == (ma, 1):
                f.write('Could not find a bell inequality in a non simple lifted case. \n')
                f.write('Lifting matrix for ALICE: {}'.format(lift_a))
            if not np.unique(lift_b, axis=1).shape == (mb, 1):
                f.write('Could not find a bell inequality in a non simple lifted case. \n')
                f.write('Lifting matrix for BOB: {}'.format(lift_a))
            # quit in this case by assertion
            assert np.unique(lift_a, axis=1).shape == (ma, 1)
            assert np.unique(lift_b, axis=1).shape == (mb, 1)
        else:
            # assert that one of the liftings actually exploits the joint/extended definition of lifting
            # so it's not just the simple case of a lifting
            assert np.unique(lift_a, axis=1).shape != (ma, 1) or np.unique(lift_b, axis=1).shape != (mb, 1)
            print('Found a bell inequality')
            if not bell_found:
                bell_found = True
                f.write('Found a Bell inequality for this case \n')
f.write('Finished iteration over all possible extended liftings')
f.close()
