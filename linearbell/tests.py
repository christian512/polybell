from itertools import product
import numpy as np
from linearbell.utils import *


def test_general_pr_box_two_input_two_output_case():
    # x = y = 0 --> symmetric
    assert general_pr_box(0, 0, 0, 0) == 1 / 2
    assert general_pr_box(0, 1, 0, 0) == 0
    assert general_pr_box(1, 0, 0, 0) == 0
    assert general_pr_box(1, 1, 0, 0) == 1 / 2
    # x = 0 and y = 1 --> symmetric
    assert general_pr_box(0, 0, 0, 1) == 1 / 2
    assert general_pr_box(0, 1, 0, 1) == 0
    assert general_pr_box(1, 0, 0, 1) == 0
    assert general_pr_box(1, 1, 0, 1) == 1 / 2
    # x = 1 and y = 0 --> symmetric
    assert general_pr_box(0, 0, 1, 0) == 1 / 2
    assert general_pr_box(0, 1, 1, 0) == 0
    assert general_pr_box(1, 0, 1, 0) == 0
    assert general_pr_box(1, 1, 1, 0) == 1 / 2
    # x = 1 and y = 1 --> anti-symmetric
    assert general_pr_box(0, 0, 1, 1) == 0
    assert general_pr_box(0, 1, 1, 1) == 1 / 2
    assert general_pr_box(1, 0, 1, 1) == 1 / 2
    assert general_pr_box(1, 1, 1, 1) == 0
    return True


def test_local_weight():
    """ Checks if the output for the local weight function is correct """
    # Define inputs, outputs and deterministic behavior
    inputs = range(3)
    outputs = range(2)
    output_input_combs = product(outputs, outputs, inputs, inputs)
    dets = get_deterministic_behaviors(inputs, inputs, outputs)
    # DETERMINISTIC
    p = np.copy(dets[0])
    bell_expression = find_local_weight(p, dets)
    assert np.abs(bell_expression @ p - 1) < 1e-9
    # EXTREMAL
    p = np.array([general_pr_box(*c) for c in output_input_combs])
    bell_expression = find_local_weight(p, dets)
    assert np.abs(bell_expression @ p) < 1e-9
    # MIXTURE
    p = 1 / 2 * p + 1 / 2 * dets[0]
    bell_expression = find_local_weight(p, dets)
    assert np.abs(bell_expression @ p - 1 / 2) < 1e-9


def test_extremal_points_binary_ns():
    """ Tests if the extremal points of the binary NS polytope have zero local weight """
    inputs_a = range(4)
    inputs_b = range(4)
    outputs = range(2)
    # get the extremal points
    extremals = extremal_ns_binary_vertices(inputs_a, inputs_b, outputs)
    dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs)
    for e in extremals:
        bell_exp = find_local_weight(e, dets)
        assert np.abs(bell_exp @ e) < 1e-8


def test_allowed_relabellings_antisymmetric():
    """ Test if we get the correct allowed relabellings back """
    inputs_a = range(3)
    inputs_b = range(2)
    outputs_a = range(2)
    outputs_b = range(2)
    # set numbers
    ma = len(inputs_a)
    mb = len(inputs_b)
    na = len(outputs_a)
    nb = len(outputs_b)

    # get the allowed relabellings
    allowed_perms = get_allowed_relabellings(inputs_a, inputs_b, outputs_a, outputs_b)

    fac = np.math.factorial
    theoretical_num = (fac(na) ** ma) * (fac(nb) ** mb) * fac(ma) * fac(mb)

    # check that the length of the found allowed permutations is correct.
    assert len(allowed_perms) == theoretical_num, 'Found {} permutation, theoretical num: {}'.format(len(allowed_perms),
                                                                                                     theoretical_num)


def test_allowed_relabellings_symmetric():
    """ Test if we get the correct allowed relabellings back """
    inputs_a = range(3)
    inputs_b = range(3)
    outputs_a = range(2)
    outputs_b = range(2)
    # set numbers
    ma = len(inputs_a)
    mb = len(inputs_b)
    na = len(outputs_a)
    nb = len(outputs_b)

    # get the allowed relabellings
    allowed_perms = get_allowed_relabellings(inputs_a, inputs_b, outputs_a, outputs_b)

    fac = np.math.factorial
    theoretical_num = 2 * (fac(na) ** ma) * (fac(nb) ** mb) * fac(ma) * fac(mb)

    # check that the length of the found allowed permutations is correct.
    assert len(allowed_perms) == theoretical_num, 'Found {} permutation, theoretical num: {}'.format(len(allowed_perms),
                                                                                                     theoretical_num)


def test_possible_liftings():
    """ Tests if the number of possible liftings is correct """
    inputs_a = range(5)
    outputs_a = range(4)
    # get the possible liftings for this case
    lifts_a = get_possible_liftings(inputs_a, outputs_a)
    # check the number of listings
    assert len(lifts_a) == len(outputs_a) ** len(inputs_a)


def test_possible_liftings_extended():
    """ Tests if the number of possible liftings for the extended case is correct """
    inputs_a = range(5)
    outputs_a = range(2)
    inputs_b = range(3)
    # get the liftings
    lifts = get_possible_liftings_extended(inputs_a, outputs_a, inputs_b)
    # check the number of liftings
    assert len(lifts) == len(outputs_a) ** (len(inputs_a) * len(inputs_b))


def test_parametrisation_configurations():
    """ Tests the configurations for a parametrisation """
    m = 5
    d = 9
    t = 2 * (d - 1) * m + ((d - 1) ** 2) * m ** 2
    inputs = range(m)
    outputs = range(d)
    # get the parametrisation
    configs_param = get_parametrisation_configs(inputs, inputs, outputs, outputs)
    # check that the length is correct
    assert len(configs_param) == t, 'Length of configs_param: {}'.format(len(configs_param))


def test_parametrisation():
    """ Tests the parametrisation of a behavior """
    m = 10
    d = 2
    t = 2 * (d - 1) * m + ((d - 1) ** 2) * m ** 2
    inputs = range(m)
    outputs = range(d)
    # get configs
    configs = get_configs(inputs, inputs, outputs, outputs)
    configs_param = get_parametrisation_configs(inputs, inputs, outputs, outputs)
    # get a behavior
    pr_box = np.array([general_pr_box(a, b, x, y) for (a, b, x, y) in configs])
    assert pr_box.shape[0] >= t
    # get parametrisation
    pr_box_param = parametrise_behavior(pr_box, configs, configs_param, inputs, inputs, outputs, outputs)
    assert pr_box_param.shape[0] == t
    assert len(np.unique(pr_box_param)) == 2


if __name__ == '__main__':
    # TODO: Write tests for reducing the PR box using liftings
    test_general_pr_box_two_input_two_output_case()
    test_local_weight()
    test_extremal_points_binary_ns()
    test_allowed_relabellings_symmetric()
    test_allowed_relabellings_antisymmetric()
    test_possible_liftings()
    test_possible_liftings_extended()
    test_parametrisation_configurations()
    test_parametrisation()
