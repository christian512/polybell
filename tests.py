from utils import general_pr_box, get_deterministic_behaviors, find_local_weight, extremal_ns_binary_vertices
from itertools import product
import numpy as np


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
    opt, bell_expression = find_local_weight(p, dets)
    assert np.abs(bell_expression @ p - 1) < 1e-9
    # EXTREMAL
    p = np.array([general_pr_box(*c) for c in output_input_combs])
    opt, bell_expression = find_local_weight(p, dets)
    assert np.abs(bell_expression @ p) < 1e-9
    # MIXTURE
    p = 1 / 2 * p + 1 / 2 * dets[0]
    opt, bell_expression = find_local_weight(p, dets)
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
        opt, bell_exp = find_local_weight(e, dets)
        assert np.abs(bell_exp @ e) < 1e-8




if __name__ == '__main__':
    test_general_pr_box_two_input_two_output_case()
    test_local_weight()
    test_extremal_points_binary_ns()
