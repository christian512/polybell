from itertools import product, permutations
from scipy.optimize import linprog
import numpy as np


def get_deterministic_behaviors(inputs_a, inputs_b, outputs):
    """
    TODO: For speed up one could use sparse matrices or just integer implementation with binaries
    Returns all deterministic behaviors corresponding to inputs for ALICE and BOB
    This assumes a binary outcome
    :param inputs_a: list of inputs for ALICE
    :param inputs_b: list of inputs for BOB
    :param outputs: possible outputs
    :return: list of np.arrays
    """
    # calculate dimension of each behavior
    dim = len(inputs_a) * len(inputs_b) * (len(outputs) ** 2)
    # define all hidden variables
    lhvs = product(outputs, repeat=len(inputs_a) + len(inputs_b))
    deterministics = []
    for lhv in lhvs:
        # counter for index within the behavior
        counter = 0
        # empty deterministic behavior
        d = np.zeros(dim)
        # iterate over the possible input and output combinations
        for a, b in product(outputs, outputs):
            for x, y in product(range(len(inputs_a)), range(len(inputs_b))):
                if lhv[x] == a and lhv[y + len(inputs_a)] == b:
                    d[counter] = 1.0
                counter += 1
        deterministics.append(d)
    assert len(deterministics) == len(outputs) ** (len(inputs_a) + len(inputs_b))
    return np.array(deterministics)


def general_pr_box(a, b, x, y):
    """
    Calculates the PR box probability for any input number and binary outcome
    :param a: Alices output
    :param b: Bobs output
    :param x: Alices input
    :param y: Bobs input
    :return: probability of PR box
    """
    # check the inputs
    assert x >= 0
    assert y >= 0
    assert a == 1 or a == 0 or a == -1
    assert b == 1 or b == 0 or b == -1
    # if outputs are -1, replace them with 0 internally for modulo addition
    if a == -1: a = 0
    if b == -1: b = 0
    # get binary arrays of the inputs
    x_bin = np.fromiter(map(int, np.binary_repr(x)), dtype=int)
    y_bin = np.fromiter(map(int, np.binary_repr(y)), dtype=int)
    # check the length -> update if not the same length
    if not len(x_bin) == len(y_bin):
        l = max(len(x_bin), len(y_bin))
        x_bin = np.fromiter(map(int, np.binary_repr(x, width=l)), dtype=int)
        y_bin = np.fromiter(map(int, np.binary_repr(y, width=l)), dtype=int)
    # multiply the two sides together
    t = np.dot(x_bin, y_bin) % 2
    # binary sum of the outputs
    s = (a + b) % 2
    return 1 / 2 * int(s == t)


def general_pr_box_extended(a, b, x, y, eta, inputs_a, inputs_b, outputs_without_failure):
    """
    Returns the probability distribution for a PR box, where the detectors have efficiency eta.
    The value 2 is used for value, i.e. a = 2 means that ALICE measurement has failed
    :param outputs_without_failure:
    :param inputs_b:

    :param a: ALICEs output
    :param b: BOBs output
    :param x: Alices input
    :param y: Bobs input
    :param eta: detection efficiency
    :param inputs_a: all possible inputs for ALICE
    :param inputs_b: all possible inputs for BOB
    :param outputs_without_failure: all possible outputs without the failure case = 2
    """
    assert x >= 0
    assert y >= 0
    # if both sides experience failure
    if a == 2 and b == 2:
        return (1 - eta) ** 2
    # if both sides have no failure
    elif a != 2 and b != 2:
        return (eta ** 2) * general_pr_box(a, b, x, y)
    # if only ALICE has a failure
    elif a == 2 and b != 2:
        s = np.sum([general_pr_box(a_new, b, x, y) for a_new in outputs_without_failure])
        return eta * (1 - eta) * s
    elif a != 2 and b == 2:
        s = np.sum([general_pr_box(a, b_new, x, y) for b_new in outputs_without_failure])
        return eta * (1 - eta) * s
    else:
        print('ERROR in calculation of general_pr_box_extended -> undefined outputs')
        return 0


def facet_inequality_check(deterministics, bell_expression, m_a, m_b, n, tol=1e-8):
    """
    Checks if a given bell inequality is a facet. This is done by getting all deterministic local behaviors,
    that equalize the inequality (rescaling might be needed). Then checking the dimensions, that these behaviors span.
    :return: is_facet, scaled_bell_expression, equalizing deterministics
    """
    equalizing_dets = []
    # rescale the bell expression
    # TODO: WHICH KIND OF RESCALING? THE SECOND One ensures that the smallest value is 1
    # bell_expression = bell_expression / np.min([d @ bell_expression for d in deterministics])
    bell_expression = bell_expression / np.min(bell_expression[bell_expression > tol])
    # iterate over the deterministics
    mask = deterministics @ bell_expression - 1 < tol
    equalizing_dets = deterministics[mask]

    # return if no equalizing dets
    if len(equalizing_dets) == 0:
        return False, bell_expression, np.array(equalizing_dets)
    # define the array of equalizing dets with the first subtracted
    equalizing_dets = np.array(equalizing_dets)
    eq_dets_new = equalizing_dets - equalizing_dets[0]
    # calculate the rank of the matrix
    rank = np.linalg.matrix_rank(eq_dets_new)
    # check if it's a facet by rank check
    is_facet = rank == (m_a * (n - 1) + 1) * (m_b * (n - 1) + 1) - 2
    # return is_facet and the rescaled bell expression
    return is_facet, bell_expression, np.array(equalizing_dets)


def find_bell_inequality(p, dets, method='interior-point'):
    """ Finds a Bell inequality that is violated if the behavior p is non local """
    # reformulate for the SciPy solver
    p = np.r_[p, [-1.0]]
    dets = np.c_[dets, -1.0 * np.ones(dets.shape[0])]
    # objective function and inequalities
    obj = -p
    lhs_ineq = np.append(dets, [p], axis=0)
    rhs_ineq = np.r_[np.zeros(dets.shape[0]), [1.0]]
    # run the optimizer
    opt = linprog(c=obj, A_ub=lhs_ineq, b_ub=rhs_ineq, method=method)
    # unpack the results
    s = opt.x[:-1]
    sl = opt.x[-1]
    # drop the -1.0's from p and dets to prevent changed p or dets when using after function call
    p = np.delete(p, -1, axis=0)
    dets = np.delete(dets, -1, axis=1)
    return opt, s, sl


def find_local_weight(p, dets, method='interior-point', options={"maxiter": 1000}):
    """ Finds the local weight for a behavior p """
    # objective function and inequalities
    obj = p
    lhs_ineq = np.copy(-1.0 * dets)
    rhs_ineq = -1.0 * np.ones(dets.shape[0])
    # run the optimizer
    opt = linprog(c=obj, A_ub=lhs_ineq, b_ub=rhs_ineq, method=method, options=options)
    # if not found -> retry with standard method
    if not opt.success:
        print('Could not find local weight, retry with standard params and return')
        opt = linprog(c=obj, A_ub=lhs_ineq, b_ub=rhs_ineq)
    return opt, opt.x


def extremal_ns_binary_vertices(inputs_a, inputs_b, outputs):
    """
    Returns a list of all extremal points of the no signalling set with 2 outputs.
    :param inputs_a: List of all inputs for ALICE
    :param inputs_b: List of all inputs for BOB
    :param outputs: List of all outputs for both
    :return: List of all extremal vertices
    """
    # change inputs and outputs to start with zero, to use as indices
    inputs_a = range(len(inputs_a))
    inputs_b = range(len(inputs_b))
    outputs = range(len(outputs))
    assert len(outputs) == 2, 'This only works for binary outputs'
    # symmetric and antisymmetric matrices
    S = 1 / 2 * np.eye(2)
    A = 1 / 2 * np.array([[0, 1.0], [1.0, 0]])
    # the first row and the first column are set to be symmetric matrices.
    # there are 2^((m_a-1)*(m_b-1)) options to set either sym or antisym
    # create a list of possible binary combinations
    # 0 is for symmetric and 1 is for anti-symmetric matrix
    settings = product([0, 1], repeat=((len(inputs_a) - 1) * (len(inputs_b) - 1)) - 1)
    # list of extremal points
    extremals = []
    # iterate through all possible settings of symmetric / anti-symmetric matrices
    for sett in settings:
        # append a one to setting, as the entry where x = 1 and y = 1 is anti symmetric
        s = [1] + list(sett)
        # define a behavior
        p = []
        # iterate through outputs
        for a, b, x, y in product(outputs, outputs, inputs_a, inputs_b):
            # if first row or first column
            if x == 0 or y == 0:
                p.append(S[a, b])
            # if not first row or column
            else:
                # set index for translation of matrix to vector
                idx = (x - 1) * (len(inputs_b) - 1) + (y - 1)
                # check if entry would be symmetric or anti-symm and add corresponding probability
                if s[idx] == 0:
                    p.append(S[a, b])
                if s[idx] == 1:
                    p.append(A[a, b])
        # append to all extremal points
        extremals.append(p)
    extremals = np.array(extremals)
    assert extremals.shape[0] == 2 ** ((len(inputs_a) - 1) * (len(inputs_b) - 1) - 1)
    return extremals


def check_diff_repr_same_ineq_vec(bell, perm_bells, dets, tol=1e-6):
    """
    Vectorized function
    :param bell:
    :param perm_bells:
    :param dets:
    :param tol:
    :return:
    """
    v1 = dets @ bell
    v2s = np.dot(dets, perm_bells.transpose()).transpose()
    v2 = v2s[0]

    # get the second largest value of v1 and v2
    # as all v2s are just permuted it's always the same second largest value
    s1 = np.min(v1[v1 > np.min(v1) + tol])
    s2 = np.min(v2[v2 > np.min(v2) + tol])

    # rescale
    v1_new = v1 / (s1 - 1) + (s1 - 2) / (s1 - 1)
    v2s_new = v2s / (s2 - 1) + (s2 - 2) / (s2 - 1)

    # check if any two vectors are the same
    return np.any(np.sum((v1 - v2s) ** 2, axis=1) < tol ** 2)


def check_diff_repr_same_ineq(bell1, bell2, dets, tol=1e-6):
    """
    Checks whether two bell expression represent the same inequality
    :param bell1: scaled first bell expression
    :param bell2: scaled second bell expression
    :param dets: deterministic behaviors
    :return: boolean
    """
    # get decimal tolerance from tolerance
    dec = int(-1 * np.log10(tol))
    # Generate Bell expressions
    v1 = dets @ bell1
    v2 = dets @ bell2
    # round the arrays
    # v1 = np.round(v1, decimals=dec)
    # v2 = np.round(v2, decimals=dec)

    # check that minimal value is 1
    # assert np.abs(np.min(v1) - 1.0) < tol, 'Minimum of v1 is not 1, did you give scaled bell expression?'
    # assert np.abs(np.min(v2) - 1.0) < tol, 'Minimum of v2 is not 1, did you give scaled bell expression?'
    """
    # check if there is more than one unique element in the array
    if len(np.unique(v1)) == 1 or len(np.unique(v2)) == 1:
        # if both have only one entry, it is one, thus the same representation
        if len(np.unique(v1)) == len(np.unique(v2)):
            return True
        # if one array has only one unique value and the other more than 1, they can not be the same
        return False
    """

    # get the second largest value
    s1 = np.min(v1[v1 > np.min(v1) + tol])
    s2 = np.min(v2[v2 > np.min(v2) + tol])
    # assert s1 > 1.0
    # assert s2 > 1.0
    # rescale
    v1_new = v1 / (s1 - 1) + (s1 - 2) / (s1 - 1)
    v2_new = v2 / (s2 - 1) + (s2 - 2) / (s2 - 1)
    # check that the second largest value is 1
    # assert np.min(v1_new[v1_new > np.min(v1_new) + tol]) == 2.0, 'Second smallest value is not 2.0'
    # assert np.min(v2_new[v2_new > np.min(v2_new) + tol]) == 2.0, 'Second smallest value is not 2.0'
    # if the two new vectors are the same, they represent the same inequality
    return np.sum((v1_new - v2_new) ** 2) < tol ** 2


def check_equiv_bell(bell1, bell2, allowed_perms, dets, tol=1e-6):
    """
    Checks if two bell inequalities are equivalent under relabelling and for each perm checks that if they are equiv
    :param bell1: first bell expression
    :param bell2: second bell expression
    :param allowed_perms: allowed permutations generated by get_allowed_permutations
    :param dets: deterministic behaviors
    :param tol: tolerance
    :return:
    """
    # check if they are the same
    if np.sum((bell1 - bell2) ** 2) < tol ** 2:
        return True

    # iterate over all allowed permutations
    if np.any(np.sum((bell1 - bell2[allowed_perms]) ** 2, axis=1) < tol ** 2):
        return True

    # check if the two bell expressions are equivalent
    if check_diff_repr_same_ineq_vec(bell1, bell2[allowed_perms], dets, tol=tol):
        return True
    """
    for perm in allowed_perms:
        if check_diff_repr_same_ineq(bell1, bell2[perm], dets):
            print('Different representation of the same bell inequality')
            return True
    """

    return False


def check_relabel(bell1, bell2, allowed_perms, tol=1e-6):
    """
    Checks if two bell expressions are equivalent under relabelling
    :param bell1:
    :param bell2:
    :param allowed_perms:
    :return:
    """
    # iterate through all possible permutations
    for perm in allowed_perms:
        # check if they bell1 is equal to permuted bell2 expression
        if np.allclose(bell1, bell2[perm], atol=tol):
            return True
    return False


def get_allowed_relabellings(inputs_a, inputs_b, outputs_a, outputs_b):
    """
    Get the allowed relabellings for given inputs
    :param inputs_a: inputs for ALICE
    :param inputs_b: inputs for BOB
    :param outputs_a: outputs for ALICE
    :param outputs_b: output for BOB
    :return: list of permutation indices
    """
    # Check if interchanging of sides is possible
    interchange = False
    if len(inputs_a) == len(inputs_b) and len(outputs_a) == len(outputs_b):
        interchange = True

    # list for all allowed permutations
    allowed_permutations = []
    # list of all configurations
    configurations = [(a, b, x, y) for a, b, x, y in product(outputs_a, outputs_b, inputs_a, inputs_b)]
    # iterate over all possible relabellings for x
    for relabel_x in permutations(range(len(inputs_a))):
        # iterate over all possible relabellings for y
        for relabel_y in permutations(range(len(inputs_b))):
            # iterate over all possible relabellings of a
            for relabel_a_list in product(list(permutations(range(len(outputs_a)))), repeat=len(inputs_a)):
                # ierate over all possible relabellings of b
                for relabel_b_list in product(list(permutations(range(len(outputs_b)))), repeat=len(inputs_b)):
                    # create empty permutation of dimension of a behavior
                    perm = list(range(len(configurations)))
                    perm_sides_interchanged = list(range(len(configurations)))
                    # iterate over all configs and relabel according to currently chosen relabelling
                    for old_idx, config in enumerate(configurations):
                        # get current config
                        a_old, b_old, x_old, y_old = config
                        # get the relabelings for the outputs according to the input
                        relabel_a = relabel_a_list[x_old]
                        relabel_b = relabel_b_list[y_old]
                        # get new labels for x and y
                        a_new, b_new, x_new, y_new = relabel_a[a_old], relabel_b[b_old], relabel_x[x_old], relabel_y[
                            y_old]
                        # find new index
                        new_idx = configurations.index((a_new, b_new, x_new, y_new))
                        # set new idx in the permutation
                        perm[old_idx] = new_idx
                        # and add permutation where Alice and Bob are interchanged
                        if interchange:
                            new_idx = configurations.index((b_new, a_new, y_new, x_new))
                            perm_sides_interchanged[old_idx] = new_idx
                    allowed_permutations.append(perm)
                    if interchange:
                        allowed_permutations.append(perm_sides_interchanged)

    return allowed_permutations
