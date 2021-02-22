from itertools import product, permutations
from scipy.optimize import linprog
import numpy as np
import gurobipy as gp
from gurobipy import GRB
from scipy import sparse

# set quiet environment for guro solver
env = gp.Env(empty=True)
env.setParam('OutputFlag', 0)
env.start()


def get_configs(inputs_a, inputs_b, outputs_a, outputs_b):
    """
    Returns the settings of a,b,x,y for each index in a probability distribution
    """
    configs = [(a, b, x, y) for a, b, x, y in product(outputs_a, outputs_b, inputs_a, inputs_b)]
    return configs


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


def get_possible_liftings(inputs, outputs):
    """
    Gets all possible liftings for given inputs and outputs.
    The return list contains for each input, which output was used to create the new output.
    """
    return list(product(outputs, repeat=len(inputs)))


def get_possible_liftings_extended(inputs, outputs, inputs_other):
    """
    This gets the possible liftings, if the output that we duplicate for the new output is dependend on both inputs.
    So if we lift the bell expression for Alice side (i.e. add one output to her measurement results) the origin location
    for copying that value depends on Alice AND on BOBs input.
    The output is a list of 2D arrays. Each array gives in the entry i,j from which output the new output was created
    """
    lifts = list(product(outputs, repeat=len(inputs) * len(inputs_other)))
    # reshape each element in the list
    for i in range(len(lifts)):
        lifts[i] = np.array(lifts[i]).reshape((len(inputs), len(inputs_other)))
    return lifts


# TODO Write function as getting determnisistic behaviors, so just input, inputs and outputs
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


def general_pr_box_extended(a, b, x, y, eta, outputs_without_failure, failure_indicator=2):
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
    assert failure_indicator not in outputs_without_failure
    # if both sides experience failure
    if a == failure_indicator and b == failure_indicator:
        return (1 - eta) ** 2
    # if both sides have no failure
    elif a != failure_indicator and b != failure_indicator:
        return (eta ** 2) * general_pr_box(a, b, x, y)
    # if only ALICE has a failure
    elif a == failure_indicator and b != failure_indicator:
        s = np.sum([general_pr_box(a_new, b, x, y) for a_new in outputs_without_failure])
        return eta * (1 - eta) * s
    elif a != failure_indicator and b == failure_indicator:
        s = np.sum([general_pr_box(a, b_new, x, y) for b_new in outputs_without_failure])
        return eta * (1 - eta) * s
    else:
        print('ERROR in calculation of general_pr_box_extended -> undefined outputs')
        return 0


def reduce_extended_pr_box(pr_box, configs_ext, configs_red, lift_a, lift_b, failure_indicator=2):
    """
    # TODO: rename this function as we can reduce a general behavior with it (also deterministic points)
    Reduces the PR box probability distribution under liftings.
    The PR-Box might have 3 outputs, where one output identifies a failure (failure_indicator).
    We assume that this output is generated by a lifting. Thus we can reduce it by inverting the lifting.
    This is done by checking where it originates from and adding the probability to the origin probability.
    """
    # check dimensions of the liftings
    assert np.ndim(
        lift_a) == 1, 'The lifting you gave is not 1D. Do you want to use the extended version of this function?'
    assert np.ndim(
        lift_b) == 1, 'The lifting you gave is not 1D. Do you want to use the extended version of this function?'

    # check length
    assert len(configs_ext) == len(pr_box)
    # create empty reduced pr box
    pr_red = np.zeros(len(configs_red))
    # iterate through the configs
    for i, c in enumerate(configs_ext):
        # get explicit configuration
        a, b, x, y = c
        # check if a or b is a failure
        if a == failure_indicator:
            a = lift_a[x]
        if b == failure_indicator:
            b = lift_b[y]
        # get the index where this new configuration is located
        idx = configs_red.index((a, b, x, y))
        # add the probability of the extended PR box to this example
        pr_red[idx] += pr_box[i]
    return pr_red


def partially_reduce_extended_pr_box(pr_box, configs_ext, lift_a, lift_b, failure_indicator=2):
    """
    Reduces a PR box under liftings. However here some entries of lift_a / lift_b might be "-1", which implies that this
    output should not be a lifting, but stay an independent variable. This way we get a weirdly shaped reduced probability
    distribution. The configs are also returned with this function
    Parameters
    ----------
    pr_box: np.array - PR box prob dist
    configs_ext: list - gives config for each entry in pr_box
    lift_a: liftings array generated with get_possible_liftings
    lift_b: liftings array generated with get_possible_liftings
    failure_indicator: indicator for failure (default 2)

    Returns
    -------
    pr_red : np.array - reduced PR box under given liftings (some may excluded)
    configs_red: list - configurations of the new PR box, might be a weird ordering.
    """
    # check dimensions of the liftings
    assert np.ndim(
        lift_a) == 1, 'The lifting you gave is not 1D. Do you want to use the extended version of this function?'
    assert np.ndim(
        lift_b) == 1, 'The lifting you gave is not 1D. Do you want to use the extended version of this function?'

    # check length
    assert len(configs_ext) == len(pr_box)
    # create empty reduced pr box
    pr_red = np.zeros(len(configs_ext))
    # indices that are copied to another value -> they have to be deleted as there is no dof anymore
    del_idx = []
    # iterate through the configs
    for i, c in enumerate(configs_ext):
        # get explicit configuration
        a, b, x, y = c
        # check if a or b is a failure
        if a == failure_indicator:
            # check that the output for this input x should be a lifting
            if not lift_a[x] == -1:
                a = lift_a[x]
        if b == failure_indicator:
            if not lift_b[y] == -1:
                b = lift_b[y]
        # get the index where this new configuration is located or append it if not yet in the list
        idx = configs_ext.index((a, b, x, y))
        # add the probability of the extended PR box to this example
        pr_red[idx] += pr_box[i]
        # append index to delete indices
        # TODO: Deleted indices are not anymore used
        if not idx == i:
            del_idx.append(i)
    return pr_red, del_idx


def update_partial_lifted_bell_expression(bell_red, configs, lift_a, lift_b, failure_indicator=2):
    """
    This function updates a Bell expression that was found for a partially reduced PR box by the function
    partially_reduce_extended_pr_box(). If the PR box is reduced under a lifting, we add probabilities from the new output
    to the old output, where it was generated from. With adding these probabilities we put the assumption that the values
    in the bell expressions are generated from a lifting. So we have to copy the value from the old output to the new output.
    Parameters
    ----------
    failure_indicator
    bell_red
    configs
    lift_a
    lift_b

    Returns
    -------
    bell_expression
    """
    assert np.ndim(lift_a) == 1
    assert np.ndim(lift_b) == 1
    assert len(configs) == len(bell_red)
    # iterate through configurations
    for i, c in enumerate(configs):
        # get the explicit configuration
        a, b, x, y = c
        # check if it's a failure indicator
        if a == failure_indicator:
            # check that the output for this input x should be a lifting
            if not lift_a[x] == -1:
                a = lift_a[x]
        # do the same for b
        if b == failure_indicator:
            if not lift_b[y] == -1:
                b = lift_b[y]
        # get the index of the new config
        idx = configs.index((a, b, x, y))
        # if the new idx is different from the old index -> it's a lifting
        # then we have to move the value in the bell expression also to the
        # lifted output
        # TODO: Do we have to do i == idx
        bell_red[i] = bell_red[idx]
    return bell_red


def reduce_extended_pr_box_extended_lifts(pr_box, configs_ext, configs_red, lift_a, lift_b, failure_indicator=2):
    """
    # TODO: rename this function as we can reduce a general behavior with it (also deterministic points) 
    This function does similar things to reduce_extended_pr_box(). However this takes liftings that are dependent on the
    input of both parties.
    """
    # check dimensions of liftings
    assert np.ndim(
        lift_a) == 2, 'The lifting you gave is not 2D. Do you want to use the reduced version of this function?'
    assert np.ndim(
        lift_b) == 2, 'The lifting you gave is not 2D. Do you want to use the reduced version of this function?'
    # check length
    assert len(configs_ext) == len(pr_box)
    # create empty reduced pr box
    pr_red = np.zeros(len(configs_red))
    # iterate through configs
    for i, c in enumerate(configs_ext):
        # get explicit configuration
        a, b, x, y = c
        # check if a or b is a failure
        if a == failure_indicator:
            # set the output to the original one -> where this output was cloned from
            a = lift_a[x, y]
        if b == failure_indicator:
            # set the output to the original one -> where this output was cloned from
            b = lift_b[y, x]
        # get idx where the configuration originated from
        idx = configs_red.index((a, b, x, y))
        pr_red[idx] += pr_box[i]
    return pr_red


def facet_inequality_check(deterministics, bell_expression, m_a, m_b, n, tol=1e-4):
    """
    Checks if a given bell inequality is a facet. This is done by getting all deterministic local behaviors,
    that equalize the inequality (rescaling might be needed). Then checking the dimensions, that these behaviors span.
    :return: is_facet, scaled_bell_expression, equalizing deterministics
    """
    equalizing_dets = []
    # factor for rescaling
    fac = np.min(deterministics @ bell_expression)
    # or fac =  np.min(bell_expression[bell_expression > tol])
    if np.abs(fac - 1) > tol:
        print('rescale factor is not close to 1: fac = {}'.format(fac))
        # rescale
        bell_expression = bell_expression / fac
    # iterate over the deterministics
    mask = np.abs(deterministics @ bell_expression - 1) < tol
    equalizing_dets = deterministics[mask]
    assert np.max(np.abs(equalizing_dets @ bell_expression - 1)) < tol, 'The equalizing dets are not equalizing ?!'
    # TODO: I wanted to do the assertion with .all() but this somehow did not work. Maybe try to make a minimal example and check out whats wrong

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
    print('WATCH OUT, YOU RUN THE SLOW SCIPY LINPROG SOLVER')
    opt = linprog(c=obj, A_ub=lhs_ineq, b_ub=rhs_ineq, method=method)
    # unpack the results
    s = opt.x[:-1]
    sl = opt.x[-1]
    # drop the -1.0's from p and dets to prevent changed p or dets when using after function call
    p = np.delete(p, -1, axis=0)
    dets = np.delete(dets, -1, axis=1)
    return opt, s, sl


def find_local_weight_scipy(p, dets, method='interior-point', options={"maxiter": 1000}, retry=True):
    """ Finds the local weight for a behavior p """
    # objective function and inequalities
    obj = p
    lhs_ineq = np.copy(-1.0 * dets)
    rhs_ineq = -1.0 * np.ones(dets.shape[0])
    # run the optimizer
    opt = linprog(c=obj, A_ub=lhs_ineq, b_ub=rhs_ineq, method=method, options=options)
    # if not found -> retry with standard method
    if not opt.success and retry:
        print('Could not find local weight, retry with standard params and return')
        opt = linprog(c=obj, A_ub=lhs_ineq, b_ub=rhs_ineq)
    return opt, opt.x


def find_local_weight_primal(p, dets, method=-1, tol=1e-6):
    """ Finds the local weights for a behavior and given deterministic points """
    # create a model
    m = gp.Model('local_weight_dual', env=env)
    m.setParam("Method", method)
    # add variable
    weight = m.addMVar(shape=dets.shape[0], name='weight')
    # set objective
    ones = np.ones(dets.shape[0])
    m.setObjective(weight @ ones, GRB.MAXIMIZE)
    # setup variables for constraints
    dets_t = np.transpose(dets)
    lhs = dets_t @ weight
    lhs = np.sum(np.transpose(lhs), axis=0)
    # add constraint
    m.addConstr(lhs <= p, name='c')
    m.addConstr(weight >= np.zeros(dets.shape[0]), name='pos')
    # optimize
    m.optimize()
    # return
    return weight.X


def find_local_weight_dual(p, dets, method=-1, tol=1e-6):
    # create a model
    m = gp.Model('local_weight_dual', env=env)
    m.setParam("Method", method)
    # add variable
    bell = m.addMVar(shape=p.shape[0], name='bell')
    # set objective
    m.setObjective(p @ bell, GRB.MINIMIZE)
    # setup variables for constraints
    rhs = np.ones(dets.shape[0]) - tol
    A = sparse.csr_matrix(dets)
    # add constraint
    m.addConstr(A @ bell >= rhs, name='c')
    # optimize
    m.optimize()
    # return
    return bell.X


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
    :return: boolean
    """

    v1 = dets @ bell
    v2s = np.dot(dets, perm_bells.transpose()).transpose()
    # get second smallest value for v1
    s1 = np.min(v1[v1 > np.min(v1) + tol])
    # Get second smallest value for each row in v2s
    # TODO: Find a more efficient way to do this for loop
    s2s = np.empty(v2s.shape[0])
    for i in range(v2s.shape[0]):
        v2 = v2s[i]
        s2s[i] = np.min(v2[v2 > np.min(v2) + tol])

    # reformat s2s to perform the scaling
    s2s = s2s[:, np.newaxis]
    # rescaling
    v1_new = v1 / (s1 - 1) + (s1 - 2) / (s1 - 1)
    v2s_new = v2s / (s2s - 1) + (s2s - 2) / (s2s - 1)

    # check if any two vectors are the same
    return np.any(np.sum((v1_new - v2s_new) ** 2, axis=1) < tol ** 2)


def check_equiv_bell(bell1, bell2, relabels_dets, dets, tol=1e-6):
    """
    Checks if two bell inequalities are equivalent under relabelling and for each perm checks that if they are equiv
    :param bell1: first bell expression
    :param bell2: second bell expression
    :param relables_dets: permutations of the deterministic behaviors -> which det maps to which det
    :param dets: deterministic behaviors
    :param tol: tolerance
    :return:
    """
    # check if they are the same
    if np.sum((bell1 - bell2) ** 2) < tol ** 2:
        return True
    # get the v vectors
    v1 = dets @ bell1
    v2 = dets @ bell2
    # get the second smallest values bigger than 1
    if np.min(np.abs(v1 - 1.0)) > tol:
        print('error in v1, minimum is not 1, but {}'.format(np.min(np.abs(v1))))
        print(v1)
        return True
    if np.min(np.abs(v2 - 1.0)) > tol:
        print('error in v2, minimum is not 1, but {}'.format(np.min(np.abs(v2))))
        print(v2)
        return True
    try:
        s1 = np.min(v1[v1 > 1.0 + tol])
        v1 = v1 / (s1 - 1) + (s1 - 2) / (s1 - 1)
    except:
        pass
    try:
        s2 = np.min(v2[v2 > 1.0 + tol])
        v2 = v2 / (s2 - 1) + (s2 - 2) / (s2 - 1)
    except:
        pass

    if np.sum((v1 - v2) ** 2) < tol: return True
    # try to see if they have the same tally
    u1, c1 = np.unique(np.round(v1, decimals=1), return_counts=True)
    u2, c2 = np.unique(np.round(v2, decimals=1), return_counts=True)
    if not u1.shape[0] == u2.shape[0]: return False
    if not np.all(u1 == u2): return False
    if not np.all(c1 == c2): return False
    # check if any relabelling is the same
    return np.any(np.sum((v1 - v2[relabels_dets]) ** 2, axis=1) < tol ** 2)


def check_equiv_bell_vertex_enum(bell1, bell2, relabels, dets, tol=1e-6):
    """
    Checks if two bell inequalities are equivalent under relabelling and for each perm checks that if they are equiv
    This function uses the relabels of the bell expressions and therefore needs less memory, but has to calculate the
    v vectors multiple times.
    We assume that the bell expressions are already rescaled by the function affine_transform_bell
    :param bell1: first bell expression
    :param bell2: second bell expression
    :param relabels: permutations of the bell expression that are allowed -> from get_allowed_relabellings
    :param dets: deterministic behaviors
    :param tol: tolerance
    :return:
    """
    # check if they are the same
    if np.sum((bell1 - bell2) ** 2) < tol ** 2:
        return True
    # get the v vectors
    v1 = dets @ bell1
    v2 = dets @ bell2
    # check that smallest value is zero, due to affine_transformation
    assert np.abs(np.min(v1)) < tol, 'min(v1): {}'.format(np.min(v1))
    assert np.abs(np.min(v1[v1 > tol]) - 1.0) < tol, 'min(v1[v1 > 0]): {}'.format(np.min(v1[v1 > 0]))
    assert np.abs(np.min(v2)) < tol, 'min(v2): {}'.format(np.min(v2))
    assert np.abs(np.min(v2[v2 > tol]) - 1.0) < tol, 'min(v1[v2 > 0]): {}'.format(np.min(v2[v2 > 0]))
    # TODO: Is this rescaling enough?

    if np.sum((v1 - v2) ** 2) < tol: return True
    # try to see if they have the same tally
    u1, c1 = np.unique(np.round(v1, decimals=3), return_counts=True)
    u2, c2 = np.unique(np.round(v2, decimals=3), return_counts=True)
    if not u1.shape[0] == u2.shape[0]: return False
    if not np.all(u1 == u2): return False
    if not np.all(c1 == c2): return False

    # check if any relabelling is the same -> we have to recalculate v2 but not do the tally check again as its just relabelled
    for relabel in relabels:
        bell_tmp = bell2[relabel]
        v2 = dets @ bell_tmp
        if np.sum((v1 - v2) ** 2) < tol ** 2: return True
    # the two expressions were not equivalent -> so they are from different classes
    return False


def check_equiv_bell_vertex_enum_non_rescale(bell1, bell2, relabels, dets, tol=1e-6):
    """
    Checks if two bell inequalities are equivalent under relabelling and for each perm checks that if they are equiv
    This function uses the relabels of the bell expressions and therefore needs less memory, but has to calculate the
    v vectors multiple times.
    :param bell1: first bell expression
    :param bell2: second bell expression
    :param relabels: permutations of the bell expression that are allowed -> from get_allowed_relabellings
    :param dets: deterministic behaviors
    :param tol: tolerance
    :return:
    """
    # check if they are the same
    if np.sum((bell1 - bell2) ** 2) < tol ** 2:
        return True
    # get the v vectors
    v1 = dets @ bell1
    v2 = dets @ bell2
    # shift to min of v_i = 0
    v1 = v1 - np.min(v1)
    v2 = v2 - np.min(v2)
    # rescale that second min is = 1
    v1 = v1 / np.min(v1[v1 > tol])
    v2 = v2 / np.min(v2[v2 > tol])
    # check that smallest value is zero, due to affine_transformation
    assert np.abs(np.min(v1)) < tol, 'min(v1): {}'.format(np.min(v1))
    assert np.abs(np.min(v1[v1 > tol]) - 1.0) < tol, 'min(v1[v1 > 0]): {}'.format(np.min(v1[v1 > 0]))
    assert np.abs(np.min(v2)) < tol, 'min(v2): {}'.format(np.min(v2))
    assert np.abs(np.min(v2[v2 > tol]) - 1.0) < tol, 'min(v1[v2 > 0]): {}'.format(np.min(v2[v2 > 0]))
    if np.sum((v1 - v2) ** 2) < tol: return True
    # try to see if they have the same tally
    u1, c1 = np.unique(np.round(v1, decimals=3), return_counts=True)
    u2, c2 = np.unique(np.round(v2, decimals=3), return_counts=True)
    if not u1.shape[0] == u2.shape[0]: return False
    if not np.all(u1 == u2): return False
    if not np.all(c1 == c2): return False

    # rescale bell expressions if relabels are given
    bell1, bell2 = affine_transform_bell([bell1, bell2], dets)

    # check if any relabelling is the same -> we have to recalculate v2 but not do the tally check again as its just relabelled
    for relabel in relabels:
        bell_tmp = bell2[relabel]
        if np.all(bell_tmp == bell2):
            continue
        if check_equiv_bell_vertex_enum(bell1, bell_tmp, [], dets):
            return True
    # the two expressions were not equivalent -> so they are from different classes
    return False


def get_relabels_dets(dets, allowed_perms, show_progress=0):
    """ Gets a list of which deterministic transforms to which under each relabelling """
    # list for relabels
    relabels_dets = []
    # iterate through possible relabels
    for j, perm in enumerate(allowed_perms):
        if show_progress:
            print('perm: {} / {}'.format(j, len(allowed_perms)))
        tmp_relabel_d = np.zeros(dets.shape[0])
        for i in range(dets.shape[0]):
            d = dets[i]
            idx = np.argmin(np.sum(np.abs(dets - d[perm]), axis=1))
            tmp_relabel_d[i] = idx
        relabels_dets.append(tmp_relabel_d)
    return np.array(relabels_dets, dtype=int)


def parametrise_behavior(p, configs, configs_param, inputs_a, inputs_b, outputs_a, outputs_b):
    """
    Parametrises a given behavior to the lower dimensional representation:
        { p(a|x), p(b|y), p(ab|xy) } for all x,y and a,b = 1,...,n-1
    """
    assert len(inputs_a) == len(inputs_b), 'Can only parametrise for equal number of inputs'
    assert len(outputs_a) == len(outputs_b), 'Can only parametrise for equal number of outputs'
    m = len(inputs_a)
    d = len(outputs_a)
    # create parametrised array of lower dimensions
    t = 2 * (d - 1) * m + (d - 1) * (d - 1) * m * m
    param = np.zeros(t)
    assert len(configs_param) == t

    for i, (a, b, x, y) in enumerate(configs_param):
        if b == -1:
            # sum all probabilities together
            s = 0
            for b_tmp in outputs_b:
                idx = configs.index((a, b_tmp, x, y))
                s += p[idx]
            # set parametrised probability to marginal probability
            param[i] += s
        elif a == -1:
            # same as for b == -1
            s = 0
            for a_tmp in outputs_a:
                idx = configs.index((a_tmp, b, x, y))
                s += p[idx]
            param[i] += s
        else:
            # copy the probability
            idx = configs.index((a, b, x, y))
            param[i] += p[idx]
    return param


def get_parametrisation_configs(inputs_a, inputs_b, outputs_a, outputs_b):
    """
    Gets the configurations for the parametrised representation of a behavior.
    We use "-1" as an identifier that this configuration is not used.
    Due to the no signalling constraint we can choose any y or x for the marginals. So we always choose the first
    """
    configs = []
    # slice the last part of the outputs
    outputs_a = np.copy(outputs_a[:-1])
    outputs_b = np.copy(outputs_b[:-1])
    # first part is Alice side
    for a, x in product(outputs_a, inputs_a):
        c = (a, -1, x, inputs_b[0])
        configs.append(c)
    # second part is Bobs side
    for b, y in product(outputs_b, inputs_b):
        c = (-1, b, inputs_a[0], y)
        configs.append(c)
    # third part is usual configs with one less output for each party
    c = get_configs(inputs_a, inputs_b, outputs_a, outputs_b)
    # return all configs
    return configs + c


def deperametrise_behavior(p_param, configs, configs_param, inputs_a, inputs_b, outputs_a, outputs_b):
    """ Deparametrise a behavior, that was parametrised by parametrise behasvior """
    assert len(inputs_a) == len(inputs_b), 'Can only parametrise for equal number of inputs'
    assert len(outputs_a) == len(outputs_b), 'Can only Parametrise for equal number of outputs'
    m = len(inputs_a)
    d = len(outputs_a)
    # define the value of the last output, which is dropped in the parametrisation
    last_out = outputs_a[-1]
    # create empty behahior
    p = np.zeros((d ** 2) * (m ** 2))
    # iterate through the configurations
    for i, (a, b, x, y) in enumerate(configs):
        # if this particular configuration is already in the parametrisation we just copy it
        if (a, b, x, y) in configs_param:
            idx = configs_param.index((a, b, x, y))
            p[i] = p_param[idx]
        # check the different cases and perform the sums
        if a == last_out and b != last_out:
            # TODO: chaage this setting of x = inputs_a[0] to just set it to -1 as it does not matter which x it is (must be done in parametrise)
            idx = configs_param.index((-1, b, inputs_a[0], y))
            p[i] += p_param[idx]
            for a_tmp in outputs_a[:-1]:
                idx = configs_param.index((a_tmp, b, x, y))
                p[i] -= p_param[idx]

        if a != last_out and b == last_out:
            # TODO: Same as in if-clause before
            idx = configs_param.index((a, -1, x, inputs_b[0]))
            p[i] += p_param[idx]
            for b_tmp in outputs_b[:-1]:
                idx = configs_param.index((a, b_tmp, x, y))
                p[i] -= p_param[idx]
    # set the probabilities for p(dd|xy) as last, as we can already use the other ones directly
    for i, (a, b, x, y) in enumerate(configs):
        if a == last_out and b == last_out:
            # value for p(dd|xy)
            val = 1
            # sum over all possible output combinations (the p(dd|xy) is still zero, so just include it)
            for a_tmp in outputs_a:
                for b_tmp in outputs_b:
                    idx = configs.index((a_tmp, b_tmp, x, y))
                    val -= p[idx]
            # set p(dd|xy)
            p[i] = val
    return p


def deparametrise_bell_expression(b_param, configs, configs_param, inputs_a, inputs_b, outputs_a, outputs_b):
    """ Deparametrisation of a bell expression that was found using parametrised behaviors """
    assert len(inputs_a) == len(inputs_b), 'Can only parametrise for equal number of inputs'
    assert len(outputs_a) == len(outputs_b), 'Can only Parametrise for equal number of outputs'
    m = len(inputs_a)
    d = len(outputs_a)
    # define the last output
    last_out = outputs_a[-1]
    # setup bell expression
    bell = np.zeros((m ** 2) * (d ** 2))
    # iterate through the configurations
    for i, (a, b, x, y) in enumerate(configs):
        # check if the configuration is in the parametrised ones
        if (a, b, x, y) in configs_param:
            # add the corresponding values to this bell expression
            idx = configs_param.index((a, b, x, y))
            bell[i] += b_param[idx]
            if x == inputs_a[0]:
                idx = configs_param.index((-1, b, x, y))
                bell[i] += b_param[idx]
            if y == inputs_b[0]:
                idx = configs_param.index((a, -1, x, y))
                bell[i] += b_param[idx]
        # if a == delta and b not
        if a == last_out and b != last_out:
            if x == inputs_a[0]:
                idx = configs_param.index((-1, b, x, y))
                bell[i] += b_param[idx]
        if a != last_out and b == last_out:
            if y == inputs_b[0]:
                idx = configs_param.index((a, -1, x, y))
                bell[i] += b_param[idx]

    # return the bell expression
    return bell


def affine_transform_bell(bell_expressions, dets):
    """ Shifts the bell expressions such that min(b @ det) = 0 and second_min(b @ det) = 1 """
    shifted_bell = []
    sum_ones = np.sum(dets[0])
    for b in bell_expressions:
        v = np.array([b @ d for d in dets])
        min_val = np.min(v)
        # shift v that lowest value is 0
        v = v - min_val
        # get second lowest value
        sec_min_pos_val = np.min(v[v > 1e-4])
        # shift the bell expression
        bell = (b - min_val / sum_ones) / sec_min_pos_val
        shifted_bell.append(bell)
    return np.array(shifted_bell)


def check_equiv_bell_vertex_enum_fast(bell1, bell2, relabels, dets, tol=1e-6):
    """
    Checks if two bell inequalities are equivalent under relabelling and for each perm checks that if they are equiv
    This function uses the relabels of the bell expressions and therefore needs less memory, but has to calculate the
    v vectors multiple times.
    We assume that the bell expressions are already rescaled by the function affine_transform_bell
    :param bell1: first bell expression
    :param bell2: second bell expression
    :param relabels: permutations of the bell expression that are allowed -> from get_allowed_relabellings
    :param dets: deterministic behaviors
    :param tol: tolerance
    :return:
    """
    # check if they are the same
    if np.sum((bell1 - bell2) ** 2) < tol ** 2:
        return True
    # get the v vectors
    v1 = dets @ bell1
    v2 = dets @ bell2

    if np.sum((v1 - v2) ** 2) < tol: return True
    # try to see if they have the same tally
    u1, c1 = np.unique(np.round(v1, decimals=3), return_counts=True)
    u2, c2 = np.unique(np.round(v2, decimals=3), return_counts=True)
    if not u1.shape[0] == u2.shape[0]: return False
    if not np.all(u1 == u2): return False
    if not np.all(c1 == c2): return False

    # check if any relabelling is the same -> we have to recalculate v2 but not do the tally check again as its just relabelled
    bell2_relabels = bell2[relabels]
    return np.any(np.sum((bell1 - bell2[relabels]) ** 2, axis=1) < tol ** 2)


def relabellings_parametrised(configs_param, inputs_a, inputs_b, outputs_a, outputs_b):
    """ Calculates the relabellings for parametrised setup """
    relabels = get_allowed_relabellings(inputs_a, inputs_b, outputs_a[:-1], outputs_b[:-1])
    configs = get_configs(inputs_a, inputs_b, outputs_a, outputs_b)
    param_relabels = []
    for relabel in relabels:
        # current parametrised relabelling
        p_relabel = []
        # indicator for valid relabel
        validRelabel = True
        for c in configs_param:

            a, b, x, y = c
            # check if current config is a marginal
            if a == -1:
                idx = configs.index((0, b, x, y))
                ridx = relabel[idx]
                _, b_new, x_new, y_new = configs[ridx]
                a_new = a
            elif b == -1:
                idx = configs.index((a, 0, x, y))
                ridx = relabel[idx]
                a_new, _, x_new, y_new = configs[ridx]
                b_new = b
            else:
                idx = configs.index((a, b, x, y))
                ridx = relabel[idx]
                a_new, b_new, x_new, y_new = configs[ridx]
            # Try to find the index (continue if not found as some relabellings do not matter, e.g. party swap)
            try:
                idx = configs_param.index((a_new, b_new, x_new, y_new))
                p_relabel.append(idx)
            except ValueError:
                validRelabel = False
        if validRelabel:
            assert len(p_relabel) == len(configs_param)
            param_relabels.append(p_relabel)
    param_relabels = np.array(param_relabels)
    return np.unique(param_relabels, axis=1)


def equiv_check_adjacency_testing(bell1, bell2, relabels, dets, tol=1e-6):
    """ Equivalence check testing for the partial adjacency """
    # get the v vectors
    v1 = dets @ bell1
    v2 = dets @ bell2

    # shift to min of v_i = 0
    v1 = v1 - np.min(v1)
    v2 = v2 - np.min(v2)

    # rescale that second min is = 1
    v1 = v1 / np.min(v1[v1 > tol])
    v2 = v2 / np.min(v2[v2 > tol])

    # round the vectors for better equality checks
    v1 = np.round(v1, 5)
    v2 = np.round(v2, 5)

    if np.all(v1 == v2): return True
    # try to see if they have the same tally
    t1 = np.bincount(v1.astype(int))
    t2 = np.bincount(v2.astype(int))
    if not np.all(t1 == t2):
        return False

    # check if any relabelling is the same -> we have to recalculate v2 but not do the tally check again as its just relabelled
    for i in range(relabels.shape[0]):
        bell_tmp = bell2[relabels[i]]
        if np.any(bell_tmp != bell2):
            if equiv_check_adjacency_testing_v(v1, bell_tmp, dets, tol=tol):
                return True
    # the two expressions were not equivalent -> so they are from different classes
    return False

def equiv_check_adjacency_testing_v(v1, bell2, dets, tol=1e-6):
    """ Performs equivalence checking, but already has the first part calculated """
    v2 = dets @ bell2
    v2 = v2 - np.min(v2)
    try:
        v2 = v2 / np.min(v2[v2 > tol])
    except ValueError as e:
        pass
    v2 = np.round(v2, 5)

    if np.all(v1 == v2): return True
    return False
