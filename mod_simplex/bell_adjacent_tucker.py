import numpy as np
from linearbell.utils import get_deterministic_behaviors, parametrise_behavior, get_configs, \
    get_parametrisation_configs, deparametrise_bell_expression, equiv_check_adjacency_testing
from itertools import combinations
from scipy.special import binom

# get deterministic behaviors
m = 2
n = 2
inputs = list(range(m))
outputs = list(range(n))
dets = get_deterministic_behaviors(inputs, inputs, outputs)
configs = get_configs(inputs, inputs, outputs, outputs)
configs_param = get_parametrisation_configs(inputs, inputs, outputs, outputs)

# define initial point as starting vertex
init_vertex = np.zeros(dets.shape[1])
# get deterministic points that don't share ones in same column
unrelated_dets = [dets[0]]
for d in dets:
    related = False
    for ud in unrelated_dets:
        if np.any(d + ud > np.max(ud)):
            related = True
            break
    if not related:
        unrelated_dets.append(d)
print('number of unrelated dets: ', len(unrelated_dets))
for i in range(len(unrelated_dets) - 1):
    idx = np.argwhere(unrelated_dets[i] > 0.0).flatten()[0]
    print(idx)
    init_vertex[int(idx)] = len(outputs) ** 2

# shift the origin to make vertex description possible
p_origin = np.sum(dets, axis=0) / dets.shape[0]
dets_unshifted = np.copy(dets)
dets = dets - p_origin
# Parametrise
dets_param = [parametrise_behavior(d, configs, configs_param, inputs, inputs, outputs, outputs) for d in dets]
dets_param = np.array(dets_param)

# load relabeling
relabels = np.loadtxt('../data/relabels/{}{}{}{}.gz'.format(m, m, n, n)).astype(int)

# set inequality system
lhs = np.copy(dets_param)
rhs = np.ones(lhs.shape[0])

# build the tableau
tableau = np.c_[lhs, rhs]
print('+++ Initial tableau +++')
print(tableau)
# setup basic and non basic variables
basic_vars = ['y' + str(i) for i in range(tableau.shape[0])]
nonbasic_vars = ['x' + str(i) for i in range(lhs.shape[1])]


def read_vertex(tab, basic):
    """ Reads a vertex from a tableau """
    vertex = np.zeros(tab.shape[1] - 1)
    for i in range(tab.shape[0]):
        if 'x' in basic[i]:
            idx = int(basic[i][1:])
            vertex[idx] = tab[i, -1]
    return vertex


def check_acceptable(tab, basic):
    """ Checks if a tableau is acceptable """
    for i in range(tab.shape[0]):
        if 'y' in basic[i] and tab[i, -1] < 0.0:
            return False
    return True


def pivot(ctab, row, col, nonbasic_x, basic_x):
    """ Pivots the Tucker tableau """
    # copy the tableau to not modify the original one
    tab = np.copy(ctab)
    # Some assertions
    assert col <= tab.shape[1] - 1, 'Can not pivot around last column'
    assert row <= tab.shape[0], 'Invalid row'
    assert 'x' not in basic_x[row], 'Pivot would make {} nonbasic again'.format(basic_x[row])
    assert tab[row, col] != 0.0, 'Pivot element is zero! Can not pivot.'
    # get pivot element and transform row
    pivot_elem = tab[row, col]
    # set pivot element to one
    tab[row, col] = 1.0
    tab[row, :] = tab[row, :] / pivot_elem
    # now every row has to be transformed
    for i in range(tab.shape[0]):
        if i == row:
            continue
        row_piv_elem = tab[i, col]
        # Iterate over columns
        for j in range(tab.shape[1]):
            if j == col:
                tab[i, j] = -tab[i, j] / pivot_elem
            else:
                tab[i, j] = tab[i, j] - row_piv_elem * tab[row, j]
    # Update basic and nonbasic variables
    basic = np.copy(basic_x)
    nonbasic = np.copy(nonbasic_x)
    basic[row], nonbasic[col] = nonbasic_x[col], basic_x[row]
    return tab, nonbasic, basic


def equal_inequalities_indices(vertex, lhs, rhs):
    """ Returns indices of rows which are equalised by vertex """
    assert lhs.shape[0] == rhs.shape[0]
    assert lhs.shape[1] == vertex.shape[0]
    indices = []
    for i in range(lhs.shape[0]):
        if lhs[i] @ vertex == rhs[i]:
            indices.append(i)
    return indices


def pivots_make_x_basic(tab, nonbasic, basic, eq_idx):
    """ Gives back possible pivots to make x basic """
    pivots = []
    for row in range(tab.shape[0]):
        if 'x' in basic[row]:
            continue
        if row not in eq_idx:
            continue
        for col in range(tab.shape[1] - 1):
            if 'y' in nonbasic[col]:
                continue
            if tab[row, col] != 0.0:
                pivots.append([row, col])
    return pivots


def pivots_non_equalities(tab, nonbasic, basic, eq_idx):
    """ Gives back possible pivots around rows that are not equalized """
    pivots = []
    for row in range(tab.shape[0]):
        if row in eq_idx:
            continue
        assert 'y' in basic[row]
        for col in range(tab.shape[1] - 1):
            assert 'y' in nonbasic[col]
            if tab[row, col] != 0.0:
                pivots.append([row, col])
    return pivots

def find_possible_pivot_in_row(tab, nonbasic, basic, row):
    for j in range(tab.shape[1] - 1):
        # we only want to make x_i's basic
        if 'y' in nonbasic[j]:
            continue
        if tab[row, j] != 0.0:
            return [row, j]
    return None


# TODO: The initial vertex is not parametrised for now

# check which inequalities are equalised by the start point
eq_indices = equal_inequalities_indices(init_vertex, dets, rhs)
assert len(eq_indices) >= lhs.shape[1]

print('number of possible combinations of subspaces: ', binom(len(eq_indices), lhs.shape[1]))
# Create a list of which subspaces could be equalized at the same time
subspace_permutations = combinations(eq_indices, lhs.shape[1])
print('generated the permutations')
classes = []

for perm in subspace_permutations:
    tab = np.copy(tableau)
    nonbasic = np.copy(nonbasic_vars)
    basic = np.copy(basic_vars)
    # generate starting tableau
    for row in perm:
        valid_tab = True
        try:
            row, col = find_possible_pivot_in_row(tab, nonbasic, basic, row)
            tab, nonbasic, basic = pivot(tab, row, col, nonbasic, basic)
        except Exception as e:
            valid_tab = False
            break
    if not valid_tab:
        continue
    # read vertex from tableau
    vertex = read_vertex(tab, basic)
    # print('initial vertex in tableau: ', vertex)
    if len(classes) == 0:
        vertex = deparametrise_bell_expression(vertex, configs, configs_param, inputs, inputs, outputs, outputs)
        classes.append(vertex)

    # Find all pivots that change an equalized with a non equalized row
    pivots = pivots_non_equalities(tab, nonbasic, basic, eq_indices)
    for p in pivots:
        row, col = p
        new_tab, new_nonbasic, new_basic = pivot(tab, row, col, nonbasic, basic)
        acceptable = check_acceptable(new_tab, new_basic)
        if acceptable:
            v_param = read_vertex(new_tab, new_basic)
            v = deparametrise_bell_expression(v_param, configs, configs_param, inputs, inputs, outputs, outputs)
            equiv = False
            for c in classes:
                if equiv_check_adjacency_testing(c, v, relabels, dets_unshifted):
                    equiv = True
                    break
            if not equiv:
                classes.append(v)
                print('nonbasic: ', nonbasic)
                print('basic: ', basic)
                print('found new class: ', v)
