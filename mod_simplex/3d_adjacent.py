import numpy as np

lhs = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0], [0, -1, 0], [0, 0, -1]])
rhs = np.array([1, 1, 1, 1, 1, 1])

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
    tab[row, :] = tab[row, :] / pivot_elem
    assert tab[row, col] == 1.0
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


# set a start point
init_vertex = 1.0 * np.ones(lhs.shape[1])
# check which inequalities are equalised by this one
eq_indices = equal_inequalities_indices(init_vertex, lhs, rhs)
assert len(eq_indices) == init_vertex.shape[0]

# pivot the tableau around the equalising inequalities
# TODO: This is where it get more difficult if the number of equalities is higher than number of dimensions
pivots = pivots_make_x_basic(tableau, nonbasic_vars, basic_vars, eq_indices)
while pivots:
    row, col = pivots[0]
    tableau, nonbasic_vars, basic_vars = pivot(tableau, row, col, nonbasic_vars, basic_vars)
    pivots = pivots_make_x_basic(tableau, nonbasic_vars, basic_vars, eq_indices)
# read the current vertex from tableau, it should be the initial vertex
vertex = read_vertex(tableau, basic_vars)
assert np.all(init_vertex == vertex)

# Find all pivots that change an equalized with a non equalized row
pivots = pivots_non_equalities(tableau, nonbasic_vars, basic_vars, eq_indices)
for p in pivots:
    row, col = p
    new_tab, new_nonbasic, new_basic = pivot(tableau, row, col, nonbasic_vars, basic_vars)
    acceptable = check_acceptable(new_tab, new_basic)
    print('acceptable: ', acceptable)
    if acceptable:
        v = read_vertex(new_tab, new_basic)
        print('vertex: ', v)
        break
