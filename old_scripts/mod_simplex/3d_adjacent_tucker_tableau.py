""" Testing a Tucker tableau approach ( tableua without cost function) and different pivot """

import numpy as np

# Linear inequality system for a cube
lhs = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
rhs = np.array([1, 1, 1])
# 0.5 x2 <= x1 is equiv to -x1 + 0.5 x2 <= 0
lhs = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 1 / 2, 0]])
rhs = np.array([1, 1, 1, 0])
number_of_vars = lhs.shape[0]


def create_tucker_tableau(lhs, rhs):
    """ Creates a Tucker tableau from a set of inequalities definied by lhs * x <= rhs """
    # TODO: Is it always that simple?
    return np.c_[lhs, rhs]


def check_feasible_tab(tab, basic_x):
    """ Checks if a tableau is feasible, i.e. if it's a vertex """
    assert tab.shape[0] == basic_x.shape[0]
    for i in range(basic_x.shape[0]):
        if basic_x[i] == -1:
            if tab[i, -1] < 0:
                return False
    return True


def pivot(tab, row, col, nonbasic_x, basic_x):
    """ Pivots the Tucker tableau """
    # Some assertions
    assert col <= tab.shape[1] - 1, 'Can not pivot around last column'
    assert nonbasic_x[col] != -1, 'Pivot would not make x_i basic'
    assert basic_x[row] == -1, 'Pivot would make x_i nonbasic again'
    assert tab[row, col] != 0.0, 'Pivot element is zero! Can not pivot.'
    # get pivot element and transform row
    pivot_elem = tab[row, col]
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
    basic_x[row] = nonbasic_x[col]
    nonbasic_x[col] = -1
    return tab, nonbasic_x, basic_x


def get_possible_pivots(tab, nonbasic_x, basic_x):
    """ Returns a list of all possible pivot points"""
    assert nonbasic_x.shape[0] == tab.shape[1] - 1
    assert basic_x.shape[0] == tab.shape[0]
    poss_pivots = []
    for i in range(basic_x.shape[0]):
        if basic_x[i] != -1:
            continue
        for j in range(nonbasic_x.shape[0]):
            if nonbasic_x[j] == -1:
                continue
            if tab[i, j] != 0.0:
                poss_pivots.append([i, j])
    return poss_pivots


def pivots_on_equalized(tab, nonbasic_x, basic_x, eq_rows):
    """ Finds pivots around rows that are equalized by initial point """
    assert nonbasic_x.shape[0] == tab.shape[1] - 1
    assert basic_x.shape[0] == tab.shape[0]
    poss_pivots = []
    for i in range(basic_x.shape[0]):
        if basic_x[i] != -1:
            continue
        if i not in eq_rows:
            continue
        for j in range(nonbasic_x.shape[0]):
            if nonbasic_x[j] == -1:
                continue
            if tab[i, j] != 0.0:
                poss_pivots.append([i, j])
    return poss_pivots


def pivots_not_on_equalized(tab, nonbasic_x, basic_x, eq_rows):
    """ Finds pivots around rows that are equalized by initial point """
    assert nonbasic_x.shape[0] == tab.shape[1] - 1
    assert basic_x.shape[0] == tab.shape[0]
    poss_pivots = []
    for i in range(basic_x.shape[0]):
        if basic_x[i] != -1:
            continue
        if i in eq_rows:
            continue
        for j in range(nonbasic_x.shape[0]):
            if nonbasic_x[j] == -1:
                continue
            if tab[i, j] != 0.0:
                poss_pivots.append([i, j])
    return poss_pivots


def read_vertex(tab, basic_x):
    """ Reads the vertex from a tableau """
    assert tab.shape[0] == basic_x.shape[0]
    # Non basic variables are set to zero
    vertex = np.zeros(tab.shape[1] - 1)
    for i in range(basic_x.shape[0]):
        idx = int(basic_x[i])
        if idx != -1:
            # set vertex value to continue
            vertex[idx] = tab[i, -1]
    return vertex


def find_adjacent_vertices(tab, nonbasic_x, basic_x, poss_pivots, eq_rows, recursive=0):
    """ Finds all adjacent vertices to current vertex in tab """
    adj_vertices = []
    # copy the start value
    start_basic_x = np.copy(basic_x)
    start_nonbasic_x = np.copy(nonbasic_x)
    start_tab = np.copy(tab)

    for p in poss_pivots:
        # load initial tableau and start basic and nonbasic variables
        print('pivot: ', p)
        basic_x = np.copy(start_basic_x)
        nonbasic_x = np.copy(start_nonbasic_x)
        tab = np.copy(start_tab)
        # load row and column from pivot
        row, col = p
        # perform transformation
        tab, nonbasic_x, basic_x = pivot(tab, row, col, nonbasic_x, basic_x)
        # check if feasible tableau
        if check_feasible_tab(tab, basic_x):
            vertex = read_vertex(tab, basic_x)
            adj_vertices.append(vertex)
            print('recursive level: ', recursive)
            print('vertex: ', vertex)
        # check if there are some possible pivots
        # TODO: Do I only need to do it with else here
        else:
            sub_poss_pivots = pivots_on_equalized(tab, nonbasic_x, basic_x, eq_rows)
            print('sub poss pivots: ', sub_poss_pivots)
            if sub_poss_pivots:
                rec_adj_vertices = find_adjacent_vertices(tab, nonbasic_x, basic_x, sub_poss_pivots, eq_rows,
                                                          recursive=recursive + 1)
                for v in rec_adj_vertices:
                    adj_vertices.append(v)
    return adj_vertices


# create tableau
tableau = create_tucker_tableau(lhs, rhs)
init_tableau = np.copy(tableau)

# initial vertex
init_vertex = np.zeros(tableau.shape[1] - 1)
# get the inequalities that are equalized by this
eq_rows = []
for i in range(lhs.shape[0]):
    if lhs[i] @ init_vertex == rhs[i]:
        eq_rows.append(i)
print('rows equalized by initial point: ', eq_rows)

# storage which row contains which x_i variable
basic_x = -1.0 * np.ones(tableau.shape[0])
nonbasic_x = np.arange(tableau.shape[1] - 1)
init_basic_x = np.copy(basic_x)
init_nonbasic_x = np.copy(nonbasic_x)

# get pivots that don't act on rows that are equalized by initial vertex
non_eq_pivots = pivots_not_on_equalized(tableau, nonbasic_x, basic_x, eq_rows)
print('pivots not on equalized: ', non_eq_pivots)

adj_vertices = find_adjacent_vertices(tableau, nonbasic_x, basic_x, non_eq_pivots, eq_rows)
print(adj_vertices)
