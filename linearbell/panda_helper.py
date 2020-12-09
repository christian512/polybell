""" Here I collect helper functions for the panda vertex enumeration code """
import numpy as np
from fractions import Fraction
import subprocess


def write_panda_input_inequalities(lhs, rhs, idx_equalities=[], symmetries=[[]], file='', denom_limit=10000):
    """
    the inequality is interpreted as lhs * x <= rhs
    lhs: left hand side of inequalities
    rhs: righ hand side of inequalities
    symmetries: the symmetries that we want to consider given in a list of lists
                each list tells us how the variables are interchanged.
    """
    # check dimensions
    if type(symmetries) == list:
        symmetries = np.array(symmetries)
    assert lhs.shape[0] == rhs.shape[0]
    # string to write out
    string = ""
    # variable names (only needed for symmetries)
    string += 'Names:\n'
    for i in range(lhs.shape[1]):
        string += 'a' + str(i) + ' '
    string += '\n'
    # add symmetries
    if len(symmetries) > 0:
        string += 'Maps:\n'
        for symm in symmetries:
            for i in range(symm.shape[0]):
                string += 'a' + str(symm[i]) + ' '
            string += '\n'
    # equalities
    if idx_equalities:
        string += 'Equations:\n'
        for i in idx_equalities:
            for x in lhs[i]:
                string += str(Fraction(x).limit_denominator(denom_limit)) + ' '
            # TODO: check if this = is correct or can I just drop it?
            string += '= ' + str(Fraction(rhs[i]).limit_denominator(denom_limit)) + '\n'
    # inequalities -> as we can not give fractions we multiply the inequality by the product of unique denominators
    string += 'Inequalities:\n'
    for i in range(lhs.shape[0]):
        if i in idx_equalities:
            continue
        # find factor to multiply inequality with, such that we have a only integers inequality
        arr = [Fraction(x).limit_denominator(denom_limit).denominator for x in lhs[i]]
        arr.append(Fraction(rhs[i]).denominator)
        factor = np.product(np.unique(arr))
        # write the inequalitiy
        for j in range(lhs.shape[1]):
            val = int(lhs[i, j] * factor)
            string += str(val) + 'a' + str(j) + ' '
        string += '<= ' + str(int(rhs[i] * factor)) + '\n'
    # write file
    if file:
        f = open(file, 'w+')
        f.write(string)
        f.close()
    return string


def run_panda(file, outfile=''):
    """ Runs panda on a file"""
    cmd = 'panda ' + file
    if outfile:
        cmd += ' > ' + outfile
    out = subprocess.run(cmd, shell=True)
    return True
