from fractions import Fraction
import numpy as np
import subprocess

def polyhedra_h_representation(lhs_ineq, rhs_ineq, name='polytope', file=''):
    """
    This function creates a file to input into lrs. The input used here is lhs_ineq * x <= rhs_ineq
    Parameters
    ----------
    lhs_ineq : left hand side of the inequality
    rhs_ineq : right hand side of the inequality
    name : name of the polytope
    file : name of the file

    Returns
    -------

    """
    # this is the string that will be written to the file in the end
    string = ''
    # setup the header
    string += name + '\n'
    string += 'H-representation \n'
    string += 'begin \n'
    # check that the inequalities have same shape
    assert lhs_ineq.shape[0] == rhs_ineq.shape[0]
    # rhs must be a one dimensional array
    assert len(rhs_ineq.shape) == 1
    string += str(rhs_ineq.shape[0]) + ' ' + str(lhs_ineq.shape[1] + 1) + ' rational \n'
    # now lrs wants input a0 + a1 * x1 + ... + an-1 * xn-1 >= 0
    # set we have to reformulate our input
    # lhs * x <= rhs is equiv to rhs - lhs * x >= 0
    lhs_ineq = -1.0 * lhs_ineq
    for i in range(lhs_ineq.shape[0]):
        b = Fraction(rhs_ineq[i])
        string += str(b)
        for a in lhs_ineq[i]:
           string += ' ' + str(Fraction(a))
        string += '\n'
    string += 'end \n'
    # write string to file
    if file:
        f = open(file, 'w+')
        f.write(string)
        f.close()
        print('wrote file: {}'.format(file))
    return string

def run_lrs_polytope(input_file, output_file='out.ext'):
    """
    Runs lrs with a given input file and stores the result and gives the vertices back.
    This will additionally check the output if it's a polytope
    Parameters
    ----------
    input_file : input file
    output_file : output file

    Returns: vertices
    -------
    """
    # TODO: Read out number of rays and vertices to check if it's a polytope (polytope -> no rays)
    cmd = 'lrs ' + input_file + ' ' + output_file
    out = subprocess.run(cmd, shell=True)
    # read output file
    f = open(output_file, 'r')
    string = f.read()
    # get the part of the output file that contains the vertices
    string = string.split('begin')[1]
    string = string.split('end')[0]
    string = string.split('rational')[1]
    # get a string list of the vertices
    vertices_string = string.split('\n')[1:-1]
    # array for the vertices
    vertices = []
    # iterate over the strings of vertices
    for v in vertices_string:
        # array for teh vertex
        vertex = []
        # remove the newline characters
        v.replace('\n', '')
        # iterate through the split
        for value in v.split()[1:]:
            # append the value to the vertex
            vertex.append(float(Fraction(value)))
        # append the vertex to the vertices
        vertices.append(vertex)
    # return the vertices as a numpy array
    return np.array(vertices)
