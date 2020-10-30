import gurobipy as gp
from gurobipy import GRB
import numpy as np
from scipy import sparse

env = gp.Env(empty=True)
env.setParam('OutputFlag', 0)
env.start()


def find_local_weight_guro(p, dets):
    # create a model
    m = gp.Model('local_weight', env=env)
    # add variable
    bell = m.addMVar(shape=p.shape[0], vtype=GRB.SEMICONT, name='bell')
    # set objective
    m.setObjective(bell @ p, GRB.MINIMIZE)
    # setup variables for constraints
    rhs = np.ones(dets.shape[0])
    A = sparse.csr_matrix(dets)
    # add constraint
    m.addConstr(A @ bell >= rhs, name='c')
    # optimize
    m.optimize()
    # return
    return bell.X
