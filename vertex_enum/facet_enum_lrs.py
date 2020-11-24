import argparse
from linearbell.utils import get_deterministic_behaviors, get_configs, general_pr_box_extended, \
    get_parametrisation_configs, parametrise_behavior, check_equiv_bell
from lrs_helper import polytope_v_representation
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument(dest='ma', help="number of inputs for ALICE")
parser.add_argument(dest='mb', help='number of inputs for BOB')
parser.add_argument(dest='n', help='number of outputs')
args = parser.parse_args()
ma = int(args.ma)
mb = int(args.mb)
n = int(args.n)

# set inputs / outputs
inputs_a = range(ma)
inputs_b = range(mb)
outputs = range(n)

# get deterministic points
dets = get_deterministic_behaviors(inputs_a, inputs_b, outputs)

# get configurations
configs = get_configs(inputs_a, inputs_b, outputs, outputs)

vrepr = polytope_v_representation(dets, file='input_v.ine')
#print configs
for i, c in enumerate(configs):
    a,b,x,y = c
    print('i: {} ||| a: {} | b: {} | x: {} | y: {}'.format(i,a,b,x,y))
