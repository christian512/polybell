import numpy as np
import argparse
fac = np.math.factorial

parser = argparse.ArgumentParser()
parser.add_argument(dest='ma', help="number of inputs for ALICE")
parser.add_argument(dest='mb', help='number of inputs for BOB')
parser.add_argument(dest='na', help='number of outputs for ALICE')
parser.add_argument(dest='nb', help='number of outputs for BOB')
args = parser.parse_args()
ma, mb, na, nb = int(args.ma), int(args.mb), int(args.na), int(args.nb)

dim = ma*(na-1)*mb*(nb-1) + ma * (na-1) + mb*(nb-1)
vert = (na**ma)*(nb**mb)
sym = fac(ma) * fac(mb) * (fac(na)**ma) * (fac(nb)**mb)
print('dimensions: ', dim)
print('vert: ', vert)
print('sym: ', sym)