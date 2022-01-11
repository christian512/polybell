import numpy as np
fac = np.math.factorial

ma = 7
mb = 3
na = 2
nb = 2


dim = ma*(na-1)*mb*(nb-1) + ma * (na-1) + mb*(nb-1)
vert = (na**ma)*(nb**mb)
sym = fac(ma) * fac(mb) * (fac(na)**ma) * (fac(nb)**mb)
print('dimensions: ', dim)
print('vert: ', vert)
print('sym: ', sym)