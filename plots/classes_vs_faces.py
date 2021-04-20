""" Here we will plot the number of classes vs the number of faces for known Bell Polytopes """

import matplotlib.pyplot as plt

nfaces = [24, 1116, 684, 1260, 12480, 71340, 36391264, 252558]
nclasses = [2, 3, 3, 5, 6, 7, 175, 25]

plt.scatter(nfaces, nclasses)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Number of Faces')
plt.ylabel('Number of Classes')
plt.title('Bell Polytope Sizes')
plt.savefig('bell_polytope_sizes.png',dpi=300)

