import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt("Nvsnum_iter.txt")
dense_dat = np.loadtxt("dense.txt")

# diag
N = data[:, 0]
num_iter = data[:, 1]

#dense
N_dense = dense_dat[:, 0]
num_iter_dense = dense_dat[:, 1]


N_poly = np.poly1d(np.polyfit(N, num_iter, 2))
N_poly_dense = np.poly1d(np.polyfit(N_dense, num_iter_dense, 2))


plt.plot(N_dense, N_poly_dense(N_dense), linestyle = 'dotted', color = 'orange', label = 'Fitted line for dense matrix')
plt.scatter(N_dense, num_iter_dense, marker = 'x', color = 'orange', label = 'Datapoints for dense matrix')

plt.plot(N, N_poly(N), color = 'royalblue', label = 'Fitted line for tridiagonal matrix')
plt.scatter(N, num_iter, marker = '.', color = 'blue', label = 'Datapoints for tridiagonal matrix')

plt.title('N vs. iterations for the Jacobi rotation method')
plt.xlabel('N')
plt.ylabel('Number of iterations')
plt.yticks(rotation = 60)
plt.legend()
plt.savefig('nvsiter.pdf')


