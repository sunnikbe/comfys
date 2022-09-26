import numpy as np
import matplotlib.pyplot as plt
import heapq as hq

v_i_vectors = np.loadtxt("eigenvectors.txt")
lambda_ = np.loadtxt("eigenvalues.txt")

# making list with the 3 smallest eigenvalues
lambdas = hq.nsmallest(3, lambda_)

# finding indices of smallest eigenvalues to get corresponding eigenvectors
v_i = []
for i in range(len(lambdas)):
    index = np.where(lambda_ == lambdas[i])
    v_i.append(v_i_vectors[:, index])


# x vector:
x_dimless = np.linspace(0, 1, len(v_i[0]) + 2)

# Boundary points:
v_0 = 0
v_n = 0


# Plotting
labels = [fr'$\lambda_1$ = {lambdas[0]:.1f}', f'$\lambda_2$ = {lambdas[1]:.1f}', f'$\lambda_3$ = {lambdas[2]:.1f}']


# This part shows how the solutions would have different lengths
def L(eigenvalue):
    # random numbers just for scaling
    F = 10
    gamma = 100
    return np.sqrt((eigenvalue*gamma)/F)


for i in range(len(v_i)):
    v_i[i] = np.insert(v_i[i], 0, v_0)
    v_i[i] = np.append(v_i[i], v_n)

    x = x_dimless*L(lambdas[i]) # plot this to see dim

    plt.plot(x_dimless, v_i[i], label = labels[i])

n = len(lambda_) + 1
plt.title(fr'Analytical eigenvectors $v_i$ vs. $\hatx_i$ for n = {n} steps')
plt.xlabel(r'x $\in$ [0, 1]')
plt.ylabel(r'$v_i$')
plt.legend()
plt.savefig(f'buckling_beams_n{n}.pdf')
