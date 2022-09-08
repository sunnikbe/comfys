import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('p2.txt')

x = data[:, 0] # x - vector
u_x = data[:, 1] # u(x) values

plt.plot(x, u_x)
plt.title(r"$u(x)=1-(1-e^{-10})x-e^{-10x}$")
plt.xlabel(r"$x\in[0, 1]$")
plt.ylabel(r"$u(x)$")
plt.savefig("p2.pdf")
