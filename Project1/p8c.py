import numpy as np
import matplotlib.pyplot as plt

# c) Plotting the maximum relative error

data = np.loadtxt("n_steps_max_rel_err.txt")
n_steps = data[:, 0]
max_rel_err = data[:, 1]

plt.plot(n_steps, np.log(max_rel_err), marker = "o")
plt.title("Max relative error (logarithmic scale)")
plt.xlabel(r"$n_{steps}$")
plt.ylabel(r"$log_{10}(Max(\epsilon_i))$")
plt.xscale("log")
plt.yscale("log")
plt.savefig('p8c.pdf')
